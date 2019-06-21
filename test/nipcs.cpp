#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "eqn.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpaceU_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR>;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR, 2>;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;

int main(int argc, char* argv[])
{
  MilliTimer t;

  ParameterDict config;

  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["nx"] = 4;
    config["ny"] = 8;
    config["dt"] = 0.1;
    config["ntime"]= 50U;
    config["nu"] = 0.1;
    config["printStep"] = 1U;
  }

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(
        *mesh,
        Vec3{0., 0., 0.},
        Vec3{1., 10., 0.},
        {config["nx"].as<uint>(), config["ny"].as<uint>(), 0});
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpaceP{*mesh};

  Eqn<FESpaceU_T, StorageType::RowMajor> eqnUstar{"uStar", feSpaceU};
  Eqn<FESpaceU_T, StorageType::RowMajor> eqnVstar{"vStar", feSpaceU};

  Eqn eqnP{"p", feSpaceP};

  Eqn eqnU{"u", feSpaceU};
  Eqn eqnV{"v", feSpaceU};
  t.stop();

  t.start("bc");
  auto const zeroS = [] (Vec3 const &) { return 0.; };
  auto const zeroV = [] (Vec3 const &) { return Vec2{0., 0.}; };
  // auto const inlet = [] (Vec3 const &) { return Vec2(0.0, 1.0); };
  auto const inlet = [] (Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p[0]*p[0])}; };
  auto const inlet0 = [&inlet] (Vec3 const & p) { return inlet(p)[0]; };
  auto const inlet1 = [&inlet] (Vec3 const & p) { return inlet(p)[1]; };
  // auto const pIn = [] (Vec3 const &) {return 3. * hy * nu;};

  BCList bcsVel{feSpaceVel};
  // last essential bc wins on corners
  bcsVel.addBC(BCEss{feSpaceVel, side::BOTTOM, inlet});
  bcsVel.addBC(BCEss{feSpaceVel, side::RIGHT, zeroV});
  bcsVel.addBC(BCEss{feSpaceVel, side::LEFT, zeroV, {0}});

  BCList bcsP{feSpaceP};
  bcsP.addBC(BCEss{feSpaceP, side::TOP, zeroS});
  // DofSet_T pinSet = {1};
  // bcsP.addBC(BCEss{feSpaceP, pinSet, zeroS});

  eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::BOTTOM, inlet0});
  eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::RIGHT, zeroS});
  eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::LEFT, zeroS});

  // last essential bc wins on corners
  eqnVstar.bcList.addBC(BCEss{eqnVstar.feSpace, side::BOTTOM, inlet1});
  eqnVstar.bcList.addBC(BCEss{eqnVstar.feSpace, side::RIGHT, zeroS});

  eqnP.bcList.addBC(BCEss{eqnP.feSpace, side::TOP, zeroS});
  // eqnP.bcList.addBC(BCEss{eqnP.feSpace, side::BOTTOM, pIn});
  t.stop();

  t.start("assembly");
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;

  auto const dt = config["dt"].as<double>();
  auto const nu = config["nu"].as<double>();

  Vec velOldMonolithic{2*dofU};
  AssemblyMass timeder(1./dt, feSpaceVel);
  AssemblyAdvection advection(1.0, velOldMonolithic, feSpaceVel, feSpaceVel);
  AssemblyTensorStiffness stiffness(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);
  AssemblyProjection timederRhs(1./dt, velOldMonolithic, feSpaceVel);
  // AssemblyBCNormal naturalBC{pIn, side::BOTTOM, feSpaceVel};
  // to apply bc on pressure
  // AssemblyMass dummy{0.0, feSpaceP, {0}, 2*dofU, 2*dofU};

  Builder<StorageType::RowMajor> builderM{dofU*FESpaceVel_T::dim + dofP};
  // builderM.buildProblem(dummy, bcsP);
  builderM.buildProblem(timeder, bcsVel);
  builderM.buildProblem(stiffness, bcsVel);
  builderM.buildProblem(grad, bcsVel, bcsP);
  builderM.buildProblem(div, bcsP, bcsVel);
  builderM.closeMatrix();
  Mat<StorageType::RowMajor> matFixed = builderM.A;
  Vec rhsFixed = builderM.b;
  builderM.clear();

  Vec vel{2 * dofU};
  setComponent(vel, feSpaceVel, eqnU.sol.data, eqnU.feSpace, 0);
  setComponent(vel, feSpaceVel, eqnV.sol.data, eqnV.feSpace, 1);

  // uStar / dt + (vel \cdot \nabla) uStar - \nabla \cdot (nu \nabla uStar) = u / dt
  // (uStar, \phi) / dt + ((vel \cdot \nabla) uStar, \phi) + nu (\nabla uStar, \nabla \phi) = (u, \phi) / dt
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceU});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, vel, feSpaceVel, feSpaceU});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyStiffness{nu, feSpaceU});
  eqnUstar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnU.sol.data, feSpaceU});

  // vStar / dt + (vel \cdot \nabla) vStar - \nabla \cdot (nu \nabla vStar) = v / dt
  // (vStar, \phi) / dt + ((vel \cdot \nabla) vStar, \phi) + nu (\nabla vStar, \nabla \phi) = (v, \phi) / dt
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceU});
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, vel, feSpaceVel, feSpaceU});
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyStiffness{nu, feSpaceU});
  eqnVstar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnV.sol.data, feSpaceU});

  // dt \nabla^2 p = \nabla \cdot velStar
  // dt (\nabla p, \nabla q) = - (\nabla \cdot velStar, q)
  eqnP.assemblyListLhs.emplace_back(new AssemblyStiffness{dt, feSpaceP});
  // eqnP lhs does not change in time, we can pre-compute and factorize it
  eqnP.buildLhs();
  eqnP.compute();
  Vec velStar{2*dofU};
  eqnP.assemblyListRhs.emplace_back(new AssemblyDivRhs{-1.0, velStar, feSpaceP, feSpaceVel});

  // u = uStar - dt dp / dx
  // (u, \phi) = (uStar, \phi) - dt (dp / dx, \phi)
  // (u, \phi) = (uStar, \phi) + dt (p, d \phi / dx)
  eqnU.assemblyListLhs.emplace_back(new AssemblyMass{1.0, feSpaceU});
  // eqnU lhs does not change in time, we can pre-compute and factorize it
  eqnU.buildLhs();
  eqnU.compute();
  eqnU.assemblyListRhs.emplace_back(new AssemblyProjection{1.0, eqnUstar.sol.data, feSpaceU});
  eqnU.assemblyListRhs.emplace_back(new AssemblyGradRhs{-dt, eqnP.sol.data, feSpaceU, feSpaceP, {0}});

  // v = vStar - dt dp / dy
  // (v, \phi) = (vStar, \phi) - dt (dp / dy, \phi)
  eqnV.assemblyListLhs.emplace_back(new AssemblyMass{1.0, feSpaceU});
  // eqnV lhs does not change in time, we can pre-compute and factorize it
  eqnV.buildLhs();
  eqnV.compute();
  eqnV.assemblyListRhs.emplace_back(new AssemblyProjection{1.0, eqnVstar.sol.data, feSpaceU});
  eqnV.assemblyListRhs.emplace_back(new AssemblyGradRhs{-dt, eqnP.sol.data, feSpaceU, feSpaceP, {1}});

  Var velM{"velM", 2*dofU + dofP};
  auto ic = [](Vec3 const &) {return Vec2(0., 1.);};
  // auto ic = zero;
  auto ic0 = [&ic] (Vec3 const & p) {return ic(p)[0];};
  auto ic1 = [&ic] (Vec3 const & p) {return ic(p)[1];};
  interpolateAnalyticFunction(ic, feSpaceVel, velM.data);
  interpolateAnalyticFunction(ic0, eqnU.feSpace, eqnU.sol.data);
  interpolateAnalyticFunction(ic1, eqnV.feSpace, eqnV.sol.data);

  eqnUstar.sol.data = eqnU.sol.data;
  eqnVstar.sol.data = eqnV.sol.data;
  t.stop();

  t.start("print");
  Var uM{"uM"};
  Var vM{"vM"};
  getComponent(uM.data, feSpaceU, velM.data, feSpaceVel, 0);
  getComponent(vM.data, feSpaceU, velM.data, feSpaceVel, 1);
  Var pM{"pM", velM.data.block(2*dofU, 0, dofP, 1)};
  IOManager ioV{feSpaceU, "output_nipcs/sol_v"};
  ioV.print({uM, vM, eqnUstar.sol, eqnVstar.sol, eqnU.sol, eqnV.sol});
  IOManager ioP{feSpaceP, "output_nipcs/sol_p"};
  ioP.print({pM, eqnP.sol});
  t.stop();

  IterSolver solverM;

  MilliTimer timerStep;
  double time = 0.0;
  auto const ntime = config["ntime"].as<uint>();
  auto const printStep = config["printStep"].as<uint>();
  for (uint itime=0; itime<ntime; itime++)
  {
    timerStep.start();
    time += dt;
    std::cout << "\n" << separator
              << "solving timestep " << itime+1
              << ", time = " << time << std::endl;
    // filelog << "\n" << separator;

    t.start("monolithic build");
    velOldMonolithic = velM.data;

    builderM.buildProblem(timederRhs, bcsVel);
    builderM.buildProblem(advection, bcsVel);
    builderM.closeMatrix();
    builderM.A += matFixed;
    builderM.b += rhsFixed;
    t.stop();

    t.start("monolithic solve");
    solverM.compute(builderM.A);
    velM.data = solverM.solve(builderM.b);
    auto res = builderM.A * velM.data - builderM.b;
    std::cout << "residual norm: " << res.norm() << std::endl;
    t.stop();

    t.start("ustar build");
    setComponent(vel, feSpaceVel, eqnU.sol.data, eqnU.feSpace, 0);
    setComponent(vel, feSpaceVel, eqnV.sol.data, eqnV.feSpace, 1);
    eqnUstar.build();
    eqnVstar.build();
    t.stop();

    t.start("ustar solve");
    eqnUstar.compute();
    eqnUstar.solve();
    std::cout << "eqnUstar residual norm: " << eqnUstar.residualNorm() << std::endl;
    eqnVstar.compute();
    eqnVstar.solve();
    std::cout << "eqnVstar residual norm: " << eqnVstar.residualNorm() << std::endl;
    t.stop();

    t.start("p build");
    setComponent(velStar, feSpaceVel, eqnUstar.sol.data, eqnUstar.feSpace, 0);
    setComponent(velStar, feSpaceVel, eqnVstar.sol.data, eqnVstar.feSpace, 1);
    eqnP.buildRhs();
    t.stop();

    t.start("p solve");
    eqnP.solve();
    std::cout << "eqnP residual norm: " << eqnP.residualNorm() << std::endl;
    t.stop();

    t.start("u build");
    eqnU.buildRhs();
    eqnV.buildRhs();
    t.stop();

    t.start("u solve");
    eqnU.solve();
    std::cout << "eqnU residual norm: " << eqnU.residualNorm() << std::endl;
    eqnV.solve();
    std::cout << "eqnV residual norm: " << eqnV.residualNorm() << std::endl;
    t.stop();

    t.start("monolithic clear");
    builderM.clear();
    t.stop();

    t.start("split clear");
    eqnUstar.builder.clear();
    eqnVstar.builder.clear();
    eqnP.builder.clear();
    eqnU.builder.clear();
    eqnV.builder.clear();
    t.stop();

    t.start("print");
    if ((itime+1) % printStep == 0)
    {
      std::cout << "printing" << std::endl;
      getComponent(uM.data, feSpaceU, velM.data, feSpaceVel, 0);
      getComponent(vM.data, feSpaceU, velM.data, feSpaceVel, 1);
      ioV.print({uM, vM, eqnUstar.sol, eqnVstar.sol, eqnU.sol, eqnV.sol}, time);

      pM.data = velM.data.block(2*dofU, 0, dofP, 1);
      ioP.print({pM, eqnP.sol}, time);
    }
    t.stop();

    std::cout << "time required: " << timerStep << " ms" << std::endl;
  }
  t.print();

  Var errorU{"errorU"};
  errorU.data = uM.data - eqnU.sol.data;
  Var errorV{"errorV"};
  errorV.data = vM.data - eqnV.sol.data;
  Var errorP{"errorP"};
  errorP.data = pM.data - eqnP.sol.data;

  std::cout << "errorU: " << std::setprecision(16) << errorU.data.norm() << std::endl;
  std::cout << "errorV: " << std::setprecision(16) << errorV.data.norm() << std::endl;
  std::cout << "errorP: " << std::setprecision(16) << errorP.data.norm() << std::endl;

  if (std::fabs(errorU.data.norm() - 0.01275069338696843) > 1.e-12 ||
      std::fabs(errorV.data.norm() - 0.1927239500306806) > 1.e-12 ||
      std::fabs(errorP.data.norm() - 0.4196491952626098) > 1.e-11 )
  {
    std::cerr << "one of the error norms does not coincide with its expected value." << std::endl;
    return 1;
  }
  return 0;
}
