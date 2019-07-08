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

// struct EqnUstarData
// {
//   uint const component;
//   double const dt;
//   Vec const & vel;
//   double const nu;
//   Vec const & pOld;
//   FESpaceP_T const & feSpaceP;
// };
//
// struct EqnUstar: Eqn<FESpaceU_T>
// {
//   EqnUstar(std::string_view const name,
//            FESpaceU_T const & fe,
//            EqnUstarData const & data):
//     Eqn<FESpaceU_T>(name, fe)
//   {
//     this->assemblyListLhs.emplace_back<std::unique_ptr<Diagonal<FESpaceU_T>>>(new AssemblyMass{1./data.dt, this->feSpace});
//     this->assemblyListLhs.emplace_back<std::unique_ptr<Diagonal<FESpaceU_T>>>(new AssemblyAdvection{1.0, data.vel, this->feSpace, this->feSpace});
//     this->assemblyListLhs.emplace_back<std::unique_ptr<Diagonal<FESpaceU_T>>>(new AssemblyStiffness{data.nu, this->feSpace});
//     this->assemblyListRhs.emplace_back<std::unique_ptr<Diagonal<FESpaceU_T>>>(new AssemblyProjection{1./data.dt, this->sol.data, this->feSpace});
//     this->assemblyListRhs.emplace_back<std::unique_ptr<Diagonal<FESpaceU_T>>>(new AssemblyGradRhs2{1.0, data.pOld, this->feSpace, data.feSpaceP, {data.component}});
//   }
// };

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
  // builderM.buildLhs(dummy, bcsP);
  builderM.buildLhs(timeder, bcsVel);
  builderM.buildLhs(stiffness, bcsVel);
  builderM.buildCoupling(grad, bcsVel, bcsP);
  builderM.buildCoupling(div, bcsP, bcsVel);
  builderM.closeMatrix();
  Mat<StorageType::RowMajor> matFixed = builderM.A;
  Vec rhsFixed = builderM.b;
  builderM.clear();

  Vec pOld = Vec::Zero(dofP);
  Vec vel{2 * dofU};
  setComponent(vel, feSpaceVel, eqnU.sol.data, eqnU.feSpace, 0);
  setComponent(vel, feSpaceVel, eqnV.sol.data, eqnV.feSpace, 1);

  // uStar / dt + (vel \cdot \nabla) uStar - \nabla \cdot (nu \nabla uStar) = u / dt - d pOld / dx
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceU});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, vel, feSpaceVel, feSpaceU});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyStiffness{nu, feSpaceU});
  eqnUstar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnU.sol.data, feSpaceU});
  eqnUstar.assemblyListRhs.emplace_back(new AssemblyGradRhs2{1.0, pOld, feSpaceU, feSpaceP, {0}});
  // EqnUstar newEqnUstar{"uStar", feSpaceU, {0, dt, vel, nu, pOld, feSpaceP}};

  // vStar / dt + (vel \cdot \nabla) vStar - \nabla \cdot (nu \nabla vStar) = v / dt - d pOld / dy
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceU});
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, vel, feSpaceVel, feSpaceU});
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyStiffness{nu, feSpaceU});
  eqnVstar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnV.sol.data, feSpaceU});
  eqnVstar.assemblyListRhs.emplace_back(new AssemblyGradRhs2{1.0, pOld, feSpaceU, feSpaceP, {1}});
  // eqnVelStar.assemblyListRhs.emplace_back(new AssemblyGradRhs{-1.0, pOld, feSpaceVel, feSpaceP});

  // dt \nabla^2 \delta p = \nabla \cdot velStar
  eqnP.assemblyListLhs.emplace_back(new AssemblyStiffness{dt, feSpaceP});
  // eqnP lhs does not change in time, we can pre-compute and factorize it
  eqnP.buildLhs();
  eqnP.compute();
  Vec velStar{2*dofU};
  eqnP.assemblyListRhs.emplace_back(new AssemblyDivRhs{-1.0, velStar, feSpaceP, feSpaceVel});
  // eqnP.assemblyListRhs.emplace_back(new AssemblyStiffnessRhs{-dt, pOld, feSpaceP, feSpaceP});

  // pOld += \delta p
  // u = uStar - dt d \delta p / dx
  eqnU.assemblyListLhs.emplace_back(new AssemblyMass{1.0, feSpaceU});
  // eqnU lhs does not change in time, we can pre-compute and factorize it
  eqnU.buildLhs();
  eqnU.compute();
  eqnU.assemblyListRhs.emplace_back(new AssemblyProjection{1.0, eqnUstar.sol.data, feSpaceU});
  eqnU.assemblyListRhs.emplace_back(new AssemblyGradRhs{-dt, eqnP.sol.data, feSpaceU, feSpaceP, {0}});

  // v = vStar - dt d \delta p / dy
  eqnV.assemblyListLhs.emplace_back(new AssemblyMass{1.0, feSpaceU});
  // eqnV lhs does not change in time, we can pre-compute and factorize it
  eqnV.buildLhs();
  eqnV.compute();
  eqnV.assemblyListRhs.emplace_back(new AssemblyProjection{1.0, eqnVstar.sol.data, feSpaceU});
  eqnV.assemblyListRhs.emplace_back(new AssemblyGradRhs{-dt, eqnP.sol.data, feSpaceU, feSpaceP, {1}});

  Var velM{"velM", 2*dofU + dofP};
  auto const ic = [](Vec3 const &) {return Vec2(0., 1.);};
  // auto const ic = [](Vec3 const & p) {return Vec2{0., 1.5 * (1. - p(0)*p(0))};};
  // auto const ic = zero;
  auto const ic0 = [&ic] (Vec3 const & p) {return ic(p)[0];};
  auto const ic1 = [&ic] (Vec3 const & p) {return ic(p)[1];};
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
  Var pSplit{"p"};
  pSplit.data = Vec::Zero(eqnP.feSpace.dof.size);
  IOManager ioV{feSpaceU, "output_ipcs/sol_v"};
  ioV.print({uM, vM, eqnUstar.sol, eqnVstar.sol, eqnU.sol, eqnV.sol});
  IOManager ioP{feSpaceP, "output_ipcs/sol_p"};
  ioP.print({pM, pSplit});
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

    builderM.buildRhs(timederRhs, bcsVel);
    builderM.buildLhs(advection, bcsVel);
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
    pOld += eqnP.sol.data;
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
    pSplit.data += eqnP.sol.data;
    if ((itime+1) % printStep == 0)
    {
      std::cout << "printing" << std::endl;
      getComponent(uM.data, feSpaceU, velM.data, feSpaceVel, 0);
      getComponent(vM.data, feSpaceU, velM.data, feSpaceVel, 1);
      ioV.print({uM, vM, eqnUstar.sol, eqnVstar.sol, eqnU.sol, eqnV.sol}, time);

      pM.data = velM.data.block(2*dofU, 0, dofP, 1);
      ioP.print({pM, pSplit}, time);
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
  errorP.data = pM.data - pOld;

  std::cout << "errorU: " << std::setprecision(16) << errorU.data.norm() << std::endl;
  std::cout << "errorV: " << std::setprecision(16) << errorV.data.norm() << std::endl;
  std::cout << "errorP: " << std::setprecision(16) << errorP.data.norm() << std::endl;

  if (std::fabs(errorU.data.norm() - 5.511440044739265e-05) > 1.e-12 ||
      std::fabs(errorV.data.norm() - 7.500752437417739e-05) > 1.e-12 ||
      std::fabs(errorP.data.norm() - 6.141568682924751e-05) > 1.e-11 )
  {
    std::cerr << "one of the error norms does not coincide with its expected value." << std::endl;
    return 1;
  }
  return 0;
}
