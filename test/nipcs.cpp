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
    config["nx"] = 4U;
    config["ny"] = 8U;
    config["dt"] = 0.1;
    config["ntime"]= 50U;
    config["nu"] = 0.1;
    config["printStep"] = 1U;
  }

  config.validate({"nx", "ny", "dt", "ntime", "nu", "printStep"});

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
  FESpaceP_T feSpaceP{*mesh, feSpaceVel.dof.size * FESpaceVel_T::dim};
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpacePSplit{*mesh};
  t.stop();

  t.start("bc");
  auto const zeroS = [] (Vec3 const &) { return 0.; };
  auto const zeroV = [] (Vec3 const &) { return Vec2{0., 0.}; };
  // auto const inlet = [] (Vec3 const &) { return Vec2(0.0, 1.0); };
  auto const inlet = [] (Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p[0]*p[0])}; };
  auto const inlet0 = [&inlet] (Vec3 const & p) { return inlet(p)[0]; };
  auto const inlet1 = [&inlet] (Vec3 const & p) { return inlet(p)[1]; };
  // auto const pIn = [] (Vec3 const &) {return 3. * hy * nu;};

  // last essential bc wins on corners
  auto bcsVel = std::make_tuple(
        BCEss{feSpaceVel, side::BOTTOM},
        BCEss{feSpaceVel, side::RIGHT},
        BCEss{feSpaceVel, side::LEFT, {0}});
  std::get<0>(bcsVel) << inlet;
  std::get<1>(bcsVel) << zeroV;
  std::get<2>(bcsVel) << zeroV;

  auto bcTopP = BCEss{feSpaceP, side::TOP};
  bcTopP << zeroS;
  auto const bcsP = std::tuple{bcTopP};
  // DofSet_T pinSet = {1};
  // bcsP.addBC(BCEss{feSpaceP, pinSet, zeroS});

  auto bcsUStar = std::make_tuple(
        BCEss{feSpaceU, side::BOTTOM},
        BCEss{feSpaceU, side::RIGHT},
        BCEss{feSpaceU, side::LEFT});
  std::get<0>(bcsUStar) << inlet0;
  std::get<1>(bcsUStar) << zeroS;
  std::get<2>(bcsUStar) << zeroS;

  // last essential bc wins on corners
  auto bcsVStar = std::make_tuple(
        BCEss{feSpaceU, side::BOTTOM},
        BCEss{feSpaceU, side::RIGHT});
  std::get<0>(bcsVStar) << inlet1;
  std::get<1>(bcsVStar) << zeroS;

  auto bcTopPSplit = BCEss{feSpacePSplit, side::TOP};
  bcTopPSplit << zeroS;
  auto const bcsPSplit = std::tuple{bcTopPSplit};
  // eqnP.bcList.addBC(BCEss{eqnP.feSpace, side::BOTTOM, pIn});
  t.stop();

  t.start("assembly");
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;

  auto const dt = config["dt"].as<double>();
  auto const nu = config["nu"].as<double>();

  Vec velOldMonolithic{2*dofU};
  AssemblyScalarMass timeder(1./dt, feSpaceVel);
  AssemblyAdvection advection(1.0, velOldMonolithic, feSpaceVel, feSpaceVel);
  AssemblyTensorStiffness stiffness(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);
  AssemblyProjection timederRhs(1./dt, velOldMonolithic, feSpaceVel);
  // AssemblyBCNormal naturalBC{pIn, side::BOTTOM, feSpaceVel};
  // to apply bc on pressure
  // AssemblyDummy dummy{feSpaceP};

  Builder<StorageType::RowMajor> builderM{dofU*FESpaceVel_T::dim + dofP};
  // builderM.buildLhs(dummy, bcsP);
  builderM.buildLhs(std::tuple{timeder, stiffness}, bcsVel);
  builderM.buildCoupling(grad, bcsVel, bcsP);
  builderM.buildCoupling(div, bcsP, bcsVel);
  builderM.closeMatrix();
  Mat<StorageType::RowMajor> matFixed = builderM.A;
  Vec rhsFixed = builderM.b;
  builderM.clear();

  Var uStar{"uStar", feSpaceU.dof.size};
  Var vStar{"vStar", feSpaceU.dof.size};
  Var p{"p", feSpacePSplit.dof.size};
  Var u{"u", feSpaceU.dof.size};
  Var v{"v", feSpaceU.dof.size};

  Vec velStar{2 * dofU};
  Vec vel{2 * dofU};

  // uStar / dt + (vel \cdot \nabla) uStar - \nabla \cdot (nu \nabla uStar) = u / dt
  // (uStar, \phi) / dt + ((vel \cdot \nabla) uStar, \phi) + nu (\nabla uStar, \nabla \phi) = (u, \phi) / dt
  auto const uStarLhs = std::make_tuple(
        AssemblyScalarMass{1./dt, feSpaceU},
        AssemblyAdvection{1.0, vel, feSpaceVel, feSpaceU},
        AssemblyStiffness{nu, feSpaceU});
  auto const uStarRhs = std::make_tuple(
        AssemblyProjection{1./dt, u.data, feSpaceU});

  // vStar / dt + (vel \cdot \nabla) vStar - \nabla \cdot (nu \nabla vStar) = v / dt
  // (vStar, \phi) / dt + ((vel \cdot \nabla) vStar, \phi) + nu (\nabla vStar, \nabla \phi) = (v, \phi) / dt
  auto const vStarLhs = std::make_tuple(
        AssemblyScalarMass{1./dt, feSpaceU},
        AssemblyAdvection{1.0, vel, feSpaceVel, feSpaceU},
        AssemblyStiffness{nu, feSpaceU});
  auto const vStarRhs = std::make_tuple(
        AssemblyProjection{1./dt, v.data, feSpaceU});

  // dt \nabla^2 p = \nabla \cdot velStar
  // dt (\nabla p, \nabla q) = - (\nabla \cdot velStar, q)
  auto const pLhs = std::make_tuple(AssemblyStiffness{dt, feSpacePSplit});
  auto const pRhs = std::make_tuple(AssemblyDivRhs{-1.0, velStar, feSpaceVel, feSpacePSplit});
  // eqnP lhs does not change in time, we can pre-compute and factorize it

  // u = uStar - dt dp / dx
  // (u, \phi) = (uStar, \phi) - dt (dp / dx, \phi)
  // (u, \phi) = (uStar, \phi) + dt (p, d \phi / dx)
  auto const uLhs = std::make_tuple(AssemblyScalarMass{1.0, feSpaceU});
  auto const uRhs = std::make_tuple(
        AssemblyProjection{1.0, uStar.data, feSpaceU},
        AssemblyGradRhs{-dt, p.data, feSpacePSplit, feSpaceU, {0}});

  // v = vStar - dt dp / dy
  // (v, \phi) = (vStar, \phi) - dt (dp / dy, \phi)
  auto const vLhs = std::make_tuple(AssemblyScalarMass{1.0, feSpaceU});
  auto const vRhs = std::make_tuple(
        AssemblyProjection{1.0, vStar.data, feSpaceU},
        AssemblyGradRhs{-dt, p.data, feSpacePSplit, feSpaceU, {1}});

  t.start("eqn");
  Eqn eqnUstar{uStar, uStarLhs, uStarRhs, bcsUStar, EqnSolver<StorageType::RowMajor>()};
  Eqn eqnVstar{vStar, vStarLhs, vStarRhs, bcsVStar, EqnSolver<StorageType::RowMajor>()};

  Eqn eqnP{p, pLhs, pRhs, bcsPSplit};
  // eqnP lhs does not change in time, we can pre-compute and factorize it
  eqnP.buildLhs();
  eqnP.compute();

  Eqn eqnU{u, uLhs, uRhs, std::make_tuple()};
  // eqnU lhs does not change in time, we can pre-compute and factorize it
  eqnU.buildLhs();
  eqnU.compute();
  Eqn eqnV{v, vLhs, vRhs, std::make_tuple()};
  // eqnV lhs does not change in time, we can pre-compute and factorize it
  eqnV.buildLhs();
  eqnV.compute();
  t.stop();

  t.start("ic");
  Var velM{"velM"};
  velM.data = Vec::Zero( 2*dofU + dofP);
  auto ic = [](Vec3 const &) {return Vec2(0., 1.);};
  // auto ic = zero;
  auto ic0 = [&ic] (Vec3 const & p) {return ic(p)[0];};
  auto ic1 = [&ic] (Vec3 const & p) {return ic(p)[1];};
  interpolateAnalyticFunction(ic, feSpaceVel, velM.data);
  interpolateAnalyticFunction(ic0, feSpaceU, u.data);
  interpolateAnalyticFunction(ic1, feSpaceU, v.data);

  uStar.data = u.data;
  vStar.data = v.data;
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

    builderM.buildRhs(std::tuple{timederRhs}, bcsVel);
    builderM.buildLhs(std::tuple{advection}, bcsVel);
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
    setComponent(vel, feSpaceVel, u.data, feSpaceU, 0);
    setComponent(vel, feSpaceVel, v.data, feSpaceU, 1);
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
    setComponent(velStar, feSpaceVel, uStar.data, feSpaceU, 0);
    setComponent(velStar, feSpaceVel, vStar.data, feSpaceU, 1);
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

  double const errorNormU = errorU.data.norm();
  double const errorNormV = errorV.data.norm();
  double const errorNormP = errorP.data.norm();
  std::cout << "errorU: " << std::setprecision(16) << errorNormU << std::endl;
  std::cout << "errorV: " << std::setprecision(16) << errorNormV << std::endl;
  std::cout << "errorP: " << std::setprecision(16) << errorNormP << std::endl;

  return checkError(
    {errorNormU, errorNormV, errorNormP},
    {0.01275069338696843, 0.1927239500306806, 0.4196491952626098},
    1.e-11);
}
