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

  uint const numElemX = (argc == 3) ? std::atoi(argv[1]) : 4;
  uint const numElemY = (argc == 3) ? std::atoi(argv[2]) : 8;

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{1., 10., 0.}, {{numElemX, numElemY, 0}});
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpaceP{*mesh};

  Eqn eqnUstar{"uStar", feSpaceU};
  Eqn eqnVstar{"vStar", feSpaceU};

  Eqn<FESpaceP_T, StorageType::ClmMajor> eqnP{"p", feSpaceP};

  Eqn<FESpaceU_T, StorageType::ClmMajor> eqnU{"u", feSpaceU};
  Eqn<FESpaceU_T, StorageType::ClmMajor> eqnV{"v", feSpaceU};
  t.stop();

  t.start("bc");
  auto zeroS = [] (Vec3 const &) {return 0.;};
  auto zeroV = [] (Vec3 const &) {return Vec2::Constant(0.);};
  // auto inlet = [] (Vec3 const &) {return Vec2(0.0, 1.0);};
  auto inlet = [] (Vec3 const &p) {return Vec2(0.0, 1.5 * (1. - p[0]*p[0]));};
  auto inlet0 = [&inlet] (Vec3 const & p) {return inlet(p)[0];};
  auto inlet1 = [&inlet] (Vec3 const & p) {return inlet(p)[1];};
  // auto pIn = [] (Vec3 const &) {return 12.;};

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

  double time = 0.0;
  double const dt = 5e-2;
  uint const ntime = 100;
  uint const printStep = 1;
  double const nu = 1.e-1;

  Vec velOldMonolithic(2*dofU);
  AssemblyMass timeder(1./dt, feSpaceVel);
  AssemblyAdvection advection(1.0, velOldMonolithic, feSpaceVel);
  AssemblyTensorStiffness stiffness(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);
  AssemblyProjection timederRhs(1./dt, velOldMonolithic, feSpaceVel);
  // AssemblyBCNormal naturalBC{pIn, side::BOTTOM, feSpaceVel};
  // to apply bc on pressure
  AssemblyMass dummy{0.0, feSpaceP, {0}, 2*dofU, 2*dofU};

  Builder builderStatic{dofU*FESpaceVel_T::dim + dofP};
  builderStatic.buildProblem(dummy, bcsP);
  builderStatic.buildProblem(timeder, bcsVel);
  builderStatic.buildProblem(stiffness, bcsVel);
  builderStatic.buildProblem(grad, bcsVel, bcsP);
  builderStatic.buildProblem(div, bcsP, bcsVel);
  builderStatic.closeMatrix();

  Vec pOld = Vec::Zero(dofP);
  Vec vel(2 * dofU);
  vel << eqnU.sol.data, eqnV.sol.data;
  // uStar / dt + (vel \cdot \nabla) uStar - \nabla \cdot (nu \nabla uStar) = u / dt - d pOld / dx
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceU});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, vel, feSpaceU});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyStiffness{nu, feSpaceU});
  eqnUstar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnU.sol.data, feSpaceU});
  eqnUstar.assemblyListRhs.emplace_back(new AssemblyGradRhs2{1.0, pOld, feSpaceU, feSpaceP, {0}});

  // vStar / dt + (vel \cdot \nabla) vStar - \nabla \cdot (nu \nabla vStar) = v / dt - d pOld / dy
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceU});
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, vel, feSpaceU});
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

  BlockVar velM{"velM", {dofU, dofU, dofP}};
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
  Var uM{"uM", velM.block(0)};
  Var vM{"vM", velM.block(1)};
  IOManager ioV{feSpaceU, "output_ipcs/sol_v"};
  ioV.print({uM, vM, eqnUstar.sol, eqnVstar.sol, eqnU.sol, eqnV.sol});
  IOManager ioP{feSpaceP, "output_ipcs/sol_p"};
  Var pM{"pM", velM.block(2)};
  ioP.print({pM, eqnP.sol});
  t.stop();

  Builder builderMonolithic{dofU*FESpaceVel_T::dim + dofP};
  LUSolver solverMonolithic;

  for (uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "\n" << separator
              << "solving timestep " << itime+1
              << ", time = " << time << std::endl;
    // filelog << "\n" << separator;

    t.start("monolithic build");
    velOldMonolithic = velM.data;
    pOld += eqnP.sol.data;

    builderMonolithic.buildProblem(timederRhs, bcsVel);
    builderMonolithic.buildProblem(advection, bcsVel);
    builderMonolithic.closeMatrix();

    Mat<StorageType::ClmMajor> const A = builderMonolithic.A + builderStatic.A;
    Vec const b = builderMonolithic.b + builderStatic.b;
    t.stop();

    t.start("monolithic solve");
    solverMonolithic.compute(A);
    velM.data = solverMonolithic.solve(b);
    auto res = A * velM.data - b;
    std::cout << "residual norm: " << res.norm() << std::endl;
    t.stop();

    t.start("ustar build");
    vel << eqnU.sol.data, eqnV.sol.data;
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
    velStar << eqnUstar.sol.data, eqnVstar.sol.data;
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
    builderMonolithic.clear();
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
      uM.data = velM.block(0);
      vM.data = velM.block(1);
      ioV.time = time;
      ioV.iter += 1;
      ioV.print({uM, vM, eqnUstar.sol, eqnVstar.sol, eqnU.sol, eqnV.sol});

      pM.data = velM.block(2);
      ioP.time = time;
      ioP.iter += 1;
      ioP.print({pM, eqnP.sol});
    }
    t.stop();
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

  if (std::fabs(errorU.data.norm() - 1.663552566403615e-05) > 1.e-12 ||
      std::fabs(errorV.data.norm() - 2.087262887058092e-05) > 1.e-12 ||
      std::fabs(errorP.data.norm() - 2.477471412883698e-05) > 1.e-12 )
  {
    std::cerr << "one of the error norms does not coincide with its expected value." << std::endl;
    return 1;
  }
  return 0;
}
