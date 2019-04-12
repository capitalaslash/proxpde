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

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{1., 10., 0.}, {{4, 8, 0}});
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpaceP{*mesh};

  Eqn eqnVelStar{"velStar", feSpaceVel};
  Eqn eqnUstar{"uStar", feSpaceU};
  Eqn eqnVstar{"vStar", feSpaceU};
  Eqn eqnP{"p", feSpaceP};
  Eqn eqnVel{"vel", feSpaceVel};
  Eqn eqnU{"u", feSpaceU};
  Eqn eqnV{"v", feSpaceU};
  t.stop();

  t.start("bc");
  auto zeroS = [] (Vec3 const &) {return 0.;};
  auto oneS = [] (Vec3 const &) {return 1.;};
  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  auto inlet = [] (Vec3 const &) {return Vec2(0.0, 1.0);};
  // auto inlet = [] (Vec3 const &p) {return Vec2(0.0, 6.0*p[0]*(1-p[0]));};
  // auto pIn = [] (Vec3 const &) {return -12.;};

  BCList bcsVel{feSpaceVel};
  bcsVel.addBC(BCEss{feSpaceVel, side::RIGHT, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::LEFT, zero, {0}});
  bcsVel.addBC(BCEss{feSpaceVel, side::BOTTOM, inlet});
  BCList bcsP{feSpaceP};
  bcsP.addBC(BCEss{feSpaceP, side::TOP, [] (Vec3 const &) {return 0.;}});
  // DofSet_T pinSet = {1};
  // bcsP.addBC(BCEss{feSpaceP, pinSet, [] (Vec3 const &) {return 0.;}});

  eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::RIGHT, zeroS});
  eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::LEFT, zeroS});
  eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::BOTTOM, zeroS});
  eqnVstar.bcList.addBC(BCEss{eqnVstar.feSpace, side::RIGHT, zeroS});
  eqnVstar.bcList.addBC(BCEss{eqnVstar.feSpace, side::BOTTOM, oneS});
  eqnVelStar.bcList.addBC(BCEss{eqnVelStar.feSpace, side::RIGHT, zero});
  eqnVelStar.bcList.addBC(BCEss{eqnVelStar.feSpace, side::LEFT, zero, {0}});
  eqnVelStar.bcList.addBC(BCEss{eqnVelStar.feSpace, side::BOTTOM, inlet});
  eqnP.bcList.addBC(BCEss{eqnP.feSpace, side::TOP, [] (Vec3 const &) {return 0.;}});
  // eqnP.bcList.addBC(BCEss{eqnP.feSpace, side::BOTTOM, [] (Vec3 const &) {return 12.;}});
  t.stop();

  t.start("assembly");
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;

  double time = 0.0;
  double const dt = 5.e-2;
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
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceU});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, eqnVel.sol.data, feSpaceU});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyTensorStiffness{nu, feSpaceU});
  eqnUstar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnU.sol.data, feSpaceU});
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceU});
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, eqnVel.sol.data, feSpaceU});
  eqnVstar.assemblyListLhs.emplace_back(new AssemblyTensorStiffness{nu, feSpaceU});
  eqnVstar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnV.sol.data, feSpaceU});
  eqnVelStar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceVel});
  eqnVelStar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, eqnVel.sol.data, feSpaceVel});
  eqnVelStar.assemblyListLhs.emplace_back(new AssemblyTensorStiffness{nu, feSpaceVel});
  eqnVelStar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnVel.sol.data, feSpaceVel});
  // eqnUstar.assemblyListRhs.emplace_back(new AssemblyGradRhs{-1.0, pOld, feSpaceVel, feSpaceP});

  eqnP.assemblyListLhs.emplace_back(new AssemblyStiffness{dt, feSpaceP});
  eqnP.assemblyListRhs.emplace_back(new AssemblyDivRhs{-1.0, eqnVelStar.sol.data, feSpaceP, feSpaceVel});
  // eqnP.assemblyListRhs.emplace_back(new AssemblyStiffnessRhs{dt, pOld, feSpaceP, feSpaceP});

  eqnVel.assemblyListLhs.emplace_back(new AssemblyMass{1.0, feSpaceVel});
  eqnVel.assemblyListRhs.emplace_back(new AssemblyProjection{1.0, eqnVelStar.sol.data, feSpaceVel});
  eqnVel.assemblyListRhs.emplace_back(new AssemblyGradRhs{-dt, eqnP.sol.data, feSpaceVel, feSpaceP});

  BlockVar velM{"velM", {dofU, dofU, dofP}};
  auto ic = [](Vec3 const &) {return Vec2(0., 1.);};
  // auto ic = zero;
  interpolateAnalyticFunction(ic, feSpaceVel, velM.data);
  interpolateAnalyticFunction(zeroS, eqnUstar.feSpace, eqnUstar.sol.data);
  interpolateAnalyticFunction(oneS, eqnVstar.feSpace, eqnVstar.sol.data);
  interpolateAnalyticFunction(ic, feSpaceVel, eqnVel.sol.data);
  eqnVelStar.sol.data = eqnVel.sol.data;
  t.stop();

  t.start("print");
  IOManager ioS{feSpaceU, "output_ipcs/sol_s"};
  ioS.print({eqnUstar.sol, eqnVstar.sol});
  IOManager ioVel{feSpaceVel, "output_ipcs/sol_v"};
  ioVel.print({velM, eqnVelStar.sol, eqnVel.sol});
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
    filelog << "\n" << separator;

    t.start("monolithic");
    velOldMonolithic = velM.data;
    // pOld += eqnP.sol.data;

    builderMonolithic.buildProblem(timederRhs, bcsVel);
    builderMonolithic.buildProblem(advection, bcsVel);
    builderMonolithic.closeMatrix();

    Mat const A = builderMonolithic.A + builderStatic.A;
    Vec const b = builderMonolithic.b + builderStatic.b;
    solverMonolithic.compute(A);
    velM.data = solverMonolithic.solve(b);
    auto res = A * velM.data - b;
    std::cout << "residual norm: " << res.norm() << std::endl;
    builderMonolithic.clear();
    t.stop();

    t.start("split");
    eqnVelStar.build();
    // filelog << "AUstar:\n" << eqnUstar.builder.A.block(dofU, dofU, dofU, dofU) << std::endl;
    // filelog << "bUstar:\n" << eqnUstar.builder.b.block(dofU,0,dofU,1) << std::endl;
    eqnVelStar.solve();
    // filelog << "ustar:\n" << eqnUstar.sol.data.block(dofU,0,dofU,1) << std::endl;
    std::cout << "eqnVelStar residual norm: " << eqnVelStar.residualNorm() << std::endl;

    eqnP.build();
    eqnP.solve();
    std::cout << "eqnP residual norm: " << eqnP.residualNorm() << std::endl;
    t.stop();

    t.start("new split");
    eqnUstar.build();
    eqnUstar.solve();
    std::cout << "eqnUstar residual norm: " << eqnUstar.residualNorm() << std::endl;
    eqnVstar.build();
    eqnVstar.solve();
    std::cout << "eqnUstar residual norm: " << eqnVstar.residualNorm() << std::endl;
    t.stop();

    t.start("split post");
    eqnVel.build();
    // filelog << "AU:\n" << eqnU.builder.A << std::endl;
    // filelog << "bU:\n" << eqnU.builder.b << std::endl;
    eqnVel.solve();
    // filelog << "u:\n" << eqnU.sol.data.block(dofU,0,dofU,1) << std::endl;
    std::cout << "eqnU residual norm: " << eqnVel.residualNorm() << std::endl;

    eqnVelStar.builder.clear();
    eqnP.builder.clear();
    eqnVel.builder.clear();
    t.stop();

    t.start("print");
    if ((itime+1) % printStep == 0)
    {
      std::cout << "printing" << std::endl;
      ioVel.time = time;
      ioVel.iter += 1;
      ioVel.print({velM, eqnVelStar.sol, eqnVel.sol});

      ioP.time = time;
      ioP.iter += 1;
      pM.data = velM.block(2);
      ioP.print({pM, eqnP.sol});
    }
    t.stop();
  }
  t.print();

  double const normM = velM.data.block(0, 0, dofU*FESpaceVel_T::dim, 1).norm();
  double const normS = eqnVel.sol.data.norm();
  Vec velStarNew(2 * dofU);
  velStarNew << eqnUstar.sol.data, eqnVstar.sol.data;
  double const normNew = velStarNew.norm();
  std::cout << "normM: " << std::setprecision(16) << normM << std::endl;
  std::cout << "normS: " << std::setprecision(16) << normS << std::endl;
  std::cout << "normStar: " << std::setprecision(16) << eqnVelStar.sol.data.norm() << std::endl;
  std::cout << "normNew: " << std::setprecision(16) << normNew << std::endl;

  if (std::fabs(normM - 13.3871233713913) > 1.e-12 ||
      std::fabs(normS - 13.33549681268924) > 1.e-12)
  {
    return 1;
  }
  return 0;
}
