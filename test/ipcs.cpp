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

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;

int main(int argc, char* argv[])
{
  MilliTimer t;

  t.start();
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{1., 10., 0.}, {{2, 1, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh};

  Eqn<Mesh_T, QuadraticRefFE, QuadraticQR, 2> eqnUstar{"uStar", *mesh};
  Eqn<Mesh_T, LinearRefFE, QuadraticQR> eqnP{"p", *mesh};
  Eqn<Mesh_T, QuadraticRefFE, QuadraticQR, 2> eqnU{"u", *mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  // auto inlet = [] (Vec3 const &) {return Vec2(0.0, 1.0);};
  // auto inlet = [] (Vec3 const &p) {return Vec2(0.0, 6.0*p[0]*(1-p[0]));};
  auto pIn = [] (Vec3 const &) {return -12.;};

  BCList bcsVel{feSpaceVel};
  bcsVel.addBC(BCEss{feSpaceVel, side::RIGHT, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::LEFT, zero});
  // bcsVel.addBC(BCEss{feSpaceVel, side::BOTTOM, inlet});
  BCList bcsP{feSpaceP};
  bcsP.addBC(BCEss{feSpaceP, side::TOP, [] (Vec3 const &) {return 0.;}});
  // DofSet_T pinSet = {1};
  // bcsP.addBC(BCEss{feSpaceP, pinSet, [] (Vec3 const &) {return 0.;}});

  eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::RIGHT, zero});
  eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::LEFT, zero});
  // eqnUstar.bcList.addBC(BCEss{eqnUstar.feSpace, side::BOTTOM, inlet});
  eqnP.bcList.addBC(BCEss{eqnP.feSpace, side::TOP, [] (Vec3 const &) {return 0.;}});
  eqnP.bcList.addBC(BCEss{eqnP.feSpace, side::BOTTOM, [] (Vec3 const &) {return 12.;}});
  eqnU.bcList.addBC(BCEss{eqnU.feSpace, side::RIGHT, zero});
  eqnU.bcList.addBC(BCEss{eqnU.feSpace, side::LEFT, zero});
  // eqnU.bcList.addBC(BCEss{eqnU.feSpace, side::BOTTOM, inlet});
  std::cout << "bcs: " << t << " ms" << std::endl;

  // t.start();
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
  AssemblyBCNormal naturalBC{pIn, side::BOTTOM, feSpaceVel};
  AssemblyMass dummy{0.0, feSpaceP, {0}, 2*dofU, 2*dofU};

  Vec pOld = Vec::Zero(dofP);
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyMass{1./dt, feSpaceVel});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyAdvection{1.0, eqnU.sol.data, feSpaceVel});
  eqnUstar.assemblyListLhs.emplace_back(new AssemblyTensorStiffness{nu, feSpaceVel});
  eqnUstar.assemblyListRhs.emplace_back(new AssemblyProjection{1./dt, eqnU.sol.data, feSpaceVel});
  // eqnUstar.assemblyListRhs.emplace_back(new AssemblyGradRhs{-1.0, pOld, feSpaceVel, feSpaceP});

  eqnP.assemblyListLhs.emplace_back(new AssemblyStiffness{dt, feSpaceP});
  eqnP.assemblyListRhs.emplace_back(new AssemblyDivRhs{-1.0, eqnUstar.sol.data, feSpaceP, feSpaceVel});
  // eqnP.assemblyListRhs.emplace_back(new AssemblyStiffnessRhs{dt, pOld, feSpaceP, feSpaceP});

  eqnU.assemblyListLhs.emplace_back(new AssemblyMass{1.0, feSpaceVel});
  eqnU.assemblyListRhs.emplace_back(new AssemblyProjection{1.0, eqnUstar.sol.data, feSpaceVel});
  eqnU.assemblyListRhs.emplace_back(new AssemblyGradRhs{-dt, eqnP.sol.data, feSpaceVel, feSpaceP});

  BlockVar solMonolithic{"sol", {dofU, dofU, dofP}};
  // auto ic = [](Vec3 const &) {return Vec2(0., 1.);};
  auto ic = zero;
  interpolateAnalyticFunction(ic, feSpaceVel, solMonolithic.data);
  interpolateAnalyticFunction(ic, feSpaceVel, eqnU.sol.data);
  eqnUstar.sol.data = eqnU.sol.data;

  IOManager ioVel{feSpaceVel, "output_ipcs/sol_v"};
  ioVel.print({solMonolithic, eqnUstar.sol, eqnU.sol});
  IOManager ioP{feSpaceP, "output_ipcs/sol_p"};
  Var pMonolithic{"pMonolithic", solMonolithic.block(2)};
  ioP.print({pMonolithic, eqnP.sol});

  Builder builderMonolithic{dofU*FESpaceVel_T::dim + dofP};
  LUSolver solverMonolithic;

  for (uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "\n" << separator
              << "solving timestep " << itime+1
              << ", time = " << time << std::endl;
    filelog << "\n" << separator;

    velOldMonolithic = solMonolithic.data;
    // pOld += eqnP.sol.data;

    builderMonolithic.buildProblem(dummy, bcsP);
    builderMonolithic.buildProblem(timeder, bcsVel);
    builderMonolithic.buildProblem(timederRhs, bcsVel);
    builderMonolithic.buildProblem(advection, bcsVel);
    builderMonolithic.buildProblem(stiffness, bcsVel);
    builderMonolithic.buildProblem(grad, bcsVel, bcsP);
    builderMonolithic.buildProblem(div, bcsP, bcsVel);
    builderMonolithic.buildProblem(naturalBC, bcsVel);
    builderMonolithic.closeMatrix();

    solverMonolithic.compute(builderMonolithic.A);
    solMonolithic.data = solverMonolithic.solve(builderMonolithic.b);
    auto res = builderMonolithic.A * solMonolithic.data - builderMonolithic.b;
    std::cout << "residual norm: " << res.norm() << std::endl;
    builderMonolithic.clear();

    eqnUstar.build();
    // filelog << "AUstar:\n" << eqnUstar.builder.A.block(dofU, dofU, dofU, dofU) << std::endl;
    // filelog << "bUstar:\n" << eqnUstar.builder.b.block(dofU,0,dofU,1) << std::endl;
    eqnUstar.solve();
    // filelog << "ustar:\n" << eqnUstar.sol.data.block(dofU,0,dofU,1) << std::endl;
    std::cout << "eqnUstar residual norm: " << eqnUstar.residualNorm() << std::endl;

    eqnP.build();
    filelog << "Ap:\n" << eqnP.builder.A << std::endl;
    filelog << "bp:\n" << eqnP.builder.b << std::endl;
    eqnP.solve();
    filelog << "p:\n" << eqnP.sol.data << std::endl;
    std::cout << "eqnP residual norm: " << eqnP.residualNorm() << std::endl;

    eqnU.build();
    // filelog << "AU:\n" << eqnU.builder.A << std::endl;
    // filelog << "bU:\n" << eqnU.builder.b << std::endl;
    eqnU.solve();
    // filelog << "u:\n" << eqnU.sol.data.block(dofU,0,dofU,1) << std::endl;
    std::cout << "eqnU residual norm: " << eqnU.residualNorm() << std::endl;

    eqnUstar.builder.clear();
    eqnP.builder.clear();
    eqnU.builder.clear();

    if ((itime+1) % printStep == 0)
    {
      std::cout << "printing" << std::endl;
      ioVel.time = time;
      ioVel.iter += 1;
      ioVel.print({solMonolithic, eqnUstar.sol, eqnU.sol});

      ioP.time = time;
      ioP.iter += 1;
      pMonolithic.data = solMonolithic.block(2);
      ioP.print({pMonolithic, eqnP.sol});
    }
  }
  return 0;
}
