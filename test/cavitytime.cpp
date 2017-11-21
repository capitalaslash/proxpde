#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "assembler.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

#include <experimental/filesystem>

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;

int main(int argc, char* argv[])
{
  MilliTimer t;
  uint const numPts_x = (argc < 3)? 3 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 3 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};

  t.start();
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpaceVel_T feSpaceVel{meshPtr};
  FESpaceP_T feSpaceP{meshPtr};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  BCList<FESpaceVel_T> bcsVel{feSpaceVel};
  bcsVel.addEssentialBC(side::RIGHT, zero, {0,1});
  bcsVel.addEssentialBC(side::LEFT, zero, {0,1});
  bcsVel.addEssentialBC(side::BOTTOM, zero, {0,1});
  bcsVel.addEssentialBC(side::TOP, [] (Vec3 const &) {return Vec2(1.0, 0.0);}, {0, 1});
  BCList<FESpaceP_T> bcsP{feSpaceP};
  // DofSet_T pinSet = {1};
  // bcsP.addEssentialBC(pinSet, [] (Vec3 const &) {return 0.;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  // t.start();
  auto const dofU = feSpaceVel.dof.totalNum;
  auto const dofP = feSpaceP.dof.totalNum;
  uint const numDOFs = dofU*FESpaceVel_T::dim + dofP;

  AssemblyStiffness<FESpaceVel_T> stiffness(1.e-3, feSpaceVel);
  AssemblyGrad<FESpaceVel_T, FESpaceP_T> grad(feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv<FESpaceP_T, FESpaceVel_T> div(feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);

  double const dt = 5.e-1;
  AssemblyMass<FESpaceVel_T> timeder(1./dt, feSpaceVel);
  Vec vel_old(2*dofU);
  AssemblyVecRhs<FESpaceVel_T> timeder_rhs(vel_old, feSpaceVel);

  uint const ntime = 10;

  Var vel{"vel"};
  vel.data = Vec::Zero(2*dofU + dofP);
  auto ic = [](Vec3 const &) {return Vec2(1., 0.);};
  interpolateAnalyticFunction(ic, feSpaceVel, vel.data);

  std::experimental::filesystem::create_directory("output");
  Eigen::UmfPackLU<Mat> solver;
  IOManager<FESpaceVel_T> ioVel{feSpaceVel, "output/sol_cavitytime_v_0.xmf", 0.0};
  ioVel.print({vel});
  IOManager<FESpaceP_T> ioP{feSpaceP, "output/sol_cavitytime_p_0.xmf", 0.0};
  Var p{"p", vel.data, 2*dofU, dofP};
  ioP.print({p});

  for (uint itime=0; itime<ntime; itime++)
  {
    std::cout << "solving timestep " << itime << std::endl;

    vel_old = vel.data / dt;

    Mat mat(numDOFs, numDOFs);
    Vec b = Vec::Zero(numDOFs);
    Builder builder(mat, b);
    builder.buildProblem(timeder, bcsVel);
    builder.buildProblem(timeder_rhs, bcsVel);
    builder.buildProblem(stiffness, bcsVel);
    builder.buildProblem(grad, bcsVel, bcsP);
    builder.buildProblem(div, bcsP, bcsVel);
    builder.closeMatrix();

    solver.compute(mat);
    vel.data = solver.solve(b);

    ioVel.fileName = "output/sol_cavitytime_v_" + std::to_string(itime) + ".xmf";
    ioVel.time = (itime+1) * dt;
    ioVel.print({vel});
    p.data = vel.data.block(2*dofU,0,dofP,1);
    ioP.fileName = "output/sol_cavitytime_p_" + std::to_string(itime) + ".xmf";
    ioP.time = (itime+1) * dt;
    ioP.print({p});
  }

  return 0;
}
