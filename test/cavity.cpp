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

#include <iostream>

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpaceU_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR>;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;

int main(int argc, char* argv[])
{
  uint const numPts_x = (argc < 3)? 21 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 21 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});

  FESpaceVel_T feSpaceVel{meshPtr};
  FESpaceP_T feSpaceP{meshPtr};

  BCList<FESpaceVel_T> bcsVel{feSpaceVel};
  bcsVel.addEssentialBC(side::RIGHT, [] (Vec3 const &) {return Vec2::Constant(0.);}, {0,1});
  bcsVel.addEssentialBC(side::LEFT, [] (Vec3 const &) {return Vec2::Constant(0.);}, {0,1});
  bcsVel.addEssentialBC(side::BOTTOM, [] (Vec3 const &) {return Vec2::Constant(0.);}, {0,1});
  bcsVel.addEssentialBC(side::TOP, [] (Vec3 const &) {return Vec2(1.0, 0.0);}, {0, 1});
  BCList<FESpaceP_T> bcsP{feSpaceP};

  auto const dofU = feSpaceVel.dof.totalNum;
  auto const dofP = feSpaceP.dof.totalNum;
  uint const numDOFs = dofU*FESpaceVel_T::dim + dofP;
  Mat A(numDOFs, numDOFs);
  Vec b = Vec::Zero(numDOFs);

  AssemblyStiffness<FESpaceVel_T> stiffness(1.0, feSpaceVel);
  AssemblyVGrad<FESpaceVel_T, FESpaceP_T> grad(feSpaceVel, feSpaceP, 0, 2*dofU);
  AssemblyVDiv<FESpaceP_T, FESpaceVel_T> div(feSpaceP, feSpaceVel, 2*dofU, 0);

  Builder builder(A, b);
  builder.buildProblem(stiffness, bcsVel);
  builder.buildProblem(grad, bcsVel, bcsP);
  builder.buildProblem(div, bcsP, bcsVel);
  builder.closeMatrix();

  Vec sol(2*dofU + dofP);
  // Eigen::SparseLU<Mat> solver;
  // solver.analyzePattern(A);
  // solver.factorize(A);
  // Eigen::GMRES<Mat> solver(A);
  Eigen::UmfPackLU<Mat> solver(A);
  // Eigen::SuperLU<Mat> solver(A);
  sol = solver.solve(b);

  // std::cout << "A:\n" << A << std::endl;
  // std::cout << "b:\n" << b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  std::cout << sol.norm() << std::endl;

  Var vel{"vel", sol};
  Var u{"u", sol, 0, dofU};
  Var v{"v", sol, dofU, dofU};
  Var p{"p", sol, 2*dofU, dofP};

  IOManager<FESpaceVel_T> ioU{feSpaceVel, "sol_cavity_vel.xmf"};
  ioU.print({vel});
  IOManager<FESpaceP_T> ioP{feSpaceP, "sol_cavity_p.xmf"};
  ioP.print({p});

  // Var exact{"exact", feSpaceU.dof.totalNum};
  // interpolateAnalyticFunction(exact_sol, feSpaceU, exact.data);
//  Var error{"e"};
//  error.data = sol /*- exact.data*/;

//  double norm = error.data.norm();
//  std::cout << "the norm of the error is " << norm << std::endl;
//  if(std::fabs(norm - 6.03323) > 1.e-4)
//  {
//    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
//    return 1;
//  }

  return 0;
}
