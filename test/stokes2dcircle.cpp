#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

#include <iostream>

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpaceU_T = FESpace<Mesh_T,
                           FEType<Elem_T,2>::RefFE_T,
                           FEType<Elem_T,2>::RecommendedQR>;
using FESpaceP_T = FESpace<Mesh_T,
                           FEType<Elem_T,1>::RefFE_T,
                           FEType<Elem_T,2>::RecommendedQR>;

int main(int argc, char* argv[])
{
  uint const numElemsX = (argc < 4)? 2 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 4)? 2 : std::stoi(argv[2]);
  uint const numElemsR = (argc < 4)? 2 : std::stoi(argv[3]);

  Vec3 const origin{0., 0., 0.};
  double const radius = 1.;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildCircleMesh(*mesh, origin, radius, {{numElemsX, numElemsY, numElemsR}});
  // std::cout << *mesh << std::endl;

  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpaceP{*mesh};

  auto zeroFun = [] (Vec3 const&) {return 0.;};
  auto oneFun = [] (Vec3 const&) {return 1.;};
  BCList bcsU{feSpaceU};
  bcsU.addBC(BCEss{feSpaceU, side::CIRCLE, zeroFun});
  bcsU.addBC(BCEss{feSpaceU, side::TOP, zeroFun});
  bcsU.addBC(BCNat<FESpaceU_T>{side::LEFT, oneFun});
  BCList bcsV{feSpaceU};
  bcsV.addBC(BCEss{feSpaceU, side::CIRCLE, zeroFun});
  bcsV.addBC(BCEss{feSpaceU, side::LEFT, zeroFun});
  BCList bcsP(feSpaceP);

  uint const numDOFs = 2*feSpaceU.dof.size + feSpaceP.dof.size;

  AssemblyStiffness stiffnessU(1.0, feSpaceU);
  AssemblyStiffness stiffnessV(1.0, feSpaceU, {1}, feSpaceU.dof.size, feSpaceU.dof.size);
  // AssemblyGrad gradU(0, feSpaceU, feSpaceP, 0, 2*feSpaceU.dof.size);
  AssemblyDiv divU(-1.0, feSpaceP, feSpaceU, {0}, 2*feSpaceU.dof.size, 0);
  AssemblyGrad gradV(-1.0, feSpaceU, feSpaceP, {1}, feSpaceU.dof.size, 2*feSpaceU.dof.size);
  AssemblyDiv divV(-1.0, feSpaceP, feSpaceU, {1}, 2*feSpaceU.dof.size, feSpaceU.dof.size);

  Builder builder{numDOFs};
  // builder.assemblies[0].push_back(&stiffness0);
  // builder.assemblies[1].push_back(&grad0);
  // builder.assemblies[1].push_back(&div0);
  // builder.assemblies[0].push_back(&stiffness1);
  builder.buildProblem(stiffnessU, bcsU);
  builder.buildProblem(make_assemblyGrad(-1.0, feSpaceU, feSpaceP, {0}, 0, 2*feSpaceU.dof.size), bcsU, bcsP);
  builder.buildProblem(divU, bcsP, bcsU);
  builder.buildProblem(stiffnessV, bcsV);
  builder.buildProblem(gradV, bcsV, bcsP);
  builder.buildProblem(divV, bcsP, bcsV);
  builder.closeMatrix();

  Vec sol{numDOFs};
  // Eigen::GMRES<Mat> solver(builder.A);
  Eigen::UmfPackLU<Mat> solver(builder.A);
  // Eigen::SuperLU<Mat> solver(builder.A);
  sol = solver.solve(builder.b);

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  std::cout << sol.norm() << std::endl;

  Var u{"u", sol, 0, feSpaceU.dof.size};
  Var v{"v", sol, feSpaceU.dof.size, feSpaceU.dof.size};
  Var p{"p", sol, 2*feSpaceU.dof.size, feSpaceP.dof.size};

  // Var exact{"exact"};
  // interpolateAnalyticFunction(exact_sol, feSpaceU, exact.data);
  Var error{"e"};
  error.data = sol /*- exact.data*/;

  IOManager ioU{feSpaceU, "output_stokes2dcircle/vel"};
  ioU.print({u, v});
  IOManager ioP{feSpaceP, "output_stokes2dcircle/p"};
  ioP.print({p});

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  // if(std::fabs(norm - 6.03323) > 1.e-4)
  // {
  //   std::cerr << "the norm of the error is not the prescribed value" << std::endl;
  //   return 1;
  // }

  return 0;
}
