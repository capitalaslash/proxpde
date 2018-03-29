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

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;

int main(int argc, char* argv[])
{
//  uint const numPts_x = (argc < 3)? 3 : std::stoi(argv[1]);
//  uint const numPts_y = (argc < 3)? 3 : std::stoi(argv[2]);

//  Vec3 const origin{0., 0., 0.};
//  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};

//  MeshBuilder<Elem_T> meshBuilder;
//  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});
  readGMSH(meshPtr, "square_uns.msh");

  FESpaceVel_T feSpaceVel{meshPtr};
  FESpaceP_T feSpaceP{meshPtr};
  // std::cout << feSpaceVel.dof << std::endl;

  auto feList = std::make_tuple(feSpaceVel, feSpaceP);
  auto assembler = make_assembler(feList);
  // auto assembler = make_assembler(std::forward_as_tuple(feSpaceU, feSpaceU, feSpaceP));

  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  auto inlet = [] (Vec3 const & p) {return Vec2(0., 0.5*(1.-p(0)*p(0)));};
  BCList<FESpaceVel_T> bcsVel{feSpaceVel};
  bcsVel.addEssentialBC(side::BOTTOM, inlet);
  bcsVel.addEssentialBC(side::RIGHT, zero);
  bcsVel.addEssentialBC(side::TOP, zero, {0});
  bcsVel.addEssentialBC(side::LEFT, zero, {0});
  // bcsVel.addNaturalBC(side::BOTTOM, [] (Point const &) {return Vec2(0.0, 1.0);});
  BCList<FESpaceP_T> bcsP{feSpaceP};

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU*FESpaceVel_T::dim + dofP;

  AssemblyTensorStiffness<FESpaceVel_T> stiffness(1.0, feSpaceVel);
  AssemblyGrad<FESpaceVel_T, FESpaceP_T> grad(feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv<FESpaceP_T, FESpaceVel_T> div(feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);

  Builder builder{numDOFs};
  // builder.assemblies[0].push_back(&stiffness0);
  // builder.assemblies[1].push_back(&grad0);
  // builder.assemblies[1].push_back(&div0);
  // builder.assemblies[0].push_back(&stiffness1);
  builder.buildProblem(stiffness, bcsVel);
  builder.buildProblem(grad, bcsVel, bcsP);
  builder.buildProblem(div, bcsP, bcsVel);
  builder.closeMatrix();

  Var sol("vel", numDOFs);
  Eigen::UmfPackLU<Mat> solver(builder.A);
  sol.data = solver.solve(builder.b);

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "builder.b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  Var exact{"exact", numDOFs};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact.data);
  interpolateAnalyticFunction([](Vec3 const & p){return p(1)-1;}, feSpaceP, exact.data, 2*dofU);

  std::cout << sol.data.norm() << std::endl;

  Var u{"u", sol.data, 0, dofU};
  Var v{"v", sol.data, dofU, dofU};
  Var p{"p", sol.data, 2*dofU, dofP};
  Var ue{"ue", exact.data, 0, dofU};
  Var ve{"ve", exact.data, dofU, dofU};
  Var pe{"pe", exact.data, 2*dofU, dofP};

  IOManager<FESpaceVel_T> ioVel{feSpaceVel, "sol_stokes2dtri_vel"};
  ioVel.print({sol, exact});
  IOManager<FESpaceP_T> ioP{feSpaceP, "sol_stokes2dtri_p"};
  ioP.print({p, pe});

  auto uNorm = (u.data - ue.data).norm();
  auto vNorm = (v.data - ve.data).norm();
  auto pNorm = (p.data - pe.data).norm();

  std::cout << "u error norm: " << uNorm << std::endl;
  std::cout << "v error norm: " << vNorm << std::endl;
  std::cout << "p error norm: " << pNorm << std::endl;

  if (uNorm + vNorm + pNorm > 2.e-13)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
