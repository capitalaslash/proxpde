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
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;

int main(int argc, char* argv[])
{
  uint const numPts_x = (argc < 3)? 4 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 4 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});

  FESpaceVel_T feSpaceVel{meshPtr};
  FESpaceP_T feSpaceP{meshPtr};

  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  BCList<FESpaceVel_T> bcsVel{feSpaceVel};
  bcsVel.addEssentialBC(side::RIGHT, zero);
  bcsVel.addEssentialBC(side::LEFT, zero);
  bcsVel.addEssentialBC(side::BOTTOM, zero);
  bcsVel.addEssentialBC(side::TOP, [] (Vec3 const &) {return Vec2(1.0, 0.0);});
  BCList<FESpaceP_T> bcsP{feSpaceP};
  // DofSet_T pinSet = {1};
  // bcsP.addEssentialBC(pinSet, [] (Vec3 const &) {return 0.;});

  auto const dofU = feSpaceVel.dof.totalNum;
  auto const dofP = feSpaceP.dof.totalNum;
  uint const numDOFs = dofU*FESpaceVel_T::dim + dofP;

  AssemblyStiffness<FESpaceVel_T> stiffness(1.0, feSpaceVel);
  AssemblyGrad<FESpaceVel_T, FESpaceP_T> grad(feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv<FESpaceP_T, FESpaceVel_T> div(feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);

  Builder builder{numDOFs};
  Var sol{"sol"};
  builder.buildProblem(stiffness, bcsVel);
  builder.buildProblem(grad, bcsVel, bcsP);
  builder.buildProblem(div, bcsP, bcsVel);
  builder.closeMatrix();

  Eigen::UmfPackLU<Mat> solver(builder.A);
  sol.data = solver.solve(builder.b);

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol << std::endl;

  Var u{"u", sol.data, 0, dofU};
  Var v{"v", sol.data, dofU, dofU};
  Var p{"p", sol.data, 2*dofU, dofP};

  auto uNorm = u.data.norm();
  auto vNorm = v.data.norm();
  auto pNorm = p.data.norm();

  std::cout << "norm of u: " << uNorm << std::endl;
  std::cout << "norm of v: " << vNorm << std::endl;
  std::cout << "norm of p: " << pNorm << std::endl;

  IOManager<FESpaceVel_T> ioVel{feSpaceVel, "sol_cavity_vel"};
  ioVel.print({sol});
  IOManager<FESpaceP_T> ioP{feSpaceP, "sol_cavity_p"};
  ioP.print({p});

  if(std::fabs(uNorm - 2.72045) > 1.e-4 ||
     std::fabs(vNorm - 0.47505) > 1.e-4 ||
     std::fabs(pNorm - 52.9371) > 1.e-4)
  {
    std::cerr << "one of the norms is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
