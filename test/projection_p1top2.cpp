#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpace1_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;
using FESpace2_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR>;

int main()
{
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, origin, length, {{2, 2, 0}});

  FESpace1_T feSpace1{*mesh};
  FESpace2_T feSpace2{*mesh};

  BCList<FESpace1_T> bc1{feSpace1};
  BCList<FESpace2_T> bc2{feSpace2};

  Var u1{"u1"};
  Var u2{"u2"};

  auto inputFun = [] (Vec3 const &p) {return 1.*p(0) + 2.*p(1);};
//  auto inputFun = [] (Vec3 const &) {return 1.;};
  interpolateAnalyticFunction(inputFun, feSpace1, u1.data);

  AssemblyMass<FESpace2_T> mass2(1.0, feSpace2);
  AssemblyProjection<FESpace2_T, FESpace1_T> mass1(1.0, u1.data, feSpace2, feSpace1);
  Builder builder{feSpace2.dof.size};
  builder.buildProblem(mass2, bc2);
  builder.buildProblem(mass1, bc2);
  builder.closeMatrix();
  std::cout << "A:\n" << builder.A << std::endl;
  std::cout << "b:\n" << builder.b << std::endl;

  LUSolver solver(builder.A);
  u2.data = solver.solve(builder.b);
  std::cout << "u1:\n" << u1.data << std::endl;
  std::cout << "u2:\n" << u2.data << std::endl;

  IOManager<FESpace1_T> io1{feSpace1, "projection1"};
  io1.print({u1});
  IOManager<FESpace2_T> io2{feSpace2, "projection2"};
  io2.print({u2});

  return 0;
}
