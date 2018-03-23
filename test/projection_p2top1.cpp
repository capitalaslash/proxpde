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
  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{2, 2, 0}});

  FESpace1_T feSpace1{meshPtr};
  FESpace2_T feSpace2{meshPtr};
  BCList<FESpace1_T> bc1{feSpace1};

  Var u1("u1", feSpace1.dof.totalNum);
  Var u2("u2", feSpace2.dof.totalNum);

  auto inputFun = [] (Vec3 const &p) {return 1.*p(0) + 2.*p(1);};
//  auto inputFun = [] (Vec3 const &) {return 1.;};
  interpolateAnalyticFunction(inputFun, feSpace2, u2.data);

  AssemblyMass<FESpace1_T> mass1(1.0, feSpace1);
  AssemblyProjection<FESpace1_T, FESpace2_T> proj2(1.0, u2.data, feSpace1, feSpace2);
  Builder builder{feSpace1.dof.totalNum};
  builder.buildProblem(mass1, bc1);
  builder.buildProblem(proj2, bc1);
  builder.closeMatrix();
  std::cout << "A:\n" << builder.A << std::endl;
  std::cout << "b:\n" << builder.b << std::endl;

  LUSolver solver(builder.A);
  u1.data = solver.solve(builder.b);
  std::cout << "u1:\n" << u1.data << std::endl;
  std::cout << "u2:\n" << u2.data << std::endl;

  IOManager<FESpace1_T> io1{feSpace1, "projection1"};
  io1.print({u1});
  IOManager<FESpace2_T> io2{feSpace2, "projection2"};
  io2.print({u2});

  return 0;
}