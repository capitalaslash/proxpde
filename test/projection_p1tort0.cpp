#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using P1RefFE = RefTriangleP1;
using RT0RefFE = RefTriangleRT0;
using QR = FEType<Elem_T, 1>::RecommendedQR;
using FESpaceP0_T = FESpace<Mesh_T, RefTriangleP0, QR, 2>;
using FESpaceP1_T = FESpace<Mesh_T, RefTriangleP1, QR, 2>;
using FESpaceRT0_T = FESpace<Mesh_T, RefTriangleRT0, QR>;

int main()
{
  std::unique_ptr<Mesh_T> triangleMesh{new Mesh_T};
  refTriangleMesh(*triangleMesh);
  addElemFacetList(*triangleMesh);

  FESpaceP0_T feSpaceP0{*triangleMesh};
  FESpaceP1_T feSpaceP1{*triangleMesh};
  FESpaceRT0_T feSpaceRT0{*triangleMesh};

  std::cout << "dofP0:\n" << feSpaceP0.dof << std::endl;

  BCList bcRT0{feSpaceRT0};
  BCList bcP0{feSpaceP0};

  Var uP0{"uP0"};
  Var uP1{"uP1"};
  Var uRT0{"uRT0"};

  auto inputFun = [] (Vec3 const &p) {return FVec<2>(1.*p(0) + 2.*p(1), 3.*p(0) - 4.*p(1));};
  // auto inputFun = [] (Vec3 const &) {return FVec<2>(1., 2.);};
  interpolateAnalyticFunction(inputFun, feSpaceP1, uP1.data);

  AssemblyVectorMass massRT0(1.0, feSpaceRT0);
  AssemblyProjection projP1RT0(1.0, uP1.data, feSpaceRT0, feSpaceP1);
  Builder builderRT0{feSpaceRT0.dof.size};
  builderRT0.buildProblem(massRT0, bcRT0);
  builderRT0.buildProblem(projP1RT0, bcRT0);
  builderRT0.closeMatrix();
  // std::cout << "A:\n" << builderRT0.A << std::endl;
  // std::cout << "b:\n" << builderRT0.b << std::endl;

  LUSolver solver(builderRT0.A);
  uRT0.data = solver.solve(builderRT0.b);
  std::cout << "uP1:\n" << uP1.data << std::endl;
  std::cout << "uRT0:\n" << uRT0.data << std::endl;

  AssemblyMass massP0(1.0, feSpaceP0);
  AssemblyProjection projRT0P0(1.0, uRT0.data, feSpaceP0, feSpaceRT0);
  Builder builderP0{feSpaceP0.dof.size * feSpaceP0.dim};
  builderP0.buildProblem(massP0, bcP0);
  builderP0.buildProblem(projRT0P0, bcP0);
  builderP0.closeMatrix();
  // std::cout << "A:\n" << builderP0.A << std::endl;
  // std::cout << "b:\n" << builderP0.b << std::endl;
  LUSolver solverP0(builderP0.A);
  uP0.data = solverP0.solve(builderP0.b);
  std::cout << "uP0:\n" << uP0.data << std::endl;

  IOManager ioP1{feSpaceP1, "output_projrt0/p1"};
  ioP1.print({uP1});
  IOManager ioP0{feSpaceP0, "output_projrt0/p0"};
  ioP0.print({uP0});

  return 0;
}
