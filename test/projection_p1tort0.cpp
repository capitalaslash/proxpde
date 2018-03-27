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
  std::shared_ptr<Mesh_T> trianglePtr{new Mesh_T};
  trianglePtr->pointList = {
    Point(Vec3(0., 0., 0.), 0),
    Point(Vec3(1., 0., 0.), 1),
    Point(Vec3(0., 1., 0.), 2)
  };
  trianglePtr->elementList = {Triangle{{&trianglePtr->pointList[0],
                                    &trianglePtr->pointList[1],
                                    &trianglePtr->pointList[2]},
                                    0}};
  trianglePtr->buildConnectivity();

  FESpaceP0_T feSpaceP0{trianglePtr};
  FESpaceP1_T feSpaceP1{trianglePtr};
  FESpaceRT0_T feSpaceRT0{trianglePtr};

  BCList<FESpaceRT0_T> bcRT0{feSpaceRT0};
  BCList<FESpaceP0_T> bcP0{feSpaceP0};

  Var uP0{"uP0"};
  Var uP1{"uP1"};
  Var uRT0{"uRT0"};

  // auto inputFun = [] (Vec3 const &p) {return 1.*p(0) + 2.*p(1);};
  auto inputFun = [] (Vec3 const &) {return FVec<2>(1., 2.);};
  interpolateAnalyticFunction(inputFun, feSpaceP1, uP1.data);

  AssemblyVectorMass<FESpaceRT0_T> massRT0(1.0, feSpaceRT0);
  AssemblyVectorProjection<FESpaceRT0_T,FESpaceP1_T> projP1RT0(1.0, uP1.data, feSpaceRT0, feSpaceP1);
  Builder builderRT0{feSpaceRT0.dof.totalNum};
  builderRT0.buildProblem(massRT0, bcRT0);
  builderRT0.buildProblem(projP1RT0, bcRT0);
  builderRT0.closeMatrix();
  std::cout << "A:\n" << builderRT0.A << std::endl;
  std::cout << "b:\n" << builderRT0.b << std::endl;

  LUSolver solver(builderRT0.A);
  uRT0.data = solver.solve(builderRT0.b);
  std::cout << "uP1:\n" << uP1.data << std::endl;
  std::cout << "uRT0:\n" << uRT0.data << std::endl;

  AssemblyMass<FESpaceP0_T> massP0(1.0, feSpaceP0);
  AssemblyProjection<FESpaceP0_T,FESpaceP1_T> projP1P0(1.0, uP1.data, feSpaceP0, feSpaceP1);
  Builder builderP0{feSpaceP0.dof.totalNum * feSpaceP0.dim};
  builderP0.buildProblem(massP0, bcP0);
  builderP0.buildProblem(projP1P0, bcP0);
  builderP0.closeMatrix();
  std::cout << "A:\n" << builderP0.A << std::endl;
  std::cout << "b:\n" << builderP0.b << std::endl;
  LUSolver solverP0(builderP0.A);
  uP0.data = solverP0.solve(builderP0.b);
  std::cout << "uP0:\n" << uP0.data << std::endl;



  //  IOManager<FESpace1_T> io1{feSpace1, "projection1"};
//  io1.print({u1});
//  IOManager<FESpace2_T> io2{feSpace2, "projection2"};
//  io2.print({u2});

  return 0;
}
