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
using FESpaceP1_T = FESpace<Mesh_T, RefTriangleP1, QR, 2>;
using FESpaceRT0_T = FESpace<Mesh_T, RefTriangleRT0, QR>;

int main()
{
  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};
  // Vec3 const origin{0., 0., 0.};
  // Vec3 const length{1., 1., 0.};
  // MeshBuilder<Elem_T> meshBuilder;
  // meshBuilder.build(meshPtr, origin, length, {{2, 2, 0}});
  meshPtr->pointList = {
    Point(Vec3(0., 0., 0.), 0),
    Point(Vec3(1., 0., 0.), 1),
    Point(Vec3(0., 1., 0.), 2)
  };
  meshPtr->elementList = {Triangle{{&meshPtr->pointList[0],
                                    &meshPtr->pointList[1],
                                    &meshPtr->pointList[2]},
                                    0}};
  meshPtr->buildConnectivity();
  // buildFacets(meshPtr);
  // markFacets2D(meshPtr, origin, length);

  FESpaceP1_T feSpaceP1{meshPtr};
  FESpaceRT0_T feSpaceRT0{meshPtr};

  BCList<FESpaceRT0_T> bc{feSpaceRT0};

  Var uP1{"uP1"};
  Var uRT0{"uRT0"};

  // auto inputFun = [] (Vec3 const &p) {return 1.*p(0) + 2.*p(1);};
  auto inputFun = [] (Vec3 const &) {return FVec<2>(1., 2.);};
  interpolateAnalyticFunction(inputFun, feSpaceP1, uP1.data);

  AssemblyVectorMass<FESpaceRT0_T> massRT0(1.0, feSpaceRT0);
  AssemblyVectorProjection<FESpaceRT0_T,FESpaceP1_T> massP1(1.0, uP1.data, feSpaceRT0, feSpaceP1);
  Builder builder{feSpaceRT0.dof.totalNum};
  builder.buildProblem(massRT0, bc);
  builder.buildProblem(massP1, bc);
  builder.closeMatrix();
  std::cout << "A:\n" << builder.A << std::endl;
  std::cout << "b:\n" << builder.b << std::endl;

  LUSolver solver(builder.A);
  uRT0.data = solver.solve(builder.b);
  std::cout << "uP1:\n" << uP1.data << std::endl;
  std::cout << "uRT0:\n" << uRT0.data << std::endl;

//  IOManager<FESpace1_T> io1{feSpace1, "projection1"};
//  io1.print({u1});
//  IOManager<FESpace2_T> io2{feSpace2, "projection2"};
//  io2.print({u2});

  return 0;
}
