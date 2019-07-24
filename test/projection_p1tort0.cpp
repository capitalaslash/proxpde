#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "feutils.hpp"

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

  Var uP0{"uP0"};
  Var uP1{"uP1"};
  Var uRT0{"uRT0"};

  double errorRefElement;
  {
    auto inputFun = [] (Vec3 const & p) { return FVec<2>(1.*p(0) + 2.*p(1), 3.*p(0) - 4.*p(1)); };
    // auto inputFun = [] (Vec3 const &) { return FVec<2>{1., 2.}; };
    interpolateAnalyticFunction(inputFun, feSpaceP1, uP1.data);
    std::cout << "uP1:\n" << uP1.data << std::endl;

    l2Projection(uRT0.data, feSpaceRT0, uP1.data, feSpaceP1);
    std::cout << "uRT0:\n" << uRT0.data << std::endl;

    l2Projection(uP0.data, feSpaceP0, uRT0.data, feSpaceRT0);
    std::cout << "uP0:\n" << uP0.data << std::endl;

    errorRefElement = (uP0.data - Vec2{1, -1./3}).norm();
  }

  double errorRealElement;
  {
    triangleMesh->pointList[0].coord = Vec3{1., 0., 0.};
    triangleMesh->pointList[1].coord = Vec3{3., 1., 0.};
    triangleMesh->pointList[2].coord = Vec3{0., 5., 0.};

    auto inputFun = [] (Vec3 const & p) {return FVec<2>{1., 2.}; };

    interpolateAnalyticFunction(inputFun, feSpaceP1, uP1.data);
    std::cout << "uP1:\n" << uP1.data << std::endl;

    l2Projection(uRT0.data, feSpaceRT0, uP1.data, feSpaceP1);
    std::cout << "uRT0:\n" << uRT0.data << std::endl;

    l2Projection(uP0.data, feSpaceP0, uRT0.data, feSpaceRT0);
    std::cout << "uP0:\n" << uP0.data << std::endl;

    errorRealElement = (uP0.data - Vec2{1., 2.}).norm();
  }

  IOManager ioP1{feSpaceP1, "output_projrt0/p1"};
  ioP1.print({uP1});
  IOManager ioP0{feSpaceP0, "output_projrt0/p0"};
  ioP0.print({uP0});

  return checkError({errorRefElement, errorRealElement}, {2.e-16, 2.e-16}, 1.e-16);
}
