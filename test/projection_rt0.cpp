#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "feutils.hpp"

template <typename Elem>
int test(
    std::array<Vec3, Elem::numPts> const & coords,
    std::function<FVec<Elem::dim> (Vec3 const &)> const & inputFun,
    long double const expectedError,
    long double const tolerance)
{
  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using QR = typename LagrangeFE<Elem_T, 1>::RecommendedQR;
  using FESpaceP0_T = FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, QR, Elem_T::dim>;
  using FESpaceP1_T = FESpace<Mesh_T, typename LagrangeFE<Elem_T, 1>::RefFE_T, QR, Elem_T::dim>;
  using FESpaceRT0_T = FESpace<Mesh_T, typename RaviartThomasFE<Elem_T, 0>::RefFE_T, QR>;

  std::unique_ptr<Mesh_T> oneElemMesh{new Mesh_T};
  referenceMesh(*oneElemMesh);

  // The new points MUST be well-sorted wrt to the reference element
  for (uint i = 0; i < Elem_T::numPts; ++i)
  {
    oneElemMesh->pointList[i].coord = coords[i];
  }

  FESpaceP0_T feSpace0{*oneElemMesh};
  FESpaceP1_T feSpace1{*oneElemMesh};
  FESpaceRT0_T feSpaceRT0{*oneElemMesh};

  Var u0{"u0"};
  Var u1{"u1"};
  Var uRT0{"uRT0"};

  interpolateAnalyticFunction(inputFun, feSpace1, u1.data);
  // std::cout << "u1:\n" << u1.data << std::endl;

  l2Projection(uRT0.data, feSpaceRT0, u1.data, feSpace1);
  // std::cout << "uRT0:\n" << uRT0.data << std::endl;

  l2Projection(u0.data, feSpace0, uRT0.data, feSpaceRT0);
  // std::cout << "u0:\n" << u0.data << std::endl;

  Vec exact0;
  integrateAnalyticFunction(inputFun, feSpace0, exact0);
  double const error = (u0.data - exact0).norm();

  IOManager io1{feSpace1, "output_projrt0/u1"};
  io1.print({u1});
  IOManager io0{feSpace0, "output_projrt0/u0"};
  io0.print({u0});

  return checkError({error}, {expectedError}, tolerance);
}

int main()
{
  std::bitset<16> tests;

  auto const constFun = [] (Vec3 const & ) { return Vec3{1.0, 2.0, 3.0}; };
  auto const linearFun = [] (Vec3 const & p) {
    return Vec3{
      3. + p[0],
      1. + 2. * p[0] + 3. * p[1],
      2. + 3. * p[0] + 4. * p[1] + 5. * p[2]};
  };

  // Triangle

  std::array<Vec3, 3> const refTrianglePts = {
    Vec3{0.0, 0.0, 0.0},
    Vec3{1.0, 0.0, 0.0},
    Vec3{0.0, 1.0, 0.0}
  };

  tests[0] = test<Triangle>(refTrianglePts,
                            [&constFun] (Vec3 const & p) { return narrow<2>(constFun(p)); },
                            2.22045e-16L,
                            1.e-17L);

  tests[1] = test<Triangle>(refTrianglePts,
                            [&linearFun] (Vec3 const & p) { return narrow<2>(linearFun(p)); },
                            9.93014e-16L,
                            1.e-17L);

  std::array<Vec3, 3> const trianglePts = {
    Vec3{1.0, 0.0, 0.0},
    Vec3{3.0, 1.0, 0.0},
    Vec3{0.0, 5.0, 0.0}
  };

  tests[2] = test<Triangle>(trianglePts,
                            [&constFun] (Vec3 const & p) { return narrow<2>(constFun(p)); },
                            9.93014e-16L,
                            1.e-17L);

  tests[3] = test<Triangle>(trianglePts,
                            [&linearFun] (Vec3 const & p) { return narrow<2>(linearFun(p)); },
                            2.51215e-15L,
                            1.e-17);

  // Quad

  std::array<Vec3, 4> const refQuadPts = {
    Vec3{-1.0, -1.0, 0.0},
    Vec3{ 1.0, -1.0, 0.0},
    Vec3{ 1.0,  1.0, 0.0},
    Vec3{-1.0,  1.0, 0.0}
  };

  tests[4] = test<Quad>(refQuadPts,
                        [&constFun] (Vec3 const & p) { return narrow<2>(constFun(p)); },
                        4.57757e-16L,
                        1.e-17L);

  tests[5] = test<Quad>(refQuadPts,
                        [&linearFun] (Vec3 const & p) { return narrow<2>(linearFun(p)); },
                        4.44089e-16L,
                        1.e-17L);

  std::array<Vec3, 4> const quadPts = {
    Vec3{0.0, 0.0, 0.0},
    Vec3{3.0, 1.0, 0.0},
    Vec3{2.0, 5.0, 0.0},
    Vec3{0.0, 5.0, 0.0}
  };

  tests[6] = test<Quad>(quadPts,
                        [&constFun] (Vec3 const & p) { return narrow<2>(constFun(p)); },
                        2.48253e-16L,
                        1.e-17L);

  tests[7] = test<Quad>(quadPts,
                        [&linearFun] (Vec3 const & p) { return narrow<2>(linearFun(p)); },
                        1.77636e-15L,
                        1.e-17L);

  // Tetrahedron

  std::array<Vec3, 4> const refTetPts = {
    Vec3{0.0, 0.0, 0.0},
    Vec3{1.0, 0.0, 0.0},
    Vec3{0.0, 1.0, 0.0},
    Vec3{0.0, 0.0, 1.0}
  };

  tests[8] = test<Tetrahedron>(refTetPts,
                               constFun,
                               1.04738e-15L,
                               1.e-17L);

  tests[9] = test<Tetrahedron>(refTetPts,
                               linearFun,
                               2.84356e-15L,
                               1.e-17L);

  std::array<Vec3, 4> const tetPts = {
    Vec3{0.0, 0.0, 0.0},
    Vec3{2.0, 0.0, 1.0},
    Vec3{1.0, 3.0, 2.0},
    Vec3{1.0, 1.0, 5.0}
  };

  tests[10] = test<Tetrahedron>(tetPts,
                                constFun,
                                0.0L,
                                1.e-17L);

  tests[11] = test<Tetrahedron>(tetPts,
                                linearFun,
                                4.35117e-15L,
                                1.e-17L);

  // Hexahedron

  std::array<Vec3, 8> const refHexaPts = {
    Vec3{-1.0, -1.0, -1.0},
    Vec3{ 1.0, -1.0, -1.0},
    Vec3{ 1.0,  1.0, -1.0},
    Vec3{-1.0,  1.0, -1.0},
    Vec3{-1.0, -1.0,  1.0},
    Vec3{ 1.0, -1.0,  1.0},
    Vec3{ 1.0,  1.0,  1.0},
    Vec3{-1.0,  1.0,  1.0},
  };

  tests[12] = test<Hexahedron>(refHexaPts,
                               constFun,
                               1.48952e-15L,
                               1.e-17L);

  tests[13] = test<Hexahedron>(refHexaPts,
                               linearFun,
                               2.38375e-15L,
                               1.e-17L);

  auto const theta = M_PI / 3.;
  auto const axis = Vec3{1.0, 2.0, 3.0}.normalized();
  auto const scale = Vec3{1.5, 2.0, 0.4};
  auto const translate = Vec3{2.0, 3.0, 1.0};

  // order of operation changes the result!
  auto const threeLinearTransform = [&axis, &theta, &scale, &translate] (Vec3 const & p)
  {
    auto const rotMatrix = (FMat<3,3>() <<
        std::cos(theta) + axis[0] * axis[0] * (1. - std::cos(theta)),
        axis[0] * axis[1] * (1. - std::cos(theta)) - axis[2] * std::sin(theta),
        axis[0] * axis[2] * (1. - std::cos(theta)) + axis[1] * std::sin(theta),
        axis[1] * axis[0] * (1. - std::cos(theta)) + axis[2] * std::sin(theta),
        std::cos(theta) + axis[1] * axis[1] * (1. - std::cos(theta)),
        axis[1] * axis[2] * (1. - std::cos(theta)) - axis[0] * std::sin(theta),
        axis[2] * axis[0] * (1. - std::cos(theta)) - axis[1] * std::sin(theta),
        axis[2] * axis[1] * (1. - std::cos(theta)) + axis[0] * std::sin(theta),
        std::cos(theta) + axis[2] * axis[2] * (1. - std::cos(theta))).finished();

    auto const scaleMatrix = (FMat<3,3>{} <<
                              scale[0], 0.0, 0.0,
                              0.0, scale[1], 0.0,
                              0.0, 0.0, scale[2]).finished();

    Vec3 const pOut = scaleMatrix * rotMatrix * (p + translate);
    return pOut;
  };

  std::array<Vec3, 8> const hexaPts = {
    threeLinearTransform(Vec3{-1.0, -1.0, -1.0}),
    threeLinearTransform(Vec3{ 1.0, -1.0, -1.0}),
    threeLinearTransform(Vec3{ 1.0,  1.0, -1.0}),
    threeLinearTransform(Vec3{-1.0,  1.0, -1.0}),
    threeLinearTransform(Vec3{-1.0, -1.0,  1.0}),
    threeLinearTransform(Vec3{ 1.0, -1.0,  1.0}),
    threeLinearTransform(Vec3{ 1.0,  1.0,  1.0}),
    threeLinearTransform(Vec3{-1.0,  1.0,  1.0}),
  };

  tests[14] = test<Hexahedron>(hexaPts,
                               constFun,
                               2.09476e-15L,
                               1.e-17L);

  tests[15] = test<Hexahedron>(hexaPts,
                               linearFun,
                               1.3293e-14L,
                               1.e-17L);

  std::cout << "test results: " << tests << std::endl;
  return tests.any();
}
