#include "def.hpp"

#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

#include <limits>

using namespace proxpde;

static int test_id = 0;

template <typename Elem>
int test(
    std::array<Vec3, Elem::numPts> const & coords,
    std::function<FVec<Elem::dim>(Vec3 const &)> const & inputFun,
    long double const expectedError,
    long double const tolerance)
{
  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using QR = typename LagrangeFE<Elem_T, 1>::RecommendedQR;
  using FESpaceP0_T =
      FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, QR, Elem_T::dim>;
  using FESpaceP0Scalar_T =
      FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, QR>;
  using FESpaceP1_T =
      FESpace<Mesh_T, typename LagrangeFE<Elem_T, 1>::RefFE_T, QR, Elem_T::dim>;
  using FESpaceRT0_T =
      FESpace<Mesh_T, typename RaviartThomasFE<Elem_T, 0>::RefFE_T, QR>;

  // flux feSpace
  using Facet_T = typename Elem_T::Facet_T;
  using MeshFacet_T = Mesh<Facet_T>;
  using FESpaceFacet_T = FESpace<
      MeshFacet_T,
      typename LagrangeFE<Facet_T, 1>::RefFE_T,
      typename LagrangeFE<Facet_T, 1>::RecommendedQR,
      2>;

  std::cout << Utils::separator << "test " << test_id << std::endl;
  test_id++;

  std::unique_ptr<Mesh_T> oneElemMesh{new Mesh_T};
  referenceMesh(*oneElemMesh);

  // The new points MUST be well-sorted wrt to the reference element
  for (uint i = 0; i < Elem_T::numPts; ++i)
  {
    oneElemMesh->pointList[i].coord = coords[i];
  }
  buildNormals(*oneElemMesh);

  oneElemMesh->facetList[0].marker = side::BOTTOM;
  if (oneElemMesh->facetList.size() == 4)
    oneElemMesh->facetList[3].marker = side::LEFT;

  std::unique_ptr<MeshFacet_T> facetMesh{new MeshFacet_T};
  buildFacetMesh(*facetMesh, *oneElemMesh);

  FESpaceP1_T feSpaceP1{*oneElemMesh};
  Var uP1{"uP1"};
  // interpolateAnalyticFunction(inputFun, feSpace1, u1.data);
  projectAnalyticFunction(inputFun, feSpaceP1, uP1.data);
  // std::cout << "u1:\n" << u1.data << std::endl;

  FESpaceFacet_T feSpaceFacet{*facetMesh};
  Vec fluxes;
  interpolateOnFacets<Component::NORMAL>(fluxes, feSpaceFacet, uP1.data, feSpaceP1);
  std::cout << "fluxes1:\n" << fluxes << std::endl;

  FESpaceRT0_T feSpaceRT0{*oneElemMesh};
  Var uRT0{"uRT0"};
  l2Projection(uRT0.data, feSpaceRT0, uP1.data, feSpaceP1);
  std::cout << "uRT0:\n" << uRT0.data << std::endl;

  Var uRT0proj{"uRT0proj"};
  auto const inputFun3d = [&inputFun](Vec3 const & p)
  { return promote<3>(inputFun(p)); };
  projectAnalyticFunction(inputFun3d, feSpaceRT0, uRT0proj.data);
  std::cout << "uRT0proj:\n" << uRT0proj.data << std::endl;

  Var uRT0div{"uRT0div"};
  uint const sizeRT0 = feSpaceRT0.dof.size;
  FESpaceP0Scalar_T feSpaceLambda{*oneElemMesh, sizeRT0};
  uint const sizeP0 = feSpaceLambda.dof.size;
  // auto bcBottom = BCEss{feSpaceRT0, side::BOTTOM};
  // bcBottom << [] (Vec3 const & ) { return -0.5; };
  auto bcLeft = BCEss{feSpaceRT0, side::LEFT};
  bcLeft << [](Vec3 const &) { return Vec3{1.0, 0.0, 0.0}; };
  auto const bcs = std::tuple{bcLeft};
  Builder builder{sizeRT0 + sizeP0};
  builder.buildLhs(std::tuple{AssemblyVectorMass{1.0, feSpaceRT0}}, bcs);
  builder.buildCoupling(
      AssemblyVectorGrad(1.0, feSpaceRT0, feSpaceLambda), bcs, std::tuple{});
  builder.buildCoupling(
      AssemblyVectorDiv(1.0, feSpaceLambda, feSpaceRT0), std::tuple{}, bcs);
  builder.buildRhs(
      std::tuple{AssemblyProjection(1.0, uP1.data, feSpaceP1, feSpaceRT0)}, bcs);
  builder.closeMatrix();
  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  Vec const sol = solver.solve(builder.b);
  std::cout << "uRT0div (+lambda):\n" << sol << std::endl;
  uRT0div.data = sol.head(sizeRT0);
  Var lambda("lambda", sizeP0);
  lambda.data = sol.tail(sizeP0);

  Var exactRT0{"exactRT0"};
  interpolateAnalyticFunction(inputFun3d, feSpaceRT0, exactRT0.data);
  std::cout << "exactRT0:\n" << exactRT0.data << std::endl;

  Var uP0{"uP0"};
  FESpaceP0_T feSpaceP0{*oneElemMesh};
  l2Projection(uP0.data, feSpaceP0, uRT0.data, feSpaceRT0);
  // std::cout << "u0:\n" << u0.data << std::endl;

  Vec exact0;
  integrateAnalyticFunction(inputFun, feSpaceP0, exact0);
  double const error = (uP0.data - exact0).norm();

  IOManager io1{feSpaceP1, "output_projrt0/u1"};
  io1.print({uP1});
  IOManager io0{feSpaceP0, "output_projrt0/u0"};
  io0.print({uP0});

  return checkError({error}, {expectedError}, tolerance);
}

int main()
{
  std::bitset<16> tests;

  auto const constFun = [](Vec3 const &) { return Vec3{1.0, 2.0, 3.0}; };
  auto const linearFun = [](Vec3 const & /*p*/)
  {
    // return Vec3{
    //   3. + p[0],
    //   1. + 2. * p[0] - p[1],
    //   2. + 3. * p[0] + 4. * p[1] + 5. * p[2]};
    return Vec3{0., -1., 0.};
  };

  auto const eps = std::numeric_limits<double>::epsilon();

  // Triangle

  std::array<Vec3, 3> const refTrianglePts = {
      Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};

  tests[0] = test<Triangle>(
      refTrianglePts,
      [&constFun](Vec3 const & p) { return narrow<2>(constFun(p)); },
      eps,
      eps);

  tests[1] = test<Triangle>(
      refTrianglePts,
      [&linearFun](Vec3 const & p) { return narrow<2>(linearFun(p)); },
      2 * eps,
      eps);

  std::array<Vec3, 3> const trianglePts = {
      Vec3{1.0, 0.0, 0.0}, Vec3{3.0, 1.0, 0.0}, Vec3{0.0, 5.0, 0.0}};

  tests[2] = test<Triangle>(
      trianglePts,
      [&constFun](Vec3 const & p) { return narrow<2>(constFun(p)); },
      4 * eps,
      eps);

  tests[3] = test<Triangle>(
      trianglePts,
      [&linearFun](Vec3 const & p) { return narrow<2>(linearFun(p)); },
      8 * eps,
      eps);

  // Quad

  std::array<Vec3, 4> const refQuadPts = {
      Vec3{-1.0, -1.0, 0.0},
      Vec3{1.0, -1.0, 0.0},
      Vec3{1.0, 1.0, 0.0},
      Vec3{-1.0, 1.0, 0.0}};

  tests[4] = test<Quad>(
      refQuadPts,
      [&constFun](Vec3 const & p) { return narrow<2>(constFun(p)); },
      2 * eps,
      eps);

  tests[5] = test<Quad>(
      refQuadPts,
      [&linearFun](Vec3 const & p) { return narrow<2>(linearFun(p)); },
      3 * eps,
      eps);

  // std::array<Vec3, 4> const quadPts = {
  //   Vec3{0.0, 0.0, 0.0},
  //   Vec3{3.0, 1.0, 0.0},
  //   Vec3{2.0, 5.0, 0.0},
  //   Vec3{0.0, 5.0, 0.0},
  // };

  // stretched and roneElemMesh->facetList[3]otated
  // std::array<Vec3, 4> const quadPts = {
  //   Vec3{0.0, 0.0, 0.0},
  //   Vec3{4.0, 2.0, 0.0},
  //   Vec3{3.0, 4.0, 0.0},
  //   Vec3{-1.0, 2.0, 0.0},
  // };

  // non-affine
  std::array<Vec3, 4> const quadPts = {
      Vec3{-1.0, -1.0, 0.0},
      Vec3{1.0, -1.0, 0.0},
      Vec3{1.0, 0.5, 0.0},
      Vec3{-1.0, 1.0, 0.0},
  };

  tests[6] = test<Quad>(
      quadPts,
      [&constFun](Vec3 const & p) { return narrow<2>(constFun(p)); },
      7.021666937153402e-16,
      eps);

  tests[7] = test<Quad>(
      quadPts,
      [&linearFun](Vec3 const & p) { return narrow<2>(linearFun(p)); },
      9.085915191948722e-02,
      eps);

  //  // Tetrahedron

  //  std::array<Vec3, 4> const refTetPts = {
  //    Vec3{0.0, 0.0, 0.0},
  //    Vec3{1.0, 0.0, 0.0},
  //    Vec3{0.0, 1.0, 0.0},
  //    Vec3{0.0, 0.0, 1.0}
  //  };

  //  tests[8] = test<Tetrahedron>(refTetPts,
  //                               constFun,
  //                               5 * eps,
  //                               eps);

  //  tests[9] = test<Tetrahedron>(refTetPts,
  //                               linearFun,
  //                               13 * eps,
  //                               eps);

  //  std::array<Vec3, 4> const tetPts = {
  //    Vec3{0.0, 0.0, 0.0},
  //    Vec3{2.0, 0.0, 1.0},
  //    Vec3{1.0, 3.0, 2.0},
  //    Vec3{1.0, 1.0, 5.0}
  //  };

  //  tests[10] = test<Tetrahedron>(tetPts,
  //                                constFun,
  //                                0.,
  //                                eps);

  //  tests[11] = test<Tetrahedron>(tetPts,
  //                                linearFun,
  //                                20 * eps,
  //                                eps);

  //  // Hexahedron

  //  std::array<Vec3, 8> const refHexaPts = {
  //    Vec3{-1.0, -1.0, -1.0},
  //    Vec3{ 1.0, -1.0, -1.0},
  //    Vec3{ 1.0,  1.0, -1.0},
  //    Vec3{-1.0,  1.0, -1.0},
  //    Vec3{-1.0, -1.0,  1.0},
  //    Vec3{ 1.0, -1.0,  1.0},
  //    Vec3{ 1.0,  1.0,  1.0},
  //    Vec3{-1.0,  1.0,  1.0},
  //  };

  //  tests[12] = test<Hexahedron>(refHexaPts,
  //                               constFun,
  //                               2 * eps,
  //                               eps);

  //  tests[13] = test<Hexahedron>(refHexaPts,
  //                               linearFun,
  //                               2 * eps,
  //                               eps);

  //  // auto const theta = M_PI / 3.;
  //  // auto const axis = Vec3{1.0, 2.0, 3.0}.normalized();
  //  // auto const scale = Vec3{1.5, 2.0, 0.4};
  //  // auto const translate = Vec3{2.0, 3.0, 1.0};

  //  // // order of operation changes the result!
  //  // auto const threeLinearTransform = [&axis, &theta, &scale, &translate] (Vec3
  //  const & p)
  //  // {
  //  //   auto const rotMatrix = (FMat<3,3>() <<
  //  //       std::cos(theta) + axis[0] * axis[0] * (1. - std::cos(theta)),
  //  //       axis[0] * axis[1] * (1. - std::cos(theta)) - axis[2] * std::sin(theta),
  //  //       axis[0] * axis[2] * (1. - std::cos(theta)) + axis[1] * std::sin(theta),
  //  //       axis[1] * axis[0] * (1. - std::cos(theta)) + axis[2] * std::sin(theta),
  //  //       std::cos(theta) + axis[1] * axis[1] * (1. - std::cos(theta)),
  //  //       axis[1] * axis[2] * (1. - std::cos(theta)) - axis[0] * std::sin(theta),
  //  //       axis[2] * axis[0] * (1. - std::cos(theta)) - axis[1] * std::sin(theta),
  //  //       axis[2] * axis[1] * (1. - std::cos(theta)) + axis[0] * std::sin(theta),
  //  //       std::cos(theta) + axis[2] * axis[2] * (1. - std::cos(theta))).finished();

  //  //   auto const scaleMatrix = (FMat<3,3>{} <<
  //  //                             scale[0], 0.0, 0.0,
  //  //                             0.0, scale[1], 0.0,
  //  //                             0.0, 0.0, scale[2]).finished();

  //  //   Vec3 const pOut = scaleMatrix * rotMatrix * (p + translate);
  //  //   return pOut;
  //  // };

  //  // // affine
  //  // std::array<Vec3, 8> const hexaPts = {
  //  //   threeLinearTransform(Vec3{-1.0, -1.0, -1.0}),
  //  //   threeLinearTransform(Vec3{ 1.0, -1.0, -1.0}),
  //  //   threeLinearTransform(Vec3{ 1.0,  1.0, -1.0}),
  //  //   threeLinearTransform(Vec3{-1.0,  1.0, -1.0}),
  //  //   threeLinearTransform(Vec3{-1.0, -1.0,  1.0}),
  //  //   threeLinearTransform(Vec3{ 1.0, -1.0,  1.0}),
  //  //   threeLinearTransform(Vec3{ 1.0,  1.0,  1.0}),
  //  //   threeLinearTransform(Vec3{-1.0,  1.0,  1.0}),
  //  // };

  //  // non-affine
  //  std::array<Vec3, 8> const hexaPts = {
  //    Vec3{-1.0, -1.0, -1.0},
  //    Vec3{ 1.0, -1.0, -1.0},
  //    Vec3{ 1.0,  1.0, -1.0},
  //    Vec3{-1.0,  1.0, -1.0},
  //    Vec3{-1.0, -1.0,  1.0},
  //    Vec3{ 1.0, -1.0,  1.0},
  //    Vec3{ 1.0,  1.0,  0.5},
  //    Vec3{-1.0,  1.0,  0.5},
  //  };

  //  tests[14] = test<Hexahedron>(hexaPts,
  //                               constFun,
  //                               6.802721088435382e-03,
  //                               eps);

  //  tests[15] = test<Hexahedron>(hexaPts,
  //                               linearFun,
  //                               2.159947638851983e-01,
  //                               eps);

  std::cout << "test results: " << tests << std::endl;
  return tests.any();
}
