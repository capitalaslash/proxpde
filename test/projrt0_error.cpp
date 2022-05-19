#include "def.hpp"

#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"

#include <limits>

template <typename FESpace>
double computeErrorL2(
    Fun<FESpace::physicalDim(), 3> const & f, Vec const & u, FESpace const & feSpace)
{
  double error = 0.0;
  for (auto const & elem: feSpace.mesh->elementList)
  {
    feSpace.curFE.reinit(elem);

    for (uint d = 0; d < FESpace::dim; ++d)
    {
      FMat<FESpace::CurFE_T::numDOFs, 1> uLocal;
      for (uint n = 0; n < FESpace::CurFE_T::numDOFs; ++n)
      {
        id_T const dofId = feSpace.dof.getId(elem.id, n, d);
        uLocal[n] = u[dofId];
      }

      for (uint q = 0; q < FESpace::QR_T::numPts; ++q)
      {
        FVec<FESpace::physicalDim()> uQ;
        if constexpr (family_v<typename FESpace::RefFE_T> == FamilyType::LAGRANGE)
        {
          uQ = uLocal.transpose() * feSpace.curFE.phi[q];
        }
        else // FamilyType::RAVIART_THOMAS
        {
          uQ = uLocal.transpose() * feSpace.curFE.phiVect[q];
        }
        FVec<FESpace::physicalDim()> const fQ = f(feSpace.curFE.qpoint[q]);
        error += feSpace.curFE.JxW[q] * (fQ - uQ).transpose() * (fQ - uQ);
      }
    }
  }
  return std::sqrt(error);
}

template <typename FESpace>
double computeErrorL2(scalarFun_T const & exact, Vec const & u, FESpace const & feSpace)
{
  static_assert(FESpace::dim == 1);
  return computeErrorL2(
      [&exact](Vec3 const & p) { return Vec1{exact(p)}; }, u, feSpace);
}

static int test_id = 0;

template <typename Elem>
int test(
    Mesh<Elem> const & meshCoarse,
    uint const refinements,
    Fun<Elem::dim, 3> const & exact,
    long double const expectedError,
    long double const tolerance)
{
  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using QR = typename LagrangeFE<Elem_T, 1>::RecommendedQR;
  using FESpaceRhs_T =
      FESpace<Mesh_T, typename LagrangeFE<Elem_T, 1>::RefFE_T, QR, Elem_T::dim>;
  using FESpaceRT0_T =
      FESpace<Mesh_T, typename RaviartThomasFE<Elem_T, 0>::RefFE_T, QR>;
  using FESpaceP0_T = FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, QR>;

  // flux feSpace
  using Facet_T = typename Elem_T::Facet_T;
  using MeshFacet_T = Mesh<Facet_T>;
  using FESpaceFacet_T = FESpace<
      MeshFacet_T,
      typename LagrangeFE<Facet_T, 0>::RefFE_T,
      typename LagrangeFE<Facet_T, 0>::RecommendedQR,
      2>;

  std::cout << Utils::separator << "test " << test_id << std::endl;
  test_id++;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T{meshCoarse}};
  for (uint r = 0U; r < refinements; ++r)
  {
    std::unique_ptr<Mesh_T> newMesh{new Mesh_T{}};
    uniformRefine2d(*mesh, *newMesh);
    mesh.reset(newMesh.release());
  }

  std::unique_ptr<MeshFacet_T> facetMesh{new MeshFacet_T};
  buildFacetMesh(*facetMesh, *mesh);

  FESpaceRhs_T feSpaceRhs{*mesh};
  Var uExact{"uExact"};
  // interpolateAnalyticFunction(inputFun, feSpace1, u1.data);
  projectAnalyticFunction(exact, feSpaceRhs, uExact.data);
  // std::cout << "uExact:\n" << uExact.data << std::endl;

  FESpaceFacet_T feSpaceFacet{*facetMesh};
  Vec fluxes;
  interpolateOnFacets<Component::NORMAL>(fluxes, feSpaceFacet, uExact.data, feSpaceRhs);
  // std::cout << "fluxesExact:\n" << fluxes << std::endl;

  FESpaceRT0_T feSpaceRT0{*mesh};
  Var uRT0{"uRT0"};
  l2Projection(uRT0.data, feSpaceRT0, uExact.data, feSpaceRhs);
  // std::cout << "uRT0:\n" << uRT0.data << std::endl;

  Var uRT0proj{"uRT0proj"};
  auto const exact3d = [&exact](Vec3 const & p) { return promote<3>(exact(p)); };
  projectAnalyticFunction(exact3d, feSpaceRT0, uRT0proj.data);
  // std::cout << "uRT0proj:\n" << uRT0proj.data << std::endl;

  Var uRT0div{"uRT0div"};
  uint const sizeRT0 = feSpaceRT0.dof.size;
  FESpaceP0_T feSpaceLambda{*mesh, sizeRT0};
  uint const sizeP0 = feSpaceLambda.dof.size;
  auto bcs = std::tuple{
      BCEss{feSpaceRT0, side::BOTTOM},
      BCEss{feSpaceRT0, side::RIGHT},
      BCEss{feSpaceRT0, side::LEFT},
      BCEss{feSpaceRT0, side::TOP},
  };
  std::get<0>(bcs) << exact3d;
  std::get<1>(bcs) << exact3d;
  std::get<2>(bcs) << exact3d;
  std::get<3>(bcs) << exact3d;

  Builder builder{sizeRT0 + sizeP0};
  builder.buildLhs(std::tuple{AssemblyVectorMass{1.0, feSpaceRT0}}, bcs);
  builder.buildCoupling(
      AssemblyVectorGrad(1.0, feSpaceRT0, feSpaceLambda), bcs, std::tuple{});
  builder.buildCoupling(
      AssemblyVectorDiv(1.0, feSpaceLambda, feSpaceRT0), std::tuple{}, bcs);
  // builder.buildLhs(std::tuple{AssemblyMass{1.0, feSpaceLambda}}, std::tuple{});
  builder.buildRhs(
      std::tuple{AssemblyProjection(1.0, uExact.data, feSpaceRhs, feSpaceRT0)}, bcs);
  builder.closeMatrix();
  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b.transpose() << std::endl;

  LUSolver solver;
  solver.compute(builder.A);
  Vec const sol = solver.solve(builder.b);
  // std::cout << "uRT0div (+lambda):\n" << sol.transpose() << std::endl;
  uRT0div.data = sol.head(sizeRT0);
  Var lambda("lambda", sizeP0);
  lambda.data = sol.tail(sizeP0);

  Var exactRT0{"exactRT0"};
  interpolateAnalyticFunction(exact3d, feSpaceRT0, exactRT0.data);
  // std::cout << "exactRT0:\n" << exactRT0.data << std::endl;

  IOManager ioRhs{feSpaceRhs, "output_projrt0err/uExact"};
  ioRhs.print({uExact});
  IOManagerP0 ioRT0{feSpaceRT0, "output_projrt0err/u"};
  ioRT0.print(std::tuple{uRT0, uRT0div, exactRT0});
  IOManager ioLambda{feSpaceLambda, "output_projrt0err/lambda"};
  ioLambda.print({lambda});

  auto const error = computeErrorL2(exact3d, uRT0div.data, feSpaceRT0);
  std::cout << "l2 error: " << error << std::endl;

  // auto const error = (uRT0div.data - exactRT0.data).norm();

  // return checkError({error}, {expectedError}, tolerance);
  return 0;
}

int main()
{
  std::bitset<16> tests;

  auto const eps = std::numeric_limits<double>::epsilon();

  // Triangle

  // std::array<Vec3, 3> const refTrianglePts = {
  //     Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}};

  // tests[0] = test<Triangle>(
  //     refTrianglePts,
  //     [&constFun](Vec3 const & p) { return narrow<2>(constFun(p)); },
  //     eps,
  //     eps);

  // tests[1] = test<Triangle>(
  //     refTrianglePts,
  //     [&linearFun](Vec3 const & p) { return narrow<2>(linearFun(p)); },
  //     2 * eps,
  //     eps);

  // std::array<Vec3, 3> const trianglePts = {
  //     Vec3{1.0, 0.0, 0.0}, Vec3{3.0, 1.0, 0.0}, Vec3{0.0, 5.0, 0.0}};

  // tests[2] = test<Triangle>(
  //     trianglePts,
  //     [&constFun](Vec3 const & p) { return narrow<2>(constFun(p)); },
  //     4 * eps,
  //     eps);

  // tests[3] = test<Triangle>(
  //     trianglePts,
  //     [&linearFun](Vec3 const & p) { return narrow<2>(linearFun(p)); },
  //     8 * eps,
  //     eps);

  // Quad

  std::unique_ptr<Mesh<Quad>> mesh{new Mesh<Quad>};

  ParameterDict meshConfig;
  meshConfig["origin"] = Vec3{-1.0, -1.0, 0.0};
  meshConfig["length"] = Vec3{2.0, 2.0, 0.0};
  meshConfig["n"] = std::array<uint, 3>{4U, 4U, 0U};
  meshConfig["flags"] =
      MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS;

  buildHyperCube(*mesh, meshConfig);

  for (uint r = 0; r < 6U; ++r)
  {
    tests[4 + r] = test<Quad>(
        *mesh,
        r,
        [](Vec3 const & p)
        {
          return Vec2{
              -(1. - p[0] * p[0]) * (1. - p[0] * p[0]) * (1. - p[1] * p[1]) * p[1],
              +(1. - p[0] * p[0]) * p[0] * (1. - p[1] * p[1]) * (1. - p[1] * p[1])};
        },
        2 * eps,
        eps);
  }

  // tests[5] = test<Quad>(
  //     refQuadPts,
  //     [&linearFun](Vec3 const & p) { return narrow<2>(linearFun(p)); },
  //     3 * eps,
  //     eps);

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
  // std::array<Vec3, 4> const quadPts = {
  //     Vec3{-1.0, -1.0, 0.0},
  //     Vec3{1.0, -1.0, 0.0},
  //     Vec3{1.0, 0.5, 0.0},
  //     Vec3{-1.0, 1.0, 0.0},
  // };

  // tests[6] = test<Quad>(
  //     quadPts,
  //     [&constFun](Vec3 const & p) { return narrow<2>(constFun(p)); },
  //     7.021666937153402e-16,
  //     eps);

  // tests[7] = test<Quad>(
  //     quadPts,
  //     [&linearFun](Vec3 const & p) { return narrow<2>(linearFun(p)); },
  //     9.085915191948722e-02,
  //     eps);

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
  //  //       axis[0] * axis[1] * (1. - std::cos(theta)) - axis[2] *
  //  std::sin(theta),
  //  //       axis[0] * axis[2] * (1. - std::cos(theta)) + axis[1] *
  //  std::sin(theta),
  //  //       axis[1] * axis[0] * (1. - std::cos(theta)) + axis[2] *
  //  std::sin(theta),
  //  //       std::cos(theta) + axis[1] * axis[1] * (1. - std::cos(theta)),
  //  //       axis[1] * axis[2] * (1. - std::cos(theta)) - axis[0] *
  //  std::sin(theta),
  //  //       axis[2] * axis[0] * (1. - std::cos(theta)) - axis[1] *
  //  std::sin(theta),
  //  //       axis[2] * axis[1] * (1. - std::cos(theta)) + axis[0] *
  //  std::sin(theta),
  //  //       std::cos(theta) + axis[2] * axis[2] * (1. -
  //  std::cos(theta))).finished();

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
