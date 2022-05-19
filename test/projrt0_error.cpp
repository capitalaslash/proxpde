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
    Fun<Elem::dim, 3> const & exactU,
    scalarFun_T const & exactLambda,
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

  fmt::print("{}test {}\n", Utils::separator, test_id);
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
  projectAnalyticFunction(exactU, feSpaceRhs, uExact.data);
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
  auto const exactU3d = [&exactU](Vec3 const & p) { return promote<3>(exactU(p)); };
  projectAnalyticFunction(exactU3d, feSpaceRT0, uRT0proj.data);
  // std::cout << "uRT0proj:\n" << uRT0proj.data << std::endl;

  Var uRT0div{"uRT0div"};
  uint const sizeRT0 = feSpaceRT0.dof.size;
  FESpaceP0_T feSpaceLambda{*mesh, sizeRT0};
  FESpaceP0_T feSpaceLambda0{*mesh};
  uint const sizeP0 = feSpaceLambda.dof.size;
  auto bcsU = std::array{
      BCEss{feSpaceRT0, side::BOTTOM},
      BCEss{feSpaceRT0, side::RIGHT},
      BCEss{feSpaceRT0, side::LEFT},
      // BCEss{feSpaceRT0, side::TOP},
  };
  std::get<0>(bcsU) << exactU3d;
  std::get<1>(bcsU) << exactU3d;
  std::get<2>(bcsU) << exactU3d;
  // std::get<3>(bcs) << exactU3d;

  auto bcsLambda = std::array{
      BCEss{feSpaceLambda, side::TOP},
  };
  std::get<0>(bcsLambda) << exactLambda;

  Builder builder{sizeRT0 + sizeP0};
  builder.buildLhs(std::tuple{AssemblyVectorMass{1.0, feSpaceRT0}}, bcsU);
  builder.buildCoupling(
      AssemblyVectorGrad(1.0, feSpaceRT0, feSpaceLambda), bcsU, bcsLambda);
  builder.buildCoupling(
      AssemblyVectorDiv(1.0, feSpaceLambda, feSpaceRT0), bcsLambda, bcsU);
  // the dummy assembly is necessary to apply Dirichlet bc for lambda
  builder.buildLhs(std::tuple{AssemblyDummy{feSpaceLambda}}, bcsLambda);
  builder.buildRhs(
      std::tuple{AssemblyProjection(1.0, uExact.data, feSpaceRhs, feSpaceRT0)}, bcsU);
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
  interpolateAnalyticFunction(exactU3d, feSpaceRT0, exactRT0.data);
  // std::cout << "exactRT0:\n" << exactRT0.data << std::endl;

  IOManager ioRhs{feSpaceRhs, "output_projrt0err/uExact"};
  ioRhs.print({uExact});
  IOManagerP0 ioRT0{feSpaceRT0, "output_projrt0err/u"};
  ioRT0.print(std::tuple{uRT0, uRT0div, exactRT0});
  IOManager ioLambda{feSpaceLambda, "output_projrt0err/lambda"};
  ioLambda.print({lambda});

  auto const errorU = computeErrorL2(exactU3d, uRT0div.data, feSpaceRT0);
  fmt::print("l2 errorU: {:.12e}\n", errorU);

  auto const errorLambda = computeErrorL2(exactLambda, lambda.data, feSpaceLambda0);
  fmt::print("l2 errorLambda: {:.12e}\n", errorLambda);

  // auto const error = (uRT0div.data - exactRT0.data).norm();

  // return checkError({error}, {expectedError}, tolerance);
  return 0;
}

int main(int argc, char * argv[])
{
  std::bitset<16> tests;

  auto const eps = std::numeric_limits<double>::epsilon();

  ParameterDict config;
  config["refinements"] = 3U;

  config["triangle"]["structured"] = true;
  config["triangle"]["gmsh"] = true;
  config["quad"]["structured"] = true;
  config["quad"]["gmsh"] = true;

  config["mesh"]["origin"] = Vec3{-1.0, -1.0, 0.0};
  config["mesh"]["length"] = Vec3{2.0, 2.0, 0.0};
  config["mesh"]["n"] = std::array<uint, 3>{4U, 4U, 0U};
  config["mesh"]["flags"] =
      MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS;

  auto const funDivZero2d = std::tuple{
      [](Vec3 const & p)
      {
        return Vec2{
            -(1. - p[0] * p[0]) * (1. - p[0] * p[0]) * (1. - p[1] * p[1]) * p[1],
            +(1. - p[0] * p[0]) * p[0] * (1. - p[1] * p[1]) * (1. - p[1] * p[1])};
      },
      [](Vec3 const &) { return 0.0; }};

  if (argc > 1)
  {
    config.override(argv[1]);
  }

  // fmt::print("configuration: {}", config);

  uint const refinements = config["refinements"].as<uint>();

  if (config["triangle"]["structured"].as<bool>())
  {
    std::unique_ptr<Mesh<Triangle>> mesh{new Mesh<Triangle>};

    buildHyperCube(*mesh, config["mesh"]);

    for (uint r = 0; r < refinements; ++r)
    {
      tests[0 + r] = test(
          *mesh, r, std::get<0>(funDivZero2d), std::get<1>(funDivZero2d), 2 * eps, eps);
    }
  }

  if (config["triangle"]["gmsh"].as<bool>())
  {
    std::unique_ptr<Mesh<Triangle>> mesh{new Mesh<Triangle>};
    readGMSH(
        *mesh,
        "square_uns.msh",
        MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS);

    // rescale to [-1,1]x[-1,1]
    for (auto & pt: mesh->pointList)
    {
      pt.coord = pt.coord * 2.0 - Vec3::Constant(1.0);
    }

    for (uint r = 0; r < refinements; ++r)
    {
      tests[0 + r] = test(
          *mesh, r, std::get<0>(funDivZero2d), std::get<1>(funDivZero2d), 2 * eps, eps);
    }
  }

  if (config["quad"]["structured"].as<bool>())
  {
    std::unique_ptr<Mesh<Quad>> mesh{new Mesh<Quad>};

    buildHyperCube(*mesh, config["mesh"]);

    for (uint r = 0; r < refinements; ++r)
    {
      tests[0 + r] = test(
          *mesh, r, std::get<0>(funDivZero2d), std::get<1>(funDivZero2d), 2 * eps, eps);
    }
  }

  if (config["quad"]["gmsh"].as<bool>())
  {
    std::unique_ptr<Mesh<Quad>> mesh{new Mesh<Quad>};
    readGMSH(
        *mesh,
        "square_skew.msh",
        MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS);

    for (uint r = 0; r < refinements; ++r)
    {
      tests[0 + r] = test(
          *mesh, r, std::get<0>(funDivZero2d), std::get<1>(funDivZero2d), 2 * eps, eps);
    }
  }

  // Tetrahedron

  // Hexahedron

  std::cout << "test results: " << tests << std::endl;
  return tests.any();
}
