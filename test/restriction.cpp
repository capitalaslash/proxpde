#include "def.hpp"

#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"
#include "multigrid.hpp"
#include "reffe.hpp"
#include "var.hpp"

using namespace proxpde;

template <typename RefFE, typename Function>
int test(Function const & f, double const expectedNorm)
{
  using Elem_T = typename RefFE::GeoElem_T;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T =
      FESpace<Mesh_T, RefFE, typename LagrangeFE<Elem_T, 1>::RecommendedQR>;
  using RefFE_T = typename FESpace_T::RefFE_T;
  using BCList_T = std::array<BCEss<FESpace_T>, 0>; // no bc considered

  std::unique_ptr<Mesh_T> meshCoarse{new Mesh_T};
  // referenceMesh(*meshCoarse);
  buildHyperCube(
      *meshCoarse,
      {0.0, 0.0, 0.0},
      {1.0, 1.0, 0.0},
      {3, 2, 0},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);
  // readGMSH(
  //     *meshCoarse,
  //     "square_uns.msh",
  //     MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);

  std::unique_ptr<Mesh_T> meshFine{new Mesh_T};
  uniformRefine(*meshCoarse, *meshFine);

  FESpace_T feSpaceCoarse{*meshCoarse};
  FESpace_T feSpaceFine{*meshFine};

  Restrictor rest{feSpaceFine, feSpaceCoarse, BCList_T{}};

  FEVar uFine{"u", feSpaceFine};
  interpolateAnalyticFunction(f, feSpaceFine, uFine.data);
  // std::cout << "uFine: " << uFine.data.transpose() << std::endl;

  FEVar uCoarse{"u", feSpaceCoarse};
  uCoarse.data = rest.mat * uFine.data;
  // std::cout << "uCoarse: " << uCoarse.data.transpose() << std::endl;

  if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
  {
    IOManager ioCoarse{
        feSpaceCoarse,
        std::string{"output_rest/coarse_"} + RefFEtoString<RefFE_T>::name};
    ioCoarse.print(std::tuple{uCoarse});

    IOManager ioFine{
        feSpaceFine, std::string{"output_rest/fine_"} + RefFEtoString<RefFE_T>::name};
    ioFine.print(std::tuple{uFine});
  }
  else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
  {
    IOManagerP0 ioCoarse{
        feSpaceCoarse,
        std::string{"output_rest/coarse_"} + RefFEtoString<RefFE_T>::name};
    ioCoarse.print(std::tuple{uCoarse});

    IOManagerP0 ioFine{
        feSpaceFine, std::string{"output_rest/fine_"} + RefFEtoString<RefFE_T>::name};
    ioFine.print(std::tuple{uFine});

    // print facet mesh
    IOManagerFacet ioFacetCoarse{
        feSpaceCoarse,
        std::string{"output_rest/flux_coarse_"} + RefFEtoString<RefFE_T>::name};
    ioFacetCoarse.print(std::tuple{uCoarse});

    IOManagerFacet ioFacetFine{
        feSpaceFine,
        std::string{"output_rest/flux_fine_"} + RefFEtoString<RefFE_T>::name};
    ioFacetFine.print(std::tuple{uFine});
  }

  return checkError({uCoarse.data.norm()}, {expectedNorm});
}

int main()
{
  std::bitset<4> tests;

  auto const fScalar = [](Vec3 const & p) { return pow(1. - p(0), 2); };
  auto const fVector = [](Vec3 const & p) { return Vec3{p[0], -p[1], 0.0}; };

  tests[0] = test<RefTriangleP1>(fScalar, 1.759875655333e+00);
  tests[1] = test<RefTriangleRT0>(fVector, 1.343709624716e+00);
  tests[2] = test<RefQuadQ1>(fScalar, 1.759929209669e+00);
  tests[3] = test<RefQuadRT0>(fVector, 1.092906420717e+00);

  fmt::print("tests: {}\n", tests.to_string());

  return tests.any();
}
