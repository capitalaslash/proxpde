#include "def.hpp"

#include "fespace.hpp"
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

  Prolongator prol{feSpaceCoarse, feSpaceFine};

  FEVar uCoarse{"u", feSpaceCoarse};
  interpolateAnalyticFunction(f, feSpaceCoarse, uCoarse.data);
  // std::cout << "uCoarse: " << uCoarse.data.transpose() << std::endl;

  FEVar uFine{"u", feSpaceFine};
  uFine.data = prol.mat * uCoarse.data;
  // std::cout << "uFine: " << uFine.data.transpose() << std::endl;

  if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
  {
    IOManager ioCoarse{
        feSpaceCoarse,
        std::string{"output_prol/coarse_"} + RefFEtoString<RefFE_T>::name};
    ioCoarse.print({uCoarse});

    IOManager ioFine{
        feSpaceFine, std::string{"output_prol/fine_"} + RefFEtoString<RefFE_T>::name};
    ioFine.print({uFine});
  }
  else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
  {
    IOManagerP0 ioCoarse{
        feSpaceCoarse,
        std::string{"output_prol/coarse_"} + RefFEtoString<RefFE_T>::name};
    ioCoarse.print(std::vector{uCoarse});

    IOManagerP0 ioFine{
        feSpaceFine, std::string{"output_prol/fine_"} + RefFEtoString<RefFE_T>::name};
    ioFine.print(std::vector{uFine});

    // print facet mesh
    IOManagerFacet ioFacetCoarse{
        feSpaceCoarse,
        std::string{"output_prol/flux_coarse_"} + RefFEtoString<RefFE_T>::name};
    ioFacetCoarse.print(std::vector{uCoarse});

    IOManagerFacet ioFacetFine{
        feSpaceFine,
        std::string{"output_prol/flux_fine_"} + RefFEtoString<RefFE_T>::name};
    ioFacetFine.print(std::vector{uFine});
  }

  return checkError({uFine.data.norm()}, {expectedNorm});
}

int main()
{
  std::bitset<4> tests;

  auto const fScalar = [](Vec3 const & p) { return pow(1. - p(0), 2); };
  auto const fVector = [](Vec3 const & p) { return Vec3{p[0], -p[1], 0.0}; };

  tests[0] = test<RefTriangleP1>(fScalar, 3.009757793463e+00);
  tests[1] = test<RefTriangleRT0>(fVector, 1.263812574009e+00);
  tests[2] = test<RefQuadQ1>(fScalar, 3.009757793463e+00);
  tests[3] = test<RefQuadRT0>(fVector, 9.718253158076e-01);

  fmt::print("tests: {}\n", tests.to_string());

  return tests.any();
}
