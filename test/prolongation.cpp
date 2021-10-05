#include "def.hpp"

#include "fespace.hpp"
#include "feutils.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"
#include "reffe.hpp"
#include "var.hpp"

template <typename RefFE, typename Function>
int test(Function const & f)
{
  using Elem_T = typename RefFE::GeoElem_T;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T =
      FESpace<Mesh_T, RefFE, typename LagrangeFE<Elem_T, 1>::RecommendedQR>;

  std::unique_ptr<Mesh_T> meshCoarse{new Mesh_T};
  // referenceMesh(*meshCoarse);
  buildHyperCube(
      *meshCoarse,
      {0.0, 0.0, 0.0},
      {1.0, 1.0, 0.0},
      {1, 2, 0},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);
  // readGMSH(
  //     *meshCoarse,
  //     "square_uns.msh",
  //     MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);

  std::unique_ptr<Mesh_T> meshFine{new Mesh_T};
  uniformRefine2d(*meshCoarse, *meshFine);

  FESpace_T feSpaceCoarse{*meshCoarse};
  FESpace_T feSpaceFine{*meshFine};

  std::vector<Triplet> triplets;
  std::set<std::pair<DOFid_T, DOFid_T>> done;
  using RefFE_T = typename FESpace_T::RefFE_T;
  for (auto const & eFine: meshFine->elementList)
  {
    auto const & eCoarse = *eFine.parent.ptr;
    auto const posCoarse = eFine.parent.corner;

    for (short_T iCoarse = 0; iCoarse < RefFE_T::numDOFs; ++iCoarse)
    {
      auto const dofCoarse = feSpaceCoarse.dof.getId(eCoarse.id, iCoarse);
      double sign = 1.0;
      if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
      {
        sign =
            (eCoarse.facets[iCoarse]->facingElem[0].ptr->id != eCoarse.id) ? -1.0 : 1.0;
      }
      for (short_T iFine = 0; iFine < RefFE_T::numDOFs; ++iFine)
      {
        auto const dofFine = feSpaceFine.dof.getId(eFine.id, iFine);
        if (!done.contains(std::pair{dofFine, dofCoarse}))
        {
          double const value = RefFE_T::embeddingMatrix[posCoarse](iFine, iCoarse);
          triplets.emplace_back(dofFine, dofCoarse, sign * value);
          done.insert(std::pair{dofFine, dofCoarse});
        }
      }
    }
  }

  Mat<StorageType::RowMajor> prol(feSpaceFine.dof.size, feSpaceCoarse.dof.size);
  prol.setFromTriplets(triplets.begin(), triplets.end());
  // std::cout << prol << std::endl;

  FEVar uCoarse{"uCoarse", feSpaceCoarse};
  interpolateAnalyticFunction(f, feSpaceCoarse, uCoarse.data);
  // std::cout << "uCoarse: " << uCoarse.data.transpose() << std::endl;

  FEVar uFine{"uFine", feSpaceFine};
  uFine.data = prol * uCoarse.data;
  // std::cout << "uFine: " << uFine.data.transpose() << std::endl;

  if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
  {
    IOManager ioCoarse{
        feSpaceCoarse,
        std::string{"output_prol/coarse_"} + RefFEtoString<RefFE_T>::name};
    ioCoarse.print(std::tuple{uCoarse});
    IOManager ioFine{
        feSpaceFine, std::string{"output_prol/fine_"} + RefFEtoString<RefFE_T>::name};
    ioFine.print(std::tuple{uFine});
  }
  else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
  {
    // project on P0
    using FESpaceP0_T = FESpace<
        Mesh_T,
        typename LagrangeFE<Elem_T, 0>::RefFE_T,
        typename LagrangeFE<Elem_T, 1>::RecommendedQR,
        Elem_T::dim>;

    FESpaceP0_T feSpaceP0Coarse{*meshCoarse};
    Var uP0Coarse{"uP0Coarse"};
    l2Projection(uP0Coarse.data, feSpaceP0Coarse, uCoarse.data, feSpaceCoarse);

    IOManager ioCoarse{
        feSpaceP0Coarse,
        std::string{"output_prol/coarse_"} + RefFEtoString<RefFE_T>::name};
    ioCoarse.print(std::tuple{uP0Coarse});

    FESpaceP0_T feSpaceP0Fine{*meshFine};
    Var uP0Fine{"uP0Fine"};
    l2Projection(uP0Fine.data, feSpaceP0Fine, uFine.data, feSpaceFine);
    // std::cout << "uP0Fine: " << uP0Fine.data.transpose() << std::endl;

    IOManager ioFine{
        feSpaceP0Fine, std::string{"output_prol/fine_"} + RefFEtoString<RefFE_T>::name};
    ioFine.print(std::tuple{uP0Fine});

    // print facet mesh
    using Facet_T = typename Elem_T::Facet_T;
    using MeshFacet_T = Mesh<Facet_T>;
    using FESpaceFacet_T = FESpace<
        MeshFacet_T,
        typename LagrangeFE<Facet_T, 0>::RefFE_T,
        typename LagrangeFE<Facet_T, 0>::RecommendedQR>;

    std::unique_ptr<MeshFacet_T> facetMeshCoarse{new MeshFacet_T};
    buildFacetMesh(*facetMeshCoarse, *meshCoarse);
    FESpaceFacet_T feSpaceFacetCoarse{*facetMeshCoarse};
    IOManager ioFluxCoarse{
        feSpaceFacetCoarse,
        std::string{"output_prol/flux_coarse_"} + RefFEtoString<RefFE_T>::name};
    Var fluxCoarse{"fluxCoarse"};
    fluxCoarse.data = Vec::Zero(static_cast<uint>(facetMeshCoarse->elementList.size()));
    for (auto const & facet: meshCoarse->facetList)
    {
      auto const [insideElemPtr, side] = facet.facingElem[0];
      auto const dofId = feSpaceCoarse.dof.getId(insideElemPtr->id, side);
      fluxCoarse.data[facet.id] = uCoarse[dofId];
    }
    ioFluxCoarse.print({fluxCoarse});

    std::unique_ptr<MeshFacet_T> facetMeshFine{new MeshFacet_T};
    buildFacetMesh(*facetMeshFine, *meshFine);
    FESpaceFacet_T feSpaceFacetFine{*facetMeshFine};
    IOManager ioFluxFine{
        feSpaceFacetFine,
        std::string{"output_prol/flux_fine_"} + RefFEtoString<RefFE_T>::name};
    Var fluxFine{"fluxFine"};
    fluxFine.data = Vec::Zero(static_cast<uint>(facetMeshFine->elementList.size()));
    for (auto const & facet: meshFine->facetList)
    {
      auto const [insideElemPtr, side] = facet.facingElem[0];
      auto const dofId = feSpaceFine.dof.getId(insideElemPtr->id, side);
      fluxFine.data[facet.id] = uFine[dofId];
    }
    ioFluxFine.print({fluxFine});
  }

  return 0;
}

int main()
{
  std::bitset<4> tests;

  auto const fScalar = [](Vec3 const & p) { return pow(1. - p(0), 2); };
  auto const fVector = [](Vec3 const & p) { return Vec3{p[0], -p[1], 0.0}; };

  tests[0] = test<RefTriangleP1>(fScalar);
  tests[1] = test<RefTriangleRT0>(fVector);
  tests[2] = test<RefQuadQ1>(fScalar);
  tests[3] = test<RefQuadRT0>(fVector);

  return tests.any();
}
