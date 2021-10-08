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
int test(Function const & f, double const expectedNorm)
{
  using Elem_T = typename RefFE::GeoElem_T;
  using Facet_T = typename Elem_T::Facet_T;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T =
      FESpace<Mesh_T, RefFE, typename LagrangeFE<Elem_T, 1>::RecommendedQR>;

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
  uniformRefine2d(*meshCoarse, *meshFine);

  FESpace_T feSpaceCoarse{*meshCoarse};
  FESpace_T feSpaceFine{*meshFine};

  std::vector<Triplet> triplets;
  std::set<std::pair<DOFid_T, DOFid_T>> done;
  using RefFE_T = typename FESpace_T::RefFE_T;
  for (auto const & eCoarse: meshCoarse->elementList)
  {
    for (short_T iChild = 0; iChild < RefFE_T::numChildren; ++iChild)
    {
      auto const & eFine = *eCoarse.children[iChild].ptr;
      assert(eFine.parent.corner == iChild);
      for (short_T iCoarse = 0; iCoarse < RefFE_T::numDOFs; ++iCoarse)
      {
        auto const dofCoarse = feSpaceCoarse.dof.getId(eCoarse.id, iCoarse);
        // don't need this on restriction?!
        double sign = 1.0;
        // if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
        // {
        //   sign = (eCoarse.facets[iChild]->facingElem[0].ptr->id != eCoarse.id) ? -1.0
        //                                                                        : 1.0;
        // }
        for (short_T iFine = 0; iFine < RefFE_T::numDOFs; ++iFine)
        {
          auto const dofFine = feSpaceFine.dof.getId(eFine.id, iFine);
          if (!done.contains({dofCoarse, dofFine}))
          {
            double const value = RefFE_T::embeddingMatrix[iChild](iFine, iCoarse);
            if (std::fabs(value) > 1.e-12)
            {
              triplets.emplace_back(dofCoarse, dofFine, sign * value);
            }
            // RT elements require to sum contribution from both side of the faces
            if constexpr (family_v<RefFE_T> != FamilyType::RAVIART_THOMAS)
            {
              done.insert({dofCoarse, dofFine});
            }
          }
        }
      }
    }
  }

  Mat<StorageType::RowMajor> rest(feSpaceCoarse.dof.size, feSpaceFine.dof.size);
  rest.setFromTriplets(triplets.begin(), triplets.end());
  // std::cout << "rest:\n" << rest << std::endl;

  // rescale each row so that it's weights sum to 1.
  // this ensures that constant solutions are restricted to the same value.
  for (int row = 0; row < rest.rows(); ++row)
  {
    auto const start = rest.outerIndexPtr()[row];
    auto const end = rest.outerIndexPtr()[row + 1];
    // std::cout << "row " << row << " (" << end - start << "): ";
    auto sum = 0.0;
    for (int clm = start; clm < end; ++clm)
    {
      sum += rest.valuePtr()[clm];
      // std::cout << rest.valuePtr()[clm] << " (" << rest.innerIndexPtr()[clm] << ") ";
    }
    // std::cout << std::endl;

    if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
    {
      sum /= Facet_T::numChildren;
    }

    // no line can be with sum 0
    assert(std::fabs(sum) > 1.e-12);
    for (int clm = start; clm < end; ++clm)
    {
      rest.valuePtr()[clm] /= sum;
    }
  }

  FEVar uFine{"u", feSpaceFine};
  interpolateAnalyticFunction(f, feSpaceFine, uFine.data);
  std::cout << "uFine: " << uFine.data.transpose() << std::endl;

  FEVar uCoarse{"u", feSpaceCoarse};
  uCoarse.data = rest * uFine.data;
  std::cout << "uCoarse: " << uCoarse.data.transpose() << std::endl;

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

  tests[0] = test<RefTriangleP1>(fScalar, 1.762571852935e+00);
  tests[1] = test<RefTriangleRT0>(fVector, 1.433720877840e+00);
  tests[2] = test<RefQuadQ1>(fScalar, 1.759929209669e+00);
  tests[3] = test<RefQuadRT0>(fVector, 1.092906420717e+00);

  return tests.any();
}
