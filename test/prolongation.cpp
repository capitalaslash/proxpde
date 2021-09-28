#include "def.hpp"

#include "fespace.hpp"
#include "feutils.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"
#include "reffe.hpp"
#include "var.hpp"

template <typename RefFE>
int test()
{
  using Elem_T = typename RefFE::GeoElem_T;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T =
      FESpace<Mesh_T, RefFE, typename LagrangeFE<Elem_T, 0>::RecommendedQR>;

  std::unique_ptr<Mesh_T> meshCoarse{new Mesh_T};
  // referenceMesh(*meshCoarse);
  buildHyperCube(*meshCoarse, {0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {2, 3, 0});

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

    for (short_T iFine = 0; iFine < RefFE_T::numDOFs; ++iFine)
    {
      auto const dofFine = feSpaceFine.dof.getId(eFine.id, iFine);
      for (short_T iCoarse = 0; iCoarse < RefFE_T::numDOFs; ++iCoarse)
      {
        auto const dofCoarse = feSpaceCoarse.dof.getId(eCoarse.id, iCoarse);
        if (!done.contains(std::pair{dofFine, dofCoarse}))
        {
          double const value = RefFE_T::embeddingMatrix[posCoarse](iFine, iCoarse);
          triplets.emplace_back(dofFine, dofCoarse, value);
          done.insert(std::pair{dofFine, dofCoarse});
        }
      }
    }
  }

  Mat<StorageType::RowMajor> prol(
      meshFine->pointList.size(), meshCoarse->pointList.size());
  prol.setFromTriplets(triplets.begin(), triplets.end());
  // std::cout << prol << std::endl;

  FEVar uCoarse{"uCoarse", feSpaceCoarse};
  interpolateAnalyticFunction(
      [](Vec3 const & p) { return pow(1. - p(0), 2); }, feSpaceCoarse, uCoarse.data);
  // std::cout << "uCoarse: " << uCoarse.data.transpose() << std::endl;

  FEVar uFine{"uFine", feSpaceFine};
  uFine.data = prol * uCoarse.data;
  // std::cout << "uFine: " << uFine.data.transpose() << std::endl;

  IOManager ioCoarse{
      feSpaceCoarse, std::string{"output_prol/coarse_"} + XDMFTraits<RefFE>::shapeName};
  ioCoarse.print(std::tuple{uCoarse});
  IOManager ioFine{
      feSpaceFine, std::string{"output_prol/fine_"} + XDMFTraits<RefFE>::shapeName};
  ioFine.print(std::tuple{uFine});

  return 0;
}

int main()
{
  std::bitset<2> tests;

  tests[0] = test<RefTriangleP1>();
  tests[1] = test<RefQuadQ1>();

  return tests.any();
}
