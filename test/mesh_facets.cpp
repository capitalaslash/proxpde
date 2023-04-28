#include "def.hpp"

#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"
#include "xdmf_traits.hpp"

using namespace proxpde;

template <typename Elem>
void test()
{
  marker_T const testMarker = 10;

  using Elem_T = Elem;
  using Mesh_T = Mesh<Elem_T>;
  using RefFE_T = typename LagrangeFE<Elem_T, 0>::RefFE_T;
  using QR_T = typename LagrangeFE<Elem_T, 0>::RecommendedQR;
  using FESpace_T = FESpace<Mesh_T, RefFE_T, QR_T>;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // referenceMesh(*mesh);
  buildHyperCube(
      *mesh, {0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {2, 1, 0}, MeshFlags::INTERNAL_FACETS);
  // readGMSH(*mesh, "square_uns.msh");
  // readGMSH(*mesh, "square_q.msh");
  for (auto & facet: mesh->facetList)
  {
    if (!facet.onBoundary())
    {
      facet.marker = testMarker;
    }
  }
  typename Mesh_T::FacetList_T facetsOrig = mesh->facetList;

  std::cout << Utils::separator << "mesh before\n" << *mesh << std::endl;

  buildFacets(*mesh, MeshFlags::INTERNAL_FACETS);

  std::cout << Utils::separator << "mesh after\n" << *mesh << std::endl;

  // check that buildFacets() has preserved an internal facet
  for ([[maybe_unused]] auto const & facet: mesh->facetList)
  {
    assert(facet == facetsOrig[facet.id]);
  }

  FESpace_T feSpace{*mesh};
  FEVar id{"id", feSpace};
  for (auto const & elem: mesh->elementList)
  {
    id.data[elem.id] = elem.id;
  }
  IOManager io{
      feSpace, std::string{"output_meshfacets/mesh_"} + XDMFTraits<RefFE_T>::shapeName};
  io.print(std::tuple{id});
}

int main()
{
  test<Triangle>();
  test<Quad>();

  return 0;
}
