#include "def.hpp"

#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"

template <typename Elem>
void test()
{
  using Elem_T = Elem;
  using Facet_T = typename Elem_T::Facet_T;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 0>::RefFE_T,
      typename LagrangeFE<Elem_T, 0>::RecommendedQR>;

  std::unique_ptr<Mesh_T> meshCoarse{new Mesh_T};

  // referenceMesh(*mesh);
  buildHyperCube(
      *meshCoarse,
      {0.0, 0.0, 0.0},
      {1.0, 1.0, 1.0},
      {2U, 3U, 4U},
      MeshFlags::INTERNAL_FACETS);
  // readGMSH(*mesh, "square_uns.msh");
  // readGMSH(*mesh, "square_q.msh");

  std::cout << Utils::separator << "mesh coarse\n" << *meshCoarse << std::endl;

  std::unique_ptr<Mesh_T> meshFine{new Mesh_T};

  uniformRefine(*meshCoarse, *meshFine);

  std::cout << Utils::separator << "mesh fine:\n" << *meshFine << std::endl;

  uint checked = 0;
  for (auto const & f: meshFine->facetList)
  {
    // check that our facing elem is the child of the facing elem of our parent
    // new internal facets cannot have a parent, so exclude them
    if (f.parent)
    {
      [[maybe_unused]] auto const inParent = f.facingElem[0].ptr->parent.ptr;
      [[maybe_unused]] auto const parentIn = f.parent.ptr->facingElem[0].ptr;
      assert(inParent->id == parentIn->id);
      checked++;
      if (f.facingElem[1])
      {
        [[maybe_unused]] auto const outParent = f.facingElem[1].ptr->parent.ptr;
        [[maybe_unused]] auto const parentOut = f.parent.ptr->facingElem[1].ptr;
        assert(outParent->id == parentOut->id);
        checked++;
      }
    }
  }
  [[maybe_unused]] uint const bdSize = std::count_if(
      meshCoarse->facetList.begin(),
      meshCoarse->facetList.end(),
      [](Facet_T const & f) { return f.onBoundary(); });
  // we do 1 check for each boundary facet on the fine mesh + 2 checks for each internal
  // facets on the fine mesh
  assert(
      checked ==
      Facet_T::numChildren * (bdSize + 2 * (meshCoarse->facetList.size() - bdSize)));
  std::cout << "facet checks performed: " << checked << std::endl;

  FESpace_T feSpace{*meshFine};
  FEVar id{"id", feSpace};
  for (auto const & elem: meshFine->elementList)
  {
    id.data[elem.id] = elem.id;
  }
  IOManager io{feSpace, "output_refine/fine"};
  io.print(std::tuple{id});
}

int main()
{
  test<Triangle>();

  test<Quad>();

  test<Tetrahedron>();

  test<Hexahedron>();

  return 0;
}
