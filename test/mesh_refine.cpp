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
      {1.0, 1.0, 0.0},
      {2, 3, 0},
      MeshFlags::INTERNAL_FACETS);
  // readGMSH(*mesh, "square_uns.msh");
  // readGMSH(*mesh, "square_q.msh");

  std::cout << Utils::separator << "mesh coarse\n" << *meshCoarse << std::endl;

  std::unique_ptr<Mesh_T> meshFine{new Mesh_T};

  uniformRefine2d(*meshCoarse, *meshFine);

  std::cout << Utils::separator << "mesh fine:\n" << *meshFine << std::endl;

  for (auto const & f: meshFine->facetList)
  {
    // check that our facing elem is the child of the facing elem of our parent
    assert(f.facingElem[0].ptr->parent.ptr->id == f.parent.ptr->facingElem[0].ptr->id);
    if (f.facingElem[1])
    {
      assert(
          f.facingElem[1].ptr->parent.ptr->id == f.parent.ptr->facingElem[1].ptr->id);
    }
  }

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

  return 0;
}
