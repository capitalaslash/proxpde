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

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  // referenceMesh(*mesh);
  buildHyperCube(
      *mesh, {0., 0., 0.}, {1., 1., 0.}, {2, 1, 0}, MeshFlags::INTERNAL_FACETS);
  // readGMSH(*mesh, "square_uns.msh");
  // readGMSH(*mesh, "square_q.msh");

  std::cout << "original mesh\n" << *mesh << std::endl;

  std::unique_ptr<Mesh_T> newMesh{new Mesh_T};

  uniformRefine2d(*mesh, *newMesh);

  std::cout << Utils::separator << "refined mesh:\n" << *newMesh << std::endl;

  FESpace_T feSpace{*newMesh};
  FEVar id{"id", feSpace};
  for (auto const & elem: newMesh->elementList)
  {
    id.data[elem.id] = elem.id;
  }
  IOManager io{feSpace, "output_refine/new_mesh"};
  io.print(std::tuple{id});
}

int main()
{
  test<Triangle>();

  test<Quad>();

  return 0;
}
