#include "def.hpp"

#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

using namespace proxpde;

template <typename Elem>
void test(uint expectedElems, uint expectedPts)
{
  using Elem_T = Elem;
  // using Facet_T = typename Elem_T::Facet_T;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 0>::RefFE_T,
      typename LagrangeFE<Elem_T, 0>::RecommendedQR>;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  buildHyperCube(
      *mesh,
      {0.0, 0.0, 0.0},
      {1.0, 1.0, 1.0},
      {4U, 4U, 4U},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS);

  FESpace_T feSpace{*mesh};
  FEVar id{"id", feSpace};
  FEVar volume{"volume", feSpace};
  std::ranges::for_each(
      mesh->elementList,
      [&id, &volume](auto const & elem)
      {
        id.data[elem.id] = elem.id;
        volume.data[elem.id] = elem.volume();
      });

  IOManager io{feSpace, "output_block/mesh"};
  io.print({id, volume});

  std::unique_ptr<Mesh_T> block{new Mesh_T};
  extractBlock(
      *mesh,
      *block,
      [](Elem_T const & elem)
      {
        auto const mp = elem.midpoint();
        if (mp[0] < 0.4 && mp[1] < 0.4)
          return true;
        return false;
      });

  FESpace_T feSpaceBlock{*block};
  FEVar idBlock{"id", feSpaceBlock};
  FEVar volumeBlock{"volume", feSpaceBlock};
  std::ranges::for_each(
      block->elementList,
      [&idBlock, &volumeBlock](auto const & elem)
      {
        idBlock.data[elem.id] = elem.id;
        volumeBlock.data[elem.id] = elem.volume();
      });

  IOManager ioBlock{feSpaceBlock, "output_block/block"};
  ioBlock.print({idBlock, volumeBlock});

  assert(block->elementList.size() == expectedElems);
  assert(block->pointList.size() == expectedPts);
}

int main()
{
  test<Triangle>(4U, 6U);
  test<Quad>(4U, 9U);
  test<Tetrahedron>(52U, 33U);
  test<Hexahedron>(16U, 45U);

  return 0;
}
