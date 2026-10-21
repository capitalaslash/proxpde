#include <algorithm> // for count_if
#include <cassert>   // for assert
#include <cmath>     // for fabs
#include <cstddef>   // for size_t
#include <iostream>  // for basic_ostream, operator<<, cout
#include <memory>    // for unique_ptr
#include <string>    // for char_traits, operator+, operator<<
#include <vector>    // for vector

#include <Eigen/Core> // for DenseBase::row, DenseBase<>::IsVe...
#include <fmt/base.h> // for println

#include <proxpde/def.hpp>         // for operator==, Utils
#include <proxpde/fe.hpp>          // for LagrangeFE
#include <proxpde/fespace.hpp>     // for FESpace
#include <proxpde/geo.hpp>         // for Hexahedron, Quad, Tetrahedron
#include <proxpde/iomanager.hpp>   // for IOManager
#include <proxpde/mesh.hpp>        // for operator<<, buildHyperCube, Mesh
#include <proxpde/mesh_refine.hpp> // for uniformRefine
#include <proxpde/var.hpp>         // for FEVar

using namespace proxpde;

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
  // readGMSH(*mesh, "cube_uns.msh");
  // readGMSH(*mesh, "cube_skew.msh", MeshFlags::CHECK_PLANAR_FACETS);

  std::cout << Utils::separator << "mesh coarse\n" << *meshCoarse << std::endl;

  std::unique_ptr<Mesh_T> meshFine{new Mesh_T};

  uniformRefine(*meshCoarse, *meshFine);

  std::cout << Utils::separator << "mesh fine:\n" << *meshFine << std::endl;

  auto checked = 0u;
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
  fmt::println("facet checks performed: {}", checked);
  [[maybe_unused]] auto const bdSize = static_cast<std::size_t>(std::count_if(
      meshCoarse->facetList.begin(),
      meshCoarse->facetList.end(),
      [](Facet_T const & f) { return f.onBoundary(); }));
  // we do 1 check for each boundary facet on the fine mesh + 2 checks for each internal
  // facets on the fine mesh
  // bd + 2 * (total - bd) = 2 * total - bd
  assert(checked == Facet_T::numChildren * (2 * meshCoarse->facetList.size() - bdSize));
  fmt::println("facet checks performed: {}", checked);

  FESpace_T feSpace{*meshFine};
  FEVar id{"id", feSpace};
  for (auto const & elem: meshFine->elementList)
  {
    id.data[elem.id] = elem.id;
  }
  IOManager io{feSpace, std::string{"output_refine/fine_"} + elemToStr_v<Elem>};
  io.print({id});
}

int main()
{
  test<Triangle>();

  test<Quad>();

  test<Tetrahedron>();

  test<Hexahedron>();

  return 0;
}
