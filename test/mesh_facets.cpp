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
      *mesh, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {2, 1, 3}, MeshFlags::INTERNAL_FACETS);
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

  using FESpaceV_T = FESpace<Mesh_T, RefFE_T, QR_T, 2U>;
  FESpaceV_T feSpaceV{*mesh};
  FEVar v{"v", feSpaceV};
  for (auto const & elem: mesh->elementList)
  {
    for (uint k = 0; k < FESpaceV_T::dim; k++)
    {
      if constexpr (feSpaceV.ordering == DofOrdering::INTERLEAVED)
      {
        v.data[elem.id * FESpaceV_T::dim + k] = 1.0 + k;
      }
      else
      {
        v.data[elem.id + k * feSpaceV.dof.size] = 1.0 + k;
      }
    }
  }
  auto const baseName = std::filesystem::path{"output_meshfacets/v_"}.concat(
      XDMFTraits<RefFE_T>::shapeName);
  IOManager ioV{feSpaceV, baseName};
  ioV.print(std::tuple{v});

  XDMFDoc doc{
      fmt::format("output_meshfacets/test_{}.xmf", XDMFTraits<RefFE_T>::shapeName)};
  doc.setTime(0.0);
  auto const meshFile = std::filesystem::path{baseName.filename()}.concat(".mesh.h5");
  doc.setTopology<typename FESpaceV_T::RefFE_T>(
      meshFile, feSpaceV.mesh->elementList.size());
  doc.setGeometry(meshFile, feSpaceV.dof.mapSize);

  auto const dataFile = std::filesystem::path{
      fmt::format("output_meshfacets/test_{}.h5", XDMFTraits<RefFE_T>::shapeName)};
  HDF5 h5{dataFile, HDF5FileMode::OVERWRITE};
  doc.setVar(
      dataFile.filename(),
      {v.name,
       XDMFType::VECTOR,
       XDMFCenter::CELL,
       XDMFNumberType::FLOAT,
       feSpaceV.dof.size,
       feSpaceV.dim});
  Table<double, FESpaceV_T::dim> data;
  data.resize(feSpaceV.dof.size, FESpaceV_T::dim);
  for (uint i = 0; i < feSpaceV.dof.size; i++)
  {
    for (uint k = 0; k < FESpaceV_T::dim; k++)
    {
      if constexpr (feSpaceV.ordering == DofOrdering::INTERLEAVED)
      {
        data(i, k) = v.data[i * FESpaceV_T::dim + k];
      }
      else
      {
        data(i, k) = v.data[i + k * feSpaceV.dof.size];
      }
    }
  }
  h5.print(data, v.name);

  // write and read again
  writeGMSH(*mesh, "output_meshfacets/mesh_mod.msh");

  // TODO: can read only meshes without internal facets
  // std::unique_ptr<Mesh_T> meshWritten{new Mesh_T};
  // readGMSH(*meshWritten, "output_meshfacets/mesh_mod.msh",
  // MeshFlags::INTERNAL_FACETS);

  // assert(mesh->elementList == meshWritten->elementList);
}

int main()
{
  test<Triangle>();
  test<Quad>();
  test<Tetrahedron>();
  test<Hexahedron>();

  return 0;
}
