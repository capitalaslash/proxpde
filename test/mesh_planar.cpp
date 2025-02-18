#include "def.hpp"

#include "fe.hpp"
#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"
#include "var.hpp"

using namespace proxpde;

int main(int argc, char * argv[])
{
  using Elem_T = Hexahedron;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 0>::RefFE_T,
      LagrangeFE<Elem_T, 0>::RecommendedQR>;

  MilliTimer t;

  ParameterDict config;

  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    config["mesh"]["length"] = Vec3{1.0, 1.0, 1.0};
    config["mesh"]["n"] = std::array{16U, 16U, 16U};
    config["mesh"]["flags"] = Bitmask{MeshFlags::NORMALS};
  }

  t.start("generation");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.start("chevron");
  uint const n = config["mesh"]["n"].as<std::array<uint, 3>>()[0];
  double const h = 1. / n;
  for (uint k = 1; k < n; k += 2)
    for (uint j = 1; j < n; j += 2)
      for (uint i = 1; i < n; i += 2)
        mesh->pointList[i + (n + 1) * j + cepow(n + 1, 2) * k].coord[1] += h / 3;
  for (uint k = 0; k <= n; k += 2)
    for (uint j = 1; j < n; j += 2)
      for (uint i = 0; i <= n; i += 2)
        mesh->pointList[i + (n + 1) * j + cepow(n + 1, 2) * k].coord[1] -= h / 3;
  t.stop();

  t.start("check");
  uint const numNonPlanar = checkPlanarFacets(*mesh);
  fmt::print("non planar facets: {}\n", numNonPlanar);
  t.stop();

  t.start("print");
  FESpace_T feSpace{*mesh};
  FEVar vol{"normalized_volume", feSpace};
  for (auto const & elem: mesh->elementList)
  {
    vol.data[elem.id] = elem.volume() * n * n * n;
  }
  IOManager io{feSpace, "output_mesh_planar/vol"};
  io.print({vol});
  t.stop();

  t.print();

  mesh.reset(new Mesh_T());
  referenceMesh(*mesh);
  mesh->pointList[6].coord += Vec3{0.1, 0.1, 0.1};
  // with MeshFlags::CHECK_FACET_PLANAR buildNormals() correctly aborts
  buildNormals(*mesh, MeshFlags::NONE);
  uint const numNonPlanarRef = checkPlanarFacets(*mesh);
  assert(numNonPlanarRef == 3);

  return 0;
}
