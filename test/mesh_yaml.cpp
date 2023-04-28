#include "def.hpp"

#include "geo.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;

  MilliTimer t;

  ParameterDict config;

  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    // config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    // config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
    // config["mesh"]["n"] = std::array{2U, 3U, 0U};
    // config["mesh"]["flags"] =
    //     MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS;
    config = YAML::Load("mesh:\n origin: [0.0, 0.0, 0.0]\n length: [1.0, 1.0, 1.0]\n "
                        "n: [2, 3, 0]\n flags: ['INTERNAL_FACETS']");
  }

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.print();

  if (!(mesh->flags & MeshFlags::INTERNAL_FACETS).any())
  {
    std::cerr << "missing flag INTERNAL_FACETS" << std::endl;
    return 1;
  }

  return 0;
}
