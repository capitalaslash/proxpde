#include "def.hpp"

#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Hexahedron;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T =
      FESpace<Mesh_T, LagrangeFE<Elem_T, 2u>::RefFE_T, GaussQR<Elem_T, 1u>>;

  auto numElems = std::array{2u, 1u, 1u};
  if (argc == 4)
  {
    numElems[0] = std::stoul(argv[1]);
    numElems[1] = std::stoul(argv[2]);
    numElems[2] = std::stoul(argv[3]);
  }

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}, numElems);

  FESpace_T feSpace{*mesh};

  Var u{"u"};

  interpolateAnalyticFunction(
      [](Vec3 const & p) { return p(0) * p(0); }, feSpace, u.data);

  IOManager io{feSpace, "output_io/u"};
  io.print({u});

  return 0;
}
