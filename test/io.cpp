#include "def.hpp"

#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

using Elem_T = Hexahedron;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T, RefHexahedronQ2, GaussQR<Hexahedron, 1>>;

int main(int argc, char * argv[])
{
  std::array<uint, 3> numElems = {{1, 1, 1}};
  if (argc == 4)
  {
    numElems[0] = static_cast<uint>(std::stoi(argv[1]));
    numElems[1] = static_cast<uint>(std::stoi(argv[2]));
    numElems[2] = static_cast<uint>(std::stoi(argv[3]));
  }

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{1., 1., 1.}, numElems);

  FESpace_T feSpace{*mesh};

  Var u{"u"};

  interpolateAnalyticFunction(
      [](Vec3 const & p) { return p(0) * p(0); }, feSpace, u.data);

  IOManager io{feSpace, "output_io/u"};
  io.print({u});

  return 0;
}
