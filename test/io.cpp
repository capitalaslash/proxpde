#include <array>      // for array
#include <cmath>      // for sqrt
#include <cstdlib>    // for uint
#include <filesystem> // for path
#include <memory>     // for unique_ptr
#include <optional>   // for optional
#include <string>     // for stoul, basic_string
#include <vector>     // for vector

#include <proxpde/def.hpp>       // for Vec3
#include <proxpde/fe.hpp>        // for LagrangeFE
#include <proxpde/fespace.hpp>   // for FESpace, interpolateAnalyticFunction
#include <proxpde/geo.hpp>       // for Hexahedron
#include <proxpde/iomanager.hpp> // for IOManager
#include <proxpde/mesh.hpp>      // for Mesh, buildHyperCube
#include <proxpde/qr.hpp>        // for GaussQR
#include <proxpde/var.hpp>       // for Var

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
    numElems[0] = static_cast<uint>(std::stoul(argv[1]));
    numElems[1] = static_cast<uint>(std::stoul(argv[2]));
    numElems[2] = static_cast<uint>(std::stoul(argv[3]));
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
