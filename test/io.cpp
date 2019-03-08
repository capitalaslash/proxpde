#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

using Elem_T =Hexahedron;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T, RefHexahedronQ2, GaussQR<Hexahedron,1>>;

int main(int argc, char* argv[])
{
  array<uint,3> numPts = {{2, 2, 2}};
  if (argc == 4)
  {
    numPts[0] = static_cast<uint>(std::stoi(argv[1]));
    numPts[1] = static_cast<uint>(std::stoi(argv[2]));
    numPts[1] = static_cast<uint>(std::stoi(argv[3]));
  }

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, Vec3{0., 0., 0.}, Vec3{1., 1., 1.}, numPts);

  FESpace_T feSpace{*mesh};

  BCList bc{feSpace};

  Var u{"u"};

  interpolateAnalyticFunction([](Vec3 const & p){ return p(0)*p(0); }, feSpace, u.data);

  IOManager io{feSpace, "output_io/u"};
  io.print({u});

  return 0;
}
