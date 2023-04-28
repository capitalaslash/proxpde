#include "builder.hpp"
#include "def.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      GaussQR<Elem_T, 1>,
      1,
      DofType::DISCONTINUOUS>;

  std::array<uint, 3> numElems = {{2, 1, 0}};
  if (argc == 3)
  {
    numElems[0] = static_cast<uint>(std::stoi(argv[1]));
    numElems[1] = static_cast<uint>(std::stoi(argv[2]));
  }

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{1., 1., 1.}, numElems);

  FESpace_T feSpace{*mesh};

  Var u{"u"};

  interpolateAnalyticFunction(
      [](Vec3 const & p) { return p(0) * p(0); }, feSpace, u.data);

  IOManager io{feSpace, "output_iodisc/u"};
  io.print({u});

  return 0;
}
