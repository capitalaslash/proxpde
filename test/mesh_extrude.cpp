#include "def.hpp"

#include "fe.hpp"
#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(int /*argc*/, char * /*argv*/[])
{
  using namespace proxpde;

  using Elem2D_T = Triangle;
  using Mesh2D_T = Mesh<Elem2D_T>;
  using Elem3D_T = Tetrahedron;
  using Mesh3D_T = Mesh<Elem3D_T>;
  using FESpace3D_T = FESpace<
      Mesh3D_T,
      LagrangeFE<Elem3D_T, 0>::RefFE_T,
      LagrangeFE<Elem3D_T, 0>::RecommendedQR>;

  MilliTimer t;

  t.start("mesh2d");
  Mesh2D_T mesh2d;
  buildHyperCube(mesh2d, {0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {3U, 2U, 0U});
  t.stop();

  t.start("extrude");
  Mesh3D_T mesh3d;
  extrude(mesh2d, mesh3d, 4U, Vec3{0.0, 0.0, 1.0}, 2.0);
  t.stop();

  t.start("io");
  FESpace3D_T feSpace3d{mesh3d};
  IOManager io3d{feSpace3d, "output_extrude/mesh3d"};
  FEVar id3d{"id", feSpace3d};
  for (auto const & e: mesh3d.elementList)
  {
    id3d.data[e.id] = e.id;
  }
  io3d.print({id3d});
  t.stop();

  t.print();

  return 0;
}
