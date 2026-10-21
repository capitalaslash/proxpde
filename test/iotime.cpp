#include <cmath>      // for sqrt
#include <cstdlib>    // for uint
#include <filesystem> // for path
#include <memory>     // for unique_ptr
#include <optional>   // for optional
#include <vector>     // for vector

#include <proxpde/def.hpp>       // for Vec3
#include <proxpde/fespace.hpp>   // for interpolateAnalyticFunction, FESpace
#include <proxpde/geo.hpp>       // for Hexahedron
#include <proxpde/iomanager.hpp> // for IOManager
#include <proxpde/mesh.hpp>      // for Mesh, buildHyperCube
#include <proxpde/qr.hpp>        // for GaussQR
#include <proxpde/reffe.hpp>     // for RefHexahedronQ2
#include <proxpde/var.hpp>       // for Var

int main()
{
  using namespace proxpde;

  using Mesh_T = Mesh<Hexahedron>;
  using FESpace_T = FESpace<Mesh_T, RefHexahedronQ2, GaussQR<Hexahedron, 1>>;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{1., 1., 1.}, {1, 1, 1});

  FESpace_T feSpace{*mesh};

  Var u{"u"};

  interpolateAnalyticFunction(
      [](Vec3 const & p) { return p(0) * p(0); }, feSpace, u.data);

  uint const nsteps = 10;
  double const dt = 0.1;
  double time = 0.0;
  IOManager io{feSpace, "output_iotime/u"};
  io.print({u});
  uint const printStep = 2;
  for (uint itime = 0; itime < nsteps; ++itime)
  {
    time += dt;
    interpolateAnalyticFunction(
        [time](Vec3 const & p) { return p(0) * p(0) * (1. - time); }, feSpace, u.data);

    if ((itime + 1) % printStep == 0)
    {
      io.print({u}, time);
    }
  }

  return 0;
}
