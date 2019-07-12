#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

int main()
{
  using Mesh_T = Mesh<Hexahedron>;
  using FESpace_T = FESpace<Mesh_T, RefHexahedronQ2, GaussQR<Hexahedron,1>>;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{1., 1., 1.}, {1, 1, 1});

  FESpace_T feSpace{*mesh};

  Var u{"u"};

  interpolateAnalyticFunction([](Vec3 const & p){ return p(0)*p(0); }, feSpace, u.data);

  uint const nsteps = 10;
  double const dt = 0.1;
  double time = 0.0;
  IOManager io{feSpace, "output_iotime/u"};
  io.print({u});
  uint const printStep = 2;
  for (uint itime=0; itime<nsteps; ++itime)
  {
    time += dt;
    interpolateAnalyticFunction([time](Vec3 const & p){ return p(0)*p(0) * (1.-time); }, feSpace, u.data);

    if ((itime+1) % printStep == 0)
    {
      io.print({u}, time);
    }
  }

  return 0;
}
