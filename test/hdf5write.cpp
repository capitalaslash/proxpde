#include <cstdlib>    // for uint
#include <filesystem> // for path
#include <memory>     // for unique_ptr
#include <optional>   // for optional
#include <string>     // for stoi, basic_string
#include <vector>     // for vector

#include <proxpde/def.hpp>       // for Vec3, Vec
#include <proxpde/fe.hpp>        // for LagrangeFE
#include <proxpde/fespace.hpp>   // for FESpace
#include <proxpde/geo.hpp>       // for Line
#include <proxpde/iomanager.hpp> // for IOManager
#include <proxpde/mesh.hpp>      // for Mesh, buildHyperCube
#include <proxpde/reffe.hpp>     // for RefLineP1
#include <proxpde/var.hpp>       // for Var

int main(int argc, char * argv[])
{
  using namespace proxpde;

  //  HDF5 h5 {"hdf5write.h5"};
  //  size_t size = 10;
  //  Var v("var", size);
  //  v.data = Vec::LinSpaced(size, 3., 7.5);
  //  h5.print(v);

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;
  uint const numElems = (argc < 2) ? 3 : std::stoi(argv[1]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {{numElems, 0, 0}});

  FESpace_T feSpace{*mesh};

  Var sol{"sol", feSpace.dof.size};
  IOManager io{feSpace, "output_io/sol"};

  for (uint itime = 0; itime < 5; ++itime)
  {
    sol.data = Vec::LinSpaced(feSpace.dof.size, itime, itime * 2);
    io.print({sol}, static_cast<double>(itime));
  }

  return 0;
}
