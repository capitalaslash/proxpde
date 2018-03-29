#include "def.hpp"
#include "iomanager.hpp"

#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"

int main(int argc, char* argv[])
{
//  HDF5 h5 {"hdf5write.h5"};
//  size_t size = 10;
//  Var v("var", size);
//  v.data = Vec::LinSpaced(size, 3., 7.5);
//  h5.print(v);

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<Mesh_T,
                            FEType<Elem_T,1>::RefFE_T,
                            FEType<Elem_T,1>::RecommendedQR>;
  uint const numPts = (argc < 2)? 4 : std::stoi(argv[1]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts, 0, 0}});

  FESpace_T feSpace{meshPtr};

  Var sol{"sol", feSpace.dof.size};
  IOManager<FESpace_T> io{feSpace, "output_io/sol"};

  for (auto itime=0; itime<5; ++itime)
  {
    io.iter = itime;
    io.time = itime;
    sol.data = Vec::LinSpaced(feSpace.dof.size, itime, itime*2);
    io.print({sol});
  }

  return 0;
}
