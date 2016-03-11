#pragma once

#include "def.hpp"
#include "xdmf_traits.hpp"

#include <fstream>

template <typename Mesh>
struct IOManager
{
  typedef XDMFTraits<typename Mesh::Elem_T> Traits_T;

  explicit IOManager(std::shared_ptr<Mesh> m):
    mesh(m)
  {}

  void print(Vec const& sol);

  std::shared_ptr<Mesh> mesh;
};

// implementation --------------------------------------------------------------

template <typename Mesh>
void IOManager<Mesh>::print(Vec const& sol)
{
  uint const numPts = mesh->pointList.size();
  uint const numElems = mesh->elementList.size();

  //  system("mkdir -p output");
  const std::string fileName = "sol.xmf";
  std::ofstream fout(fileName.c_str());

  fout << "<?xml version=\"1.0\" ?>" << std::endl;
  fout << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
  fout << "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">" << std::endl;
  fout << "  <Domain>" << std::endl;
  fout << "    <Grid GridType=\"Uniform\">" << std::endl;
  fout << "      <Time Value=\"" << 0.0 << "\" />" << std::endl;
  fout << "      <Topology TopologyType=\"" << Traits_T::shape_name << "\" Dimensions=\"" << numElems << "\">" << std::endl;
  //  fout << "        <DataItem Dimensions=\"" << numElems << " " << Elem_T::numPts << "\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\">" << std::endl;
  fout << "        <DataItem Dimensions=\"" << numElems << " " << Mesh::Elem_T::numPts << "\" NumberType=\"Int\" Precision=\"8\" Format=\"XML\">" << std::endl;
  for(auto& e: mesh->elementList)
  {
    for(auto n: e.pointList)
      fout << n->id << " ";
    fout << std::endl;
  }
  //  fout << _data.name << ".mesh.h5:connectivity " << std::endl;
  fout << "        </DataItem>" << std::endl;
  fout << "      </Topology>" << std::endl;
  fout << "      <Geometry GeometryType=\"XYZ\">" << std::endl;
  //  fout << "        <DataItem Dimensions=\"" << numPts << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
  fout << "        <DataItem Dimensions=\"" << numPts << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">" << std::endl;
  for(auto &p: mesh->pointList)
  {
    fout << p.coord[0] << " " << p.coord[1] << " " << p.coord[2] << " " << std::endl;
  }
  //  fout << _data.name << ".mesh.h5:coords " << std::endl;
  fout << "        </DataItem>" << std::endl;
  fout << "      </Geometry>" << std::endl;

  for(uint v=0; v<1; ++v)
  {
    fout << "      <Attribute Name=\"" << "sol" << "\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
    // fout << "        <DataItem Dimensions=\"" << numPts << " 1\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
    fout << "        <DataItem Dimensions=\"" << numPts << " 1\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">" << std::endl;
    for(uint p=0; p<numPts; ++p)
      fout << sol(p) << "\n";
    //      fout << _data.name << "." << step << ".h5:" << _varNames[v] << std::endl;
    fout << "        </DataItem>" << std::endl;
    fout << "      </Attribute>" << std::endl;
  }
  fout << "    </Grid>" << std::endl;
  fout << "  </Domain>" << std::endl;
  fout << "</Xdmf>" << std::endl;
}
