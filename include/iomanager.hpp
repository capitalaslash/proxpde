#pragma once

#include "def.hpp"
#include "xdmf_traits.hpp"
#include "var.hpp"

#include <tinyxml2.h>
#include <hdf5.h>
#include <fstream>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

template <typename T, unsigned long I>
using Table = Eigen::Matrix<T,Eigen::Dynamic,I,Eigen::RowMajor>;

template<typename T>
struct HDF5Var
{
  static hid_t constexpr type = 0;
};

template <>
struct HDF5Var<uint>
{
  static hid_t type;
};
hid_t HDF5Var<uint>::type = H5T_STD_I32LE;

template <>
struct HDF5Var<double>
{
  static hid_t type;
};
hid_t HDF5Var<double>::type = H5T_IEEE_F64LE;

class HDF5
{
public:

  HDF5(fs::path const & fn):
    filepath(fn),
    file_id(H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)),
    status(0)
  {}

  void print(Var const & var)
  {
    hsize_t dimsf[2] = {1, static_cast<hsize_t>(var.data.size())};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    // for(uint v = 0; v < varNames.size(); v++)
    {
      dataset = H5Dcreate(file_id, var.name.c_str(), H5T_NATIVE_DOUBLE,
                          dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, var.data.data());
      H5Dclose(dataset);
    }
    H5Sclose(dspace);
  }

  void print(Vec const & vec, std::string const name)
  {
    hsize_t dimsf[2] = {1, static_cast<hsize_t>(vec.size())};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    // for(uint v = 0; v < varNames.size(); v++)
    {
      dataset = H5Dcreate(file_id, name.c_str(), H5T_NATIVE_DOUBLE,
                          dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, vec.data());
      H5Dclose(dataset);
    }
    H5Sclose(dspace);
  }

  template <typename T, unsigned long I>
  void print(Table<T, I> const & tab, std::string const name)
  {
    hsize_t dimsf[2] = {static_cast<hsize_t>(tab.rows()), static_cast<hsize_t>(tab.cols())};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    dataset = H5Dcreate(file_id, name.c_str(), HDF5Var<T>::type,
                        dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, HDF5Var<T>::type, H5S_ALL, H5S_ALL,H5P_DEFAULT, tab.data());
    H5Dclose(dataset);
    H5Sclose(dspace);
  }

  ~HDF5()
  {
    status = H5Fclose(file_id);
  }

protected:
  fs::path filepath;
  hid_t file_id;
  herr_t status;
};

template <typename FESpace>
struct IOManager
{
  using Elem_T = typename FESpace::Mesh_T::Elem_T;
  using MeshPtr_T = std::shared_ptr<typename FESpace::Mesh_T>;
  using Traits_T = XDMFTraits<typename FESpace::RefFE_T>;

  IOManager(FESpace const & fes,
            fs::path const & fn,
            double tin = 0.0,
            uint it = 0):
    feSpace(fes),
    filepath(fn),
    // h5Time{fs::path{filepath} += ".time.h5"},
    time(tin),
    iter(it)
  {
    fs::create_directory(filepath.parent_path());
    _printMeshData();
  }

  void print(std::vector<Var> const & data);

protected:
  void _printMeshData()
  {
    typename FESpace::Mesh_T const & mesh = *(feSpace.meshPtr);
    HDF5 h5Mesh{fs::path{filepath} += ".mesh.h5"};
    Table<id_T, FESpace::RefFE_T::numFuns> conn(mesh.elementList.size(), FESpace::RefFE_T::numFuns);
    for (auto const & e: mesh.elementList)
      for (uint p=0; p<FESpace::RefFE_T::numFuns; ++p)
        conn(e.id, p) = feSpace.dof.elemMap[e.id][p];
    h5Mesh.print<id_T, FESpace::RefFE_T::numFuns>(conn, "connectivity");

    Table<double, 3> coords(feSpace.dof.totalNum, 3);
    for (auto const & e: mesh.elementList)
    {
      auto localPts = FESpace::RefFE_T::dofPts(e);
      for (uint p=0; p<localPts.size(); ++p)
      {
        coords.row(feSpace.dof.elemMap[e.id][p]) = localPts[p];
      }
    }
    h5Mesh.print<double, 3>(coords, "coords");
  }

public:
  FESpace const & feSpace;
  fs::path filepath;
  // HDF5 h5Time;
  double time;
  uint iter;
};

// implementation --------------------------------------------------------------

template <typename FESpace>
void IOManager<FESpace>::print(std::vector<Var> const & data)
{
  typename FESpace::Mesh_T const & mesh = *(feSpace.meshPtr);
  uint const numPts = feSpace.dof.totalNum;
  uint const numElems = mesh.elementList.size();

  tinyxml2::XMLDocument doc;
  doc.InsertEndChild(doc.NewDeclaration());
  doc.InsertEndChild(doc.NewUnknown("DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []"));

  auto xdmf_el = doc.NewElement("Xdmf");
  xdmf_el->SetAttribute("xmlns:xi", "http://www.w3.org/2003/XInclude");
  xdmf_el->SetAttribute("Version", "2.2");
  auto xdmf = doc.InsertEndChild(xdmf_el);

  auto domain = xdmf->InsertEndChild(doc.NewElement("Domain"));

  auto grid_el = doc.NewElement("Grid");
  grid_el->SetAttribute("GridType", "Uniform");
  auto grid = domain->InsertEndChild(grid_el);

  auto time_el = doc.NewElement("Time");
  time_el->SetAttribute("Value", time);
  grid->InsertEndChild(time_el);

  auto topo_el = doc.NewElement("Topology");
  topo_el->SetAttribute("TopologyType", Traits_T::shape_name);
  topo_el->SetAttribute("Dimensions", numElems);
  auto topo = grid->InsertEndChild(topo_el);

  auto topodata_el = doc.NewElement("DataItem");
  topodata_el->SetAttribute("Dimensions",
    (std::to_string(numElems) + " " + std::to_string(FESpace::RefFE_T::numFuns)).c_str());
  topodata_el->SetAttribute("NumberType", "Int");
  topodata_el->SetAttribute("Precision", 8);
  topodata_el->SetAttribute("Format", "HDF");
  std::stringstream buf;
  buf << std::endl;
  buf << filepath.filename().string() << ".mesh.h5:/connectivity " << std::endl;
  topodata_el->SetText(buf.str().c_str());
  topo->InsertEndChild(topodata_el);

  auto geometry_el = doc.NewElement("Geometry");
  geometry_el->SetAttribute("GeometryType", "XYZ");
  auto geometry = grid->InsertEndChild(geometry_el);

  auto geodata_el = doc.NewElement("DataItem");
  geodata_el->SetAttribute("Dimensions",
    (std::to_string(numPts) + " 3").c_str());
  geodata_el->SetAttribute("NumberType", "Float");
  geodata_el->SetAttribute("Precision", 8);
  geodata_el->SetAttribute("Format", "HDF");
  buf.str("");
  buf << std::endl;
  buf << filepath.filename().string() << ".mesh.h5:/coords " << std::endl;
  geodata_el->SetText(buf.str().c_str());
  geometry->InsertEndChild(geodata_el);

  HDF5 h5Iter{filepath.string() + "." + std::to_string(iter) + ".h5"};
  for(auto& v: data)
  {
    for (uint d=0; d<feSpace.dim; ++d)
    {
      std::string name = v.name;
      if (feSpace.dim > 1)
      {
        name += "_" + std::to_string(d);
      }

      auto var_el = doc.NewElement("Attribute");
      var_el->SetAttribute("Name", name.c_str());
      var_el->SetAttribute("Active", 1);
      var_el->SetAttribute("AttributeType", "Scalar");
      var_el->SetAttribute("Center", "Node");
      auto var = grid->InsertEndChild(var_el);

      auto vardata_el = doc.NewElement("DataItem");
      vardata_el->SetAttribute("Dimensions",
                               ("1 " + std::to_string(numPts)).c_str());
      vardata_el->SetAttribute("NumberType", "Float");
      vardata_el->SetAttribute("Precision", 8);
      vardata_el->SetAttribute("Format", "HDF");
      buf.str("");
      buf << std::endl;
      buf << filepath.filename().string() << "." << iter << ".h5:/" << name << std::endl;
      // buf << filepath.filename().string() << ".time.h5:/" << name << "." << iter << std::endl;
      vardata_el->SetText(buf.str().c_str());
      var->InsertEndChild(vardata_el);
      Vec const compdata = v.data.block(d*feSpace.dof.totalNum, 0, feSpace.dof.totalNum, 1);
      h5Iter.print(compdata, name);
      // h5Time.print(compdata, name + "." + std::to_string(iter));
    }
  }

  auto xmlFilepath = filepath;
  xmlFilepath += std::string(".") + std::to_string(iter) + ".xmf";
  doc.SaveFile(xmlFilepath.c_str());
}
