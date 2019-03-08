#pragma once

#include "def.hpp"
#include "xdmf_traits.hpp"
#include "fe.hpp"
#include "var.hpp"

#include <pugixml.hpp>
#include <hdf5.h>
#include <fstream>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

template <typename T>
std::string join(std::vector<T> const & v)
{
  std::string buf;
  for (auto const & i: v)
  {
    buf += " " + std::to_string(i);
  }
  buf = buf.substr(1, buf.size()-1);
  return buf;
}

enum class XDMFNumberType : int8_t
{
    INT,
    FLOAT
};

static const std::map<XDMFNumberType, char const *> XDMFNumberTypeToString =
{
    {XDMFNumberType::INT, "Int"},
    {XDMFNumberType::FLOAT, "Float"},
};

enum class XDMFFormat : int8_t
{
    HDF,
    //INLINE
};

static const std::map<XDMFFormat, char const *> XDMFFormatToString =
{
    {XDMFFormat::HDF, "HDF"},
};

template <typename RefFE>
class XDMFDoc
{
public:
  explicit XDMFDoc(fs::path const & fp,
                   std::string_view s,
                   std::string_view geoS = "mesh"):
    filepath(fp),
    suffix(s),
    geoSuffix(geoS)
  {
    auto decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version") = "1.0";
    decl.append_attribute("encoding") = "UTF-8";

    doc.append_child(pugi::node_doctype).set_value("Xdmf SYSTEM \"Xdmf.dtd\" []");

    auto xdmfNode = doc.append_child("Xdmf");
    xdmfNode.append_attribute("xmlns:xi") = "http://www.w3.org/2003/XInclude";
    xdmfNode.append_attribute("Version") = "2.2";

    auto domainNode = xdmfNode.append_child("Domain");

    gridNode = domainNode.append_child("Grid");
    gridNode.append_attribute("GridType") = "Uniform";
  }

  ~XDMFDoc()
  {
    auto xmlFilepath = filepath;
    if (!suffix.empty())
    {
      xmlFilepath += "." + suffix;
    }
    xmlFilepath += ".xmf";
    doc.save_file(xmlFilepath.c_str());
  }

  void setTime(double const time)
  {
    auto timeNode = gridNode.append_child("Time");
    timeNode.append_attribute("Value") = time;
  }

  void setTopology(uint const numElems)
  {
    auto topoNode = gridNode.append_child("Topology");
    topoNode.append_attribute("TopologyType") = XDMFTraits<RefFE>::shapeName;
    topoNode.append_attribute("Dimensions") = numElems;

    auto const buf = "\n" + filepath.filename().string() + "." + geoSuffix
        + ".h5:/connectivity\n";
    createDataItem(
            topoNode,
            {numElems, RefFE::numGeoFuns},
            XDMFNumberType::INT,
            8,
            XDMFFormat::HDF,
            buf);
  }

  void setGeometry(uint const mapSize)
  {
    auto geoNode = gridNode.append_child("Geometry");
    geoNode.append_attribute("GeometryType") = "XYZ";

    auto const buf = "\n" + filepath.filename().string() + "." + geoSuffix
        + ".h5:/coords\n";
    createDataItem(
            geoNode,
            {mapSize, 3},
            XDMFNumberType::FLOAT,
            8,
            XDMFFormat::HDF,
            buf);
  }

  // TODO: set up struct with all relevant data (name, type, size, format, number type)
  void setVar(std::string_view const name, std::string_view const type, uint const size)
  {
    auto varNode = gridNode.append_child("Attribute");
    varNode.append_attribute("Name") = name.data();
    varNode.append_attribute("Active") = 1;
    varNode.append_attribute("AttributeType") = "Scalar";
    varNode.append_attribute("Center") = type.data();

    auto const buf = "\n" + filepath.filename().string() + "." + suffix
        + ".h5:/" + name.data() + "\n";
    // auto const buf = "\n" + filepath.filename().string() + ".time.h5:/"
    //     + name + "." + iter + "\n";
    createDataItem(
            varNode,
            {1, size},
            XDMFNumberType::FLOAT,
            8,
            XDMFFormat::HDF,
            buf);
  }

private:
  pugi::xml_node createDataItem(
          pugi::xml_node & parent,
          std::vector<uint> const & dims,
          XDMFNumberType const type,
          uint const precision,
          XDMFFormat const format,
          std::string_view const content)
  {
    auto node = parent.append_child("DataItem");
    node.append_attribute("Dimensions") = join(dims).c_str();
    node.append_attribute("NumberType") = XDMFNumberTypeToString.at(type);
    node.append_attribute("Precision") = precision;
    node.append_attribute("Format") = XDMFFormatToString.at(format);
    node.text() = content.data();
    return node;
  }

  fs::path filepath;
  std::string const suffix;
  std::string const geoSuffix;
  pugi::xml_document doc;
  pugi::xml_node gridNode;
};

template<typename T>
struct HDF5Var {};

template <>
struct HDF5Var<uint>
{
  static hid_t value;
};
hid_t HDF5Var<uint>::value = H5T_STD_I32LE;

template <>
struct HDF5Var<double>
{
  static hid_t value;
};
hid_t HDF5Var<double>::value = H5T_IEEE_F64LE;

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
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, var.data.data());
      H5Dclose(dataset);
    }
    H5Sclose(dspace);
  }

  void print(Vec const & vec, std::string_view const name)
  {
    hsize_t dimsf[2] = {1, static_cast<hsize_t>(vec.size())};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    // for(uint v = 0; v < varNames.size(); v++)
    {
      dataset = H5Dcreate(file_id, name.data(), H5T_NATIVE_DOUBLE,
                          dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, vec.data());
      H5Dclose(dataset);
    }
    H5Sclose(dspace);
  }

  template <typename T, unsigned long I>
  void print(Table<T, I> const & tab, std::string_view const name)
  {
    hsize_t dimsf[2] = {static_cast<hsize_t>(tab.rows()), static_cast<hsize_t>(tab.cols())};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    dataset = H5Dcreate(file_id, name.data(), HDF5Var<T>::value,
                        dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, HDF5Var<T>::value, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, tab.data());
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
  using FESpace_T = FESpace;
  using Mesh_T = typename FESpace::Mesh_T;
  using Elem_T = typename Mesh_T::Elem_T;
  using Traits_T = XDMFTraits<typename FESpace::RefFE_T>;

  IOManager(FESpace_T const & fe,
            fs::path const & fp,
            double const tin = 0.0,
            uint const it = 0):
    feSpace(fe),
    filepath(fp),
    // h5Time{fs::path{filepath} += ".time.h5"},
    time(tin),
    iter(it)
  {
    if (filepath.parent_path() != fs::path(""))
      fs::create_directory(filepath.parent_path());
    printMeshData();
    if constexpr (Elem_T::dim > 1)
    {
      printBoundary();
    }
  }

  void print(std::vector<Var> const & data);

protected:
  void printMeshData()
  {
    using FE_T = typename FESpace_T::RefFE_T;
    Mesh_T const & mesh = feSpace.mesh;
    HDF5 h5Mesh{fs::path{filepath} += ".mesh.h5"};

    if constexpr (Traits_T::needsMapping == true)
    {
      Table<DOFid_T, FE_T::numGeoFuns> mappedConn(feSpace.mesh.elementList.size(), FE_T::numGeoFuns);
      for (uint i=0; i<feSpace.mesh.elementList.size(); ++i)
      {
        for (uint j=0; j<FE_T::numGeoFuns; ++j)
        {
          mappedConn(i, j) = feSpace.dof.geoMap(i, Traits_T::mapping[j]);
        }
      }
      h5Mesh.print<id_T, FE_T::numGeoFuns>(mappedConn, "connectivity");
    }
    else
    {
      h5Mesh.print<id_T, FE_T::numGeoFuns>(feSpace.dof.geoMap, "connectivity");
    }

    Table<double, 3> coords(feSpace.dof.mapSize, 3);
    for (auto const & e: mesh.elementList)
    {
      for (uint p=0; p<FE_T::numGeoFuns; ++p)
      {
        coords.row(feSpace.dof.geoMap(e.id, p)) = FE_T::mappingPts(e)[p];
      }
    }
    h5Mesh.print<double, 3>(coords, "coords");
  }

  void printBoundary()
  {
    using Facet_T = typename Mesh_T::Facet_T;
    Mesh_T const & mesh = feSpace.mesh;

    // we always use linear elements for boundary facets
    XDMFDoc<typename FEType<Facet_T, 1>::RefFE_T> doc{filepath, "meshb", "meshb"};
    doc.setTopology(mesh.facetList.size());
    doc.setGeometry(mesh.pointList.size());
    doc.setVar("facetMarker", "Cell", mesh.facetList.size());
    doc.setVar("nodeMarker", "Node", mesh.pointList.size());

    HDF5 h5Mesh{fs::path{filepath} += ".meshb.h5"};

    Table<id_T, Facet_T::numPts> conn(mesh.facetList.size(), Facet_T::numPts);
    for (auto const & f: mesh.facetList)
    {
      for (uint p=0; p<Facet_T::numPts; ++p)
      {
        conn(f.id, p) = f.pointList[p]->id;
      }
    }
    h5Mesh.print<id_T, Facet_T::numPts>(conn, "connectivity");

    Table<double, 3> coords(mesh.pointList.size(), 3);
    for (auto const & p: mesh.pointList)
    {
      coords.row(p.id) = p.coord;
    }
    h5Mesh.print<double, 3>(coords, "coords");

    Vec facetMarkerData(mesh.facetList.size());
    for (auto const & f: mesh.facetList)
    {
      facetMarkerData[f.id] = static_cast<double>(f.marker);
    }
    h5Mesh.print(facetMarkerData, "facetMarker");

    Vec nodeMarkerData(mesh.pointList.size());
    for (auto const & p: mesh.pointList)
    {
      nodeMarkerData[p.id] = static_cast<double>(p.marker);
    }
    h5Mesh.print(nodeMarkerData, "nodeMarker");
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
  XDMFDoc<typename FESpace::RefFE_T> doc{filepath, std::to_string(iter)};
  doc.setTime(time);
  doc.setTopology(feSpace.mesh.elementList.size());
  doc.setGeometry(feSpace.dof.mapSize);

  HDF5 h5Iter{filepath.string() + "." + std::to_string(iter) + ".h5"};
  for(auto& v: data)
  {
    for (uint d=0; d<feSpace.dim; ++d)
    {
      auto name = v.name;
      if constexpr (FESpace_T::dim > 1)
      {
        name += "_" + std::to_string(d);
      }

      doc.setVar(name, Traits_T::attributeType, feSpace.dof.size);

      // this works only with Lagrange elements
      Vec const compdata = v.data.block(d*feSpace.dof.size, 0, feSpace.dof.size, 1);
      h5Iter.print(compdata, name);
      // h5Time.print(compdata, name + "." + std::to_string(iter));
    }
  }
}
