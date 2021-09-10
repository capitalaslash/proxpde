#pragma once

#include "def.hpp"

#include "fe.hpp"
#include "var.hpp"
#include "xdmf_traits.hpp"

#include <experimental/filesystem>

#include <hdf5.h>
#include <pugixml.hpp>

namespace fs = std::experimental::filesystem;

template <typename T>
std::string join(std::vector<T> const & v, std::string_view const divider = " ")
{
  std::string buf;
  for (auto const & i: v)
  {
    buf += std::to_string(i);
    buf += divider;
  }
  return buf.substr(0, buf.size() - 1);
}

template <typename RefFE>
class XDMFDoc
{
public:
  explicit XDMFDoc(
      fs::path const fp,
      std::string_view s = "",
      std::string_view meshS = "mesh",
      std::string_view dataS = "",
      XDMFGridType gridType = XDMFGridType::SINGLE):
      filePath(std::move(fp)),
      suffix(s),
      meshSuffix(meshS),
      dataSuffix(dataS)
  {
    if (!suffix.empty())
    {
      suffix = "." + suffix;
    }
    if (!meshSuffix.empty())
    {
      meshSuffix = "." + meshSuffix;
    }
    if (!dataSuffix.empty())
    {
      dataSuffix = "." + dataSuffix;
    }
    else
    {
      dataSuffix = suffix;
    }

    auto decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version") = "1.0";
    decl.append_attribute("encoding") = "UTF-8";

    doc.append_child(pugi::node_doctype).set_value("Xdmf SYSTEM \"Xdmf.dtd\" []");

    auto xdmfNode = doc.append_child("Xdmf");
    xdmfNode.append_attribute("xmlns:xi") = "http://www.w3.org/2003/XInclude";
    xdmfNode.append_attribute("Version") = "2.2";

    auto domainNode = xdmfNode.append_child("Domain");

    gridNode = domainNode.append_child("Grid");
    gridNode.append_attribute("GridType") = XDMFGridTypeToString.at(gridType).data();
    if (gridType == XDMFGridType::COLLECTION)
    {
      gridNode.append_attribute("CollectionType") = "Temporal";
    }
  }

  ~XDMFDoc()
  {
    auto xmlFilepath = filePath.c_str() + suffix + ".xmf";
    doc.save_file(xmlFilepath.c_str());
  }

  void setTime(double const time)
  {
    auto timeNode = gridNode.append_child("Time");
    timeNode.append_attribute("Value") = time;
  }

  void
  setTopology(uint const numElems, std::string_view const connName = "connectivity")
  {
    auto topoNode = gridNode.append_child("Topology");
    topoNode.append_attribute("TopologyType") = XDMFTraits<RefFE>::shapeName;
    topoNode.append_attribute("Dimensions") = numElems;

    auto const buf =
        filePath.filename().string() + meshSuffix + ".h5:/" + connName.data();
    createDataItem(
        topoNode,
        {numElems, RefFE::numGeoFuns},
        XDMFNumberType::INT,
        8,
        XDMFFormat::HDF,
        buf);
  }

  void setGeometry(uint const mapSize, std::string_view const coordName = "coords")
  {
    auto geoNode = gridNode.append_child("Geometry");
    geoNode.append_attribute("GeometryType") = "XYZ";

    auto const buf =
        filePath.filename().string() + meshSuffix + ".h5:/" + coordName.data();
    createDataItem(
        geoNode, {mapSize, 3}, XDMFNumberType::FLOAT, 8, XDMFFormat::HDF, buf);
  }

  void setVar(XDMFVar const & var)
  {
    auto varNode = gridNode.append_child("Attribute");
    varNode.append_attribute("Name") = var.name.data();
    varNode.append_attribute("Active") = 1;
    varNode.append_attribute("AttributeType") = "Scalar";
    varNode.append_attribute("Center") = XDMFCenterToString.at(var.center).data();

    auto const buf = filePath.filename().string() + dataSuffix + ".h5:/" + var.name;
    // auto const buf = filepath.filename().string() + ".time.h5:/" + name + "." + iter;
    createDataItem(
        varNode, {var.size, 1}, XDMFNumberType::FLOAT, 8, XDMFFormat::HDF, buf);
  }

  void setTimeSeries(std::vector<std::pair<uint, double>> const & timeSeries)
  {
    std::string timesStr = "";
    for (auto const & step: timeSeries)
    {
      timesStr += std::to_string(step.second) + " ";
      auto stepNode = gridNode.append_child("xi:include");
      stepNode.append_attribute("href") =
          (filePath.filename().string() + "." + std::to_string(step.first) + ".xmf")
              .data();
      stepNode.append_attribute("xpointer") = "element(/1/1/1)";
    }
    auto timeNode = gridNode.append_child("Time");
    timeNode.append_attribute("TimeType") = "List";
    createDataItem(
        timeNode,
        {1, timeSeries.size()},
        XDMFNumberType::FLOAT,
        8,
        XDMFFormat::INLINE,
        timesStr);
  }

private:
  pugi::xml_node createDataItem(
      pugi::xml_node & parent,
      std::vector<ulong> const & dims,
      XDMFNumberType const type,
      uint const precision,
      XDMFFormat const format,
      std::string const & content)
  {
    auto node = parent.append_child("DataItem");
    node.append_attribute("Dimensions") = join(dims).c_str();
    node.append_attribute("NumberType") = XDMFNumberTypeToString.at(type);
    node.append_attribute("Precision") = precision;
    node.append_attribute("Format") = XDMFFormatToString.at(format);
    node.text() = ("\n" + content + "\n").data();
    return node;
  }

  fs::path filePath;
  std::string suffix;
  std::string meshSuffix;
  std::string dataSuffix;
  pugi::xml_document doc;
  pugi::xml_node gridNode;
};

template <typename T>
struct HDF5Var
{};

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

enum class HDF5FileMode : int8_t
{
  OVERWRITE,
  APPEND
};

class HDF5
{
public:
  HDF5(fs::path const fp, HDF5FileMode const mode): filePath(std::move(fp)), status(0)
  {
    if (mode == HDF5FileMode::OVERWRITE)
    {
      fileId = H5Fcreate(filePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    else if (mode == HDF5FileMode::APPEND)
    {
      fileId = H5Fopen(filePath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
  }

  // template <typename FESpace>
  // void print(Var const & var)
  // {
  //   print(var.data, var.name);
  //   // hsize_t dimsf[2] = {static_cast<hsize_t>(var.data.size()), 1};
  //   // hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
  //   // hid_t dataset;
  //   // // for(uint v = 0; v < varNames.size(); v++)
  //   // {
  //   //   dataset = H5Dcreate(file_id, var.name.c_str(), H5T_NATIVE_DOUBLE,
  //   //                       dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //   //   status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
  //   //                     H5P_DEFAULT, var.data.data());
  //   //   H5Dclose(dataset);
  //   // }
  //   // H5Sclose(dspace);
  // }

  void print(Vec const & vec, std::string_view const name)
  {
    hsize_t dimsf[2] = {static_cast<hsize_t>(vec.size()), 1};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    dataset = H5Dcreate(
        fileId,
        name.data(),
        H5T_NATIVE_DOUBLE,
        dspace,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT);
    status =
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
    H5Dclose(dataset);
    H5Sclose(dspace);
  }

  template <typename T, unsigned long I>
  void print(Table<T, I> const & tab, std::string_view const name)
  {
    hsize_t dimsf[2] = {
        static_cast<hsize_t>(tab.rows()), static_cast<hsize_t>(tab.cols())};
    hid_t dspace = H5Screate_simple(2, dimsf, nullptr);
    hid_t dataset;
    dataset = H5Dcreate(
        fileId,
        name.data(),
        HDF5Var<T>::value,
        dspace,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT);
    status =
        H5Dwrite(dataset, HDF5Var<T>::value, H5S_ALL, H5S_ALL, H5P_DEFAULT, tab.data());
    H5Dclose(dataset);
    H5Sclose(dspace);
  }

  ~HDF5() { status = H5Fclose(fileId); }

protected:
  fs::path filePath;
  hid_t fileId;
  herr_t status;
};

template <typename FESpace>
struct IOManager
{
  using FESpace_T = FESpace;
  using Mesh_T = typename FESpace_T::Mesh_T;
  using RefFE_T = typename FESpace_T::RefFE_T;
  using FESpaceScalar_T = Scalar_T<FESpace>;
  using Elem_T = typename Mesh_T::Elem_T;
  using Traits_T = XDMFTraits<RefFE_T>;

  IOManager(FESpace_T const & fe, fs::path const fp, uint const it = 0):
      feSpace{fe},
      feSpaceScalar{fe.mesh},
      filePath(std::move(fp)),
      // h5Time{fs::path{filepath} += ".time.h5"},
      iter(it)
  {
    // create subfolder if not saving in the current directory
    if (filePath.parent_path() != fs::path(""))
    {
      fs::create_directories(filePath.parent_path());
    }

    printMeshData();

    if constexpr (Elem_T::dim > 1)
    {
      printBoundary();
    }
  }

  ~IOManager() { printTimeSeries(); }

  void print(std::vector<Var> const && data, double const t = 0.0);

  template <typename VarTup>
  void print(VarTup const && data, double const t = 0.0);

protected:
  void printTimeSeries()
  {
    XDMFDoc<RefFE_T> doc{filePath, "time", "", "", XDMFGridType::COLLECTION};
    doc.setTimeSeries(timeSeries);
  }

  void printMeshData()
  {
    HDF5 h5Mesh{fs::path{filePath} += ".mesh.h5", HDF5FileMode::OVERWRITE};

    if constexpr (Traits_T::needsMapping == true)
    {
      Table<DOFid_T, RefFE_T::numGeoFuns> mappedConn(
          feSpace.mesh.elementList.size(), RefFE_T::numGeoFuns);
      for (uint i = 0; i < feSpace.mesh.elementList.size(); ++i)
      {
        for (uint j = 0; j < RefFE_T::numGeoFuns; ++j)
        {
          mappedConn(i, j) = feSpace.dof.geoMap(i, Traits_T::mapping[j]);
        }
      }
      h5Mesh.print<id_T, RefFE_T::numGeoFuns>(mappedConn, "connectivity");
    }
    else
    {
      h5Mesh.print<id_T, RefFE_T::numGeoFuns>(feSpace.dof.geoMap, "connectivity");
    }

    Table<double, 3> coords(feSpace.dof.mapSize, 3);
    for (auto const & e: feSpace.mesh.elementList)
    {
      for (uint p = 0; p < RefFE_T::numGeoFuns; ++p)
      {
        coords.row(feSpace.dof.geoMap(e.id, p)) = RefFE_T::mappingPts(e)[p];
      }
    }
    h5Mesh.print<double, 3>(coords, "coords");
  }

  void printBoundary()
  {
    using Facet_T = typename Mesh_T::Facet_T;
    Mesh_T const & mesh = feSpace.mesh;

    // we always use linear elements for boundary facets
    XDMFDoc<typename LagrangeFE<Facet_T, 1>::RefFE_T> doc{
        filePath, "boundary", "mesh", "mesh"};
    doc.setTopology(mesh.facetList.size(), "connectivity_bd");
    doc.setGeometry(mesh.pointList.size(), "coords_bd");
    doc.setVar({"facetMarker", XDMFCenter::CELL, mesh.facetList.size()});
    doc.setVar({"nodeMarker", XDMFCenter::NODE, mesh.pointList.size()});

    HDF5 h5Mesh{fs::path{filePath} += ".mesh.h5", HDF5FileMode::APPEND};

    Table<id_T, Facet_T::numPts> conn(mesh.facetList.size(), Facet_T::numPts);
    for (auto const & f: mesh.facetList)
    {
      for (uint p = 0; p < Facet_T::numPts; ++p)
      {
        conn(f.id, p) = f.pointList[p]->id;
      }
    }
    h5Mesh.print<id_T, Facet_T::numPts>(conn, "connectivity_bd");

    Table<double, 3> coords(mesh.pointList.size(), 3);
    for (auto const & p: mesh.pointList)
    {
      coords.row(p.id) = p.coord;
    }
    h5Mesh.print<double, 3>(coords, "coords_bd");

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
  FESpace_T const & feSpace;
  FESpaceScalar_T feSpaceScalar;
  fs::path filePath;
  // HDF5 h5Time;
  uint iter;
  std::vector<std::pair<uint, double>> timeSeries;
};

// implementation --------------------------------------------------------------

template <typename FESpace>
void IOManager<FESpace>::print(std::vector<Var> const && vars, double const t)
{
  timeSeries.push_back({iter, t});
  XDMFDoc<typename FESpace::RefFE_T> doc{filePath, std::to_string(iter)};
  doc.setTime(t);
  doc.setTopology(feSpace.mesh.elementList.size());
  doc.setGeometry(feSpace.dof.mapSize);

  HDF5 h5Iter{
      filePath.string() + "." + std::to_string(iter) + ".h5", HDF5FileMode::OVERWRITE};
  for (auto const & v: vars)
  {
    // mixed variable vectors can be longer than current fespace
    assert(v.data.size() >= feSpace.dof.size * FESpace_T::dim);

    for (uint d = 0; d < FESpace_T::dim; ++d)
    {
      auto name = v.name;
      if constexpr (FESpace_T::dim > 1)
      {
        name += "_" + std::to_string(d);
      }

      doc.setVar({name, Traits_T::attributeType, feSpace.dof.size});

      // this works only with Lagrange elements
      Vec compdata{feSpace.dof.size};
      // TODO: pass data as const &
      if constexpr (FESpace_T::dim > 1)
      {
        // TODO: print vector variable as vector xdmf data
        getComponent(compdata, feSpaceScalar, v.data, feSpace, d);
      }
      else
      {
        compdata = v.data;
      }
      h5Iter.print(compdata, name);
      // h5Time.print(compdata, name + "." + std::to_string(iter));
    }
  }
  iter++;
}

template <typename VarT, typename FESpace>
inline constexpr uint getDim()
{
  if constexpr (std::is_same_v<VarT, Var>)
  {
    return FESpace::dim;
  }
  else
  {
    return VarT::FESpace_T::dim;
  }
}

template <typename FESpace>
template <typename VarTup>
void IOManager<FESpace>::print(VarTup const && vars, double const t)
{
  timeSeries.push_back({iter, t});
  XDMFDoc<typename FESpace::RefFE_T> doc{filePath, std::to_string(iter)};
  doc.setTime(t);
  doc.setTopology(feSpace.mesh.elementList.size());
  doc.setGeometry(feSpace.dof.mapSize);

  HDF5 h5Iter{
      filePath.string() + "." + std::to_string(iter) + ".h5", HDF5FileMode::OVERWRITE};
  static_for(
      vars,
      [&](auto const /*i*/, auto const & v)
      {
        assert(v.name != "");
        assert(v.data.size() == feSpace.dof.size * FESpace_T::dim);
        using Var_T = std::decay_t<decltype(v)>;
        constexpr uint dim = getDim<Var_T, FESpace>();
        if constexpr (!std::is_same_v<Var_T, Var>)
        {
          // we need only the mesh and reffe to be the same, no need to check the
          // qr
          static_assert(std::is_same_v<
                        typename Var_T::FESpace_T::Mesh_T,
                        typename FESpace::Mesh_T>);
          static_assert(std::is_same_v<
                        typename Var_T::FESpace_T::RefFE_T,
                        typename FESpace::RefFE_T>);
        }

        for (uint d = 0; d < dim; ++d)
        {
          auto name = v.name;
          if constexpr (dim > 1)
          {
            name += "_" + std::to_string(d);
          }

          doc.setVar({name, Traits_T::attributeType, feSpace.dof.size});

          // this works only with Lagrange elements
          Vec compData{feSpace.dof.size};
          // TODO: pass data as const &
          if constexpr (dim > 1)
          {
            // TODO: print vector variable as vector xdmf data
            if constexpr (std::is_same_v<Var_T, Var>)
            {
              getComponent(compData, feSpaceScalar, v.data, feSpace, d);
            }
            else
            {
              getComponent(compData, feSpaceScalar, v.data, v.feSpace, d);
            }
          }
          else
          {
            compData = v.data;
          }
          h5Iter.print(compData, name);
          // h5Time.print(compdata, name + "." + std::to_string(iter));
        }
      });
  iter++;
}
