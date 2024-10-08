#pragma once

#include "def.hpp"

#include "fe.hpp"
#include "feutils.hpp"
#include "var.hpp"
#include "xdmf_traits.hpp"

#include <filesystem>

#include <hdf5.h>
#include <pugixml.hpp>

namespace proxpde
{

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

class XDMFDoc
{
public:
  explicit XDMFDoc(
      std::filesystem::path const fp, XDMFGridType gridType = XDMFGridType::SINGLE):
      filePath{fp}
  {
    assert(fp.extension() == ".xmf");

    auto decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version") = "1.0";
    decl.append_attribute("encoding") = "UTF-8";

    doc.append_child(pugi::node_doctype).set_value("Xdmf SYSTEM \"Xdmf.dtd\" []");

    auto xdmfNode = doc.append_child("Xdmf");
    xdmfNode.append_attribute("xmlns:xi") = "http://www.w3.org/2003/XInclude";
    xdmfNode.append_attribute("Version") = "2.2";

    auto domainNode = xdmfNode.append_child("Domain");

    gridNode = domainNode.append_child("Grid");
    gridNode.append_attribute("GridType") = to_string(gridType).data();
    if (gridType == XDMFGridType::COLLECTION)
    {
      gridNode.append_attribute("CollectionType") = "Temporal";
    }
  }

  ~XDMFDoc() { doc.save_file(filePath.c_str()); }

  void setTime(double const time)
  {
    auto timeNode = gridNode.append_child("Time");
    timeNode.append_attribute("Value") = time;
  }

  template <typename RefFE>
  void setTopology(
      std::filesystem::path meshPath,
      uint const numElems,
      std::string_view const connName = "connectivity")
  {
    assert(meshPath.extension() == ".h5");
    auto topoNode = gridNode.append_child("Topology");
    topoNode.append_attribute("TopologyType") = XDMFTraits<RefFE>::shapeName;
    topoNode.append_attribute("Dimensions") = numElems;

    createDataItem(
        topoNode,
        {numElems, RefFE::numGeoDOFs},
        XDMFNumberType::INT,
        8,
        XDMFFormat::HDF,
        meshPath.string() + ":/" + connName.data());
  }

  void setGeometry(
      std::filesystem::path meshPath,
      uint const mapSize,
      std::string_view const coordName = "coords")
  {
    assert(meshPath.extension() == ".h5");
    auto geoNode = gridNode.append_child("Geometry");
    geoNode.append_attribute("GeometryType") = "XYZ";

    createDataItem(
        geoNode,
        {mapSize, 3},
        XDMFNumberType::FLOAT,
        8,
        XDMFFormat::HDF,
        meshPath.string() + ":/" + coordName.data());
  }

  void setVar(std::filesystem::path dataPath, XDMFVar const & var)
  {
    assert(dataPath.extension() == ".h5");
    auto varNode = gridNode.append_child("Attribute");
    varNode.append_attribute("Name") = var.name.data();
    varNode.append_attribute("Active") = 1;
    varNode.append_attribute("AttributeType") = to_string(var.type).data();
    varNode.append_attribute("Center") = to_string(var.center).data();

    // auto const buf = filepath.filename().string() + ".time.h5:/" + name + "." + iter;
    createDataItem(
        varNode,
        {var.size, var.dim},
        var.numberType,
        8,
        XDMFFormat::HDF,
        dataPath.string() + ":/" + var.name);
  }

  void setTimeSeries(
      std::filesystem::path basePath,
      std::vector<std::pair<uint, double>> const & timeSeries)
  {
    // TODO: validate that the timeSeries contains unique time values
    std::string timesStr = "";
    for (auto const & step: timeSeries)
    {
      timesStr += std::to_string(step.second) + " ";
      auto stepNode = gridNode.append_child("xi:include");
      std::filesystem::path p = basePath;
      p += fmt::format(".{}.xmf", step.first);
      stepNode.append_attribute("href") = p.c_str();
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
    node.append_attribute("NumberType") = to_string(type).data();
    node.append_attribute("Precision") = precision;
    node.append_attribute("Format") = to_string(format).data();
    node.text() = ("\n" + content + "\n").data();
    return node;
  }

  std::filesystem::path filePath;
  pugi::xml_document doc;
  pugi::xml_node gridNode;
};

// =====================================================================================

template <typename T>
struct HDF5Var
{};

template <>
struct HDF5Var<int>
{
  static hid_t const value;
  static hid_t const type = H5T_INTEGER;
};

template <>
struct HDF5Var<uint>
{
  static hid_t const value;
  static hid_t const type = H5T_INTEGER;
};

template <>
struct HDF5Var<double>
{
  static hid_t const value;
  static hid_t const type = H5T_FLOAT;
};

enum class HDF5FileMode : int8_t
{
  OVERWRITE,
  APPEND,
  READ,
};

class HDF5
{
public:
  HDF5(std::filesystem::path const fp, HDF5FileMode const mode): filePath{fp}, status(0)
  {
    if (mode == HDF5FileMode::OVERWRITE)
    {
      fileId = H5Fcreate(filePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    else if (mode == HDF5FileMode::APPEND)
    {
      fileId = H5Fopen(filePath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else if (mode == HDF5FileMode::READ)
    {
      fileId = H5Fopen(filePath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    }
    else
    {
      fmt::print(stderr, "mode {} not recognized, aborting!\n", static_cast<int>(mode));
    }
    assert(fileId != H5I_INVALID_HID);
  }

  ~HDF5() { status = H5Fclose(fileId); }

  int readAttribute(
      std::string_view const datasetName, std::string_view const attributeName)
  {
    hid_t const dataset = H5Dopen(fileId, datasetName.data(), H5P_DEFAULT);
    const hid_t attributeId = H5Aopen_name(dataset, attributeName.data());

    int numVals = -1;
    H5Aread(attributeId, H5T_NATIVE_INT, &numVals);
    H5Aclose(attributeId);
    assert(numVals != -1);

    H5Dclose(dataset);

    return numVals;
  }

  template <typename T>
  void readDataset(std::string_view const datasetName, std::vector<T> & buf)
  {
    hid_t const dataset = H5Dopen(fileId, datasetName.data(), H5P_DEFAULT);
    const hid_t attributeId = H5Aopen_name(dataset, "NBR");

    int numVals = -1;
    H5Aread(attributeId, H5T_NATIVE_INT, &numVals);
    H5Aclose(attributeId);
    assert(numVals != -1);

    hid_t const dataType = H5Dget_type(dataset);
    const hid_t type = H5Tget_class(dataType);
    assert(type == HDF5Var<T>::type);

    if constexpr (HDF5Var<T>::type == H5T_FLOAT)
    {
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    }
    else if constexpr (HDF5Var<T>::type == H5T_INTEGER)
    {
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    }
    else
    {
      std::abort();
    }

    H5Dclose(dataset);
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

protected:
  std::filesystem::path filePath;
  hid_t fileId;

public:
  herr_t status;
};

// =====================================================================================

template <typename FESpace>
struct IOManager
{
  static DofType constexpr type = FESpace::type;
  using FESpace_T = FESpace;
  using Mesh_T = typename FESpace_T::Mesh_T;
  using RefFE_T = typename FESpace_T::RefFE_T;
  using FESpaceScalar_T = Scalar_T<FESpace>;
  using Elem_T = typename Mesh_T::Elem_T;
  using Traits_T = XDMFTraits<RefFE_T>;

  IOManager() = default;

  explicit IOManager(
      FESpace_T const & fe,
      std::optional<std::filesystem::path> const fp = std::nullopt,
      uint const it = 0):
      feSpace{&fe},
      feSpaceScalar{*fe.mesh},
      // h5Time{fs::path{filepath} += ".time.h5"},
      iter(it)
  {
    // init only when path is specified
    if (fp.has_value())
    {
      setPath(fp.value());
    }
  }

  void init(FESpace_T const & fe, std::filesystem::path const fp, uint const it = 0)
  {
    feSpace = &fe;
    feSpaceScalar = FESpaceScalar_T{*fe.mesh};
    filePath = fp;
    iter = it;

    setPath(fp);
  }

  void setPath(std::filesystem::path const fp)
  {
    filePath = fp;
    // create subfolder if not saving in the current directory
    if (filePath->parent_path() != std::filesystem::path(""))
    {
      std::filesystem::create_directories(filePath->parent_path());
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
  void printTimeSeries();
  void printMeshData();
  void printBoundary();

public:
  FESpace_T const * feSpace;
  FESpaceScalar_T feSpaceScalar;
  std::optional<std::filesystem::path> filePath = std::filesystem::current_path();
  // HDF5 h5Time;
  uint iter;
  std::vector<std::pair<uint, double>> timeSeries;
};

// implementation ----------------------------------------------------------------------

template <typename FESpace>
void IOManager<FESpace>::print(std::vector<Var> const && vars, double const t)
{
  timeSeries.push_back({iter, t});
  printTimeSeries();

  std::filesystem::path fp = filePath.value();
  fp += fmt::format(".{}.xmf", iter);
  std::filesystem::path dp = filePath->filename();
  dp += fmt::format(".{}.h5", iter);
  std::filesystem::path mp = filePath->filename();
  mp += ".mesh.h5";
  XDMFDoc doc{fp};
  doc.setTime(t);
  doc.setTopology<typename FESpace::RefFE_T>(mp, feSpace->mesh->elementList.size());
  doc.setGeometry(mp, feSpace->dof.mapSize);

  std::filesystem::path mh = filePath.value();
  mh += fmt::format(".{}.h5", iter);
  HDF5 h5Iter{mh, HDF5FileMode::OVERWRITE};

  for (auto const & v: vars)
  {
    // mixed variable vectors can be longer than current fespace
    assert(v.data.size() >= feSpace->dof.size * FESpace_T::dim);

    for (uint d = 0; d < FESpace_T::dim; ++d)
    {
      std::string name;
      if constexpr (FESpace_T::dim == 1)
      {
        name = v.name;
      }
      else if constexpr (FESpace_T::dim > 1)
      {
        name = v.name + "_" + std::to_string(d);
      }

      doc.setVar(
          dp,
          {name,
           XDMFType::SCALAR,
           Traits_T::center,
           XDMFNumberType::FLOAT,
           feSpace->dof.size,
           1});

      // this works only with Lagrange elements
      Vec compdata{feSpace->dof.size};
      // TODO: pass data as const &
      if constexpr (FESpace_T::dim > 1)
      {
        // TODO: print vector variable as vector xdmf data
        getComponent(compdata, feSpaceScalar, v.data, *feSpace, d);
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
  printTimeSeries();

  std::filesystem::path fp = filePath.value();
  fp += fmt::format(".{}.xmf", iter);
  XDMFDoc doc{fp};
  doc.setTime(t);
  std::filesystem::path mp = filePath->filename();
  mp += ".mesh.h5";
  doc.setTopology<typename FESpace::RefFE_T>(mp, feSpace->mesh->elementList.size());
  doc.setGeometry(mp, feSpace->dof.mapSize);

  std::filesystem::path mh = filePath.value();
  mh += fmt::format(".{}.h5", iter);
  HDF5 h5Iter{mh, HDF5FileMode::OVERWRITE};
  std::filesystem::path dp = filePath->filename();
  dp += fmt::format(".{}.h5", iter);

  static_for(
      vars,
      [&](auto const /*i*/, auto const & v)
      {
        assert(v.name != "");
        assert(v.data.size() == feSpace->dof.size * FESpace_T::dim);
        using Var_T = std::decay_t<decltype(v)>;
        if constexpr (!std::is_same_v<Var_T, Var>)
        {
          // no need to check the qr
          static_assert(std::is_same_v<
                        typename Var_T::FESpace_T::Mesh_T,
                        typename FESpace::Mesh_T>);
          static_assert(std::is_same_v<
                        typename Var_T::FESpace_T::RefFE_T,
                        typename FESpace::RefFE_T>);
        }

        constexpr uint dim = getDim<Var_T, FESpace>();
        for (uint d = 0; d < dim; ++d)
        {
          std::string name;
          if constexpr (dim == 1)
          {
            name = v.name;
          }
          else if constexpr (dim > 1)
          {
            name = v.name + "_" + std::to_string(d);
          }

          doc.setVar(
              dp,
              {name,
               XDMFType::SCALAR,
               Traits_T::center,
               XDMFNumberType::FLOAT,
               feSpace->dof.size,
               1});

          // this works only with Lagrange elements
          // TODO: pass data as const &
          if constexpr (dim > 1)
          {
            // TODO: print vector variable as vector xdmf data
            Vec compData{feSpace->dof.size};
            if constexpr (std::is_same_v<Var_T, Var>)
            {
              getComponent(compData, feSpaceScalar, v.data, *feSpace, d);
            }
            else
            {
              getComponent(compData, feSpaceScalar, v.data, *v.feSpace, d);
            }
            h5Iter.print(compData, name);
            // h5Time.print(compdata, name + "." + std::to_string(iter));
          }
          else
          {
            h5Iter.print(v.data, name);
            // h5Time.print(v.data, name + "." + std::to_string(iter));
          }
        }
      });
  iter++;
}

template <typename FESpace>
void IOManager<FESpace>::printTimeSeries()
{
  std::filesystem::path fp = filePath.value();
  fp += ".time.xmf";
  XDMFDoc doc{fp, XDMFGridType::COLLECTION};

  doc.setTimeSeries(filePath->filename(), timeSeries);
}

template <typename FESpace>
void IOManager<FESpace>::printMeshData()
{
  std::filesystem::path fp = filePath.value();
  fp += ".mesh.h5";
  HDF5 h5Mesh{fp, HDF5FileMode::OVERWRITE};

  if constexpr (Traits_T::needsMapping == true)
  {
    Table<DOFid_T, RefFE_T::numGeoDOFs> mappedConn(
        feSpace->mesh->elementList.size(), RefFE_T::numGeoDOFs);
    for (uint i = 0; i < feSpace->mesh->elementList.size(); ++i)
    {
      for (uint j = 0; j < RefFE_T::numGeoDOFs; ++j)
      {
        mappedConn(i, j) = feSpace->dof.geoMap(i, Traits_T::mapping[j]);
      }
    }
    h5Mesh.print(mappedConn, "connectivity");
  }
  else
  {
    h5Mesh.print(feSpace->dof.geoMap, "connectivity");
  }

  Table<double, 3> coords(feSpace->dof.mapSize, 3);
  for (auto const & e: feSpace->mesh->elementList)
  {
    for (uint p = 0; p < RefFE_T::numGeoDOFs; ++p)
    {
      coords.row(feSpace->dof.geoMap(e.id, p)) = RefFE_T::mappingPts(e)[p];
    }
  }
  h5Mesh.print(coords, "coords");
}

template <typename FESpace>
void IOManager<FESpace>::printBoundary()
{
  using Facet_T = typename Mesh_T::Facet_T;
  Mesh_T const & mesh = *feSpace->mesh;

  std::filesystem::path fp = filePath.value();
  fp += ".boundary.xmf";
  std::filesystem::path mp = filePath->filename();
  mp += ".mesh.h5";
  XDMFDoc doc{fp};
  // we always use linear elements for boundary facets
  doc.setTopology<typename LagrangeFE<Facet_T, 1>::RefFE_T>(
      mp, mesh.facetList.size(), "connectivity_bd");
  doc.setGeometry(mp, mesh.pointList.size(), "coords_bd");
  doc.setVar(
      mp,
      {"facetMarker",
       XDMFType::SCALAR,
       XDMFCenter::CELL,
       XDMFNumberType::FLOAT,
       mesh.facetList.size(),
       1});
  doc.setVar(
      mp,
      {"nodeMarker",
       XDMFType::SCALAR,
       XDMFCenter::NODE,
       XDMFNumberType::FLOAT,
       mesh.pointList.size(),
       1});

  std::filesystem::path mh = filePath.value();
  mh += ".mesh.h5";
  HDF5 h5Mesh{mh, HDF5FileMode::APPEND};

  Table<id_T, Facet_T::numPts> conn(mesh.facetList.size(), Facet_T::numPts);
  for (auto const & f: mesh.facetList)
  {
    for (uint p = 0; p < Facet_T::numPts; ++p)
    {
      conn(f.id, p) = f.pts[p]->id;
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

// =====================================================================================

template <typename FESpaceOrig>
struct ToP0
{
  using FESpaceOrig_T = FESpaceOrig;
  using Mesh_T = typename FESpaceOrig_T::Mesh_T;
  using Elem_T = typename Mesh_T::Elem_T;
  using type = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 0>::RefFE_T,
      typename FESpaceOrig_T::QR_T,
      Elem_T::dim>;
};

template <typename FESpaceOrig>
using ToP0_T = typename ToP0<FESpaceOrig>::type;

template <typename T>
struct ToVar
{
  using type = Var;
};

template <typename T>
using ToVar_T = typename ToVar<T>::type;

template <typename FESpaceOrig>
struct IOManagerP0
{
  using FESpaceOrig_T = FESpaceOrig;
  using FESpaceP0_T = ToP0_T<FESpaceOrig>;

  IOManagerP0() = default;
  ~IOManagerP0() = default;

  IOManagerP0(
      FESpaceOrig_T const & fe, std::filesystem::path const fp, uint const it = 0):
      feSpaceOrig{&fe},
      feSpaceP0{*fe.mesh},
      io{feSpaceP0, fp, it},
      l2Projector{feSpaceP0, *feSpaceOrig}
  {}

  void init(FESpaceOrig_T const & fe, std::filesystem::path const fp, uint const it = 0)
  {
    feSpaceOrig->init(fe);
    feSpaceP0.init(*fe.mesh);
    io.init(feSpaceP0, fp, it);
    l2Projector.init(feSpaceP0, *feSpaceOrig);
  }

  void print(std::vector<Var> const && data, double const t = 0.0)
  {
    std::vector<Var> dataP0(data.size());
    for (short_T i = 0; i < data.size(); ++i)
    {
      dataP0[i].name = data[i].name + "P0";
      l2Projector.setRhs(data[i].data);
      dataP0[i].data = l2Projector.apply();
    }
    io.print(dataP0, t);
  }

  template <typename... Vars>
  void print(std::tuple<Vars...> const && data, double const t = 0.0)
  {
    // std::tuple<FEVar<ToP0_T<typename Vars::FESpace_T>...>> dataP0{
    //     FEVar<ToP0_T<typename Vars::FESpace_T>>{"none"}...};
    std::tuple<ToVar_T<Vars>...> dataP0{};
    static_for(
        dataP0,
        data,
        [this](uint const /*i*/, auto & vP0, auto const & v)
        {
          using Var_T = std::decay_t<decltype(v)>;
          if constexpr (!std::is_same_v<Var_T, Var>)
          {
            // we need only the mesh and reffe to be the same, no need to check the
            // qr
            static_assert(std::is_same_v<
                          typename Var_T::FESpace_T::Mesh_T,
                          typename FESpaceOrig_T::Mesh_T>);
            static_assert(std::is_same_v<
                          typename Var_T::FESpace_T::RefFE_T,
                          typename FESpaceOrig_T::RefFE_T>);
          }
          vP0.name = v.name + "P0";
          l2Projector.setRhs(v.data);
          vP0.data = l2Projector.apply();
        });

    io.print(std::move(dataP0), t);
  }

  FESpaceOrig_T const * feSpaceOrig;
  FESpaceP0_T feSpaceP0;
  IOManager<FESpaceP0_T> io;
  L2Projector<FESpaceP0_T, FESpaceOrig_T> l2Projector;
};

// =====================================================================================

template <typename FESpaceOrig>
struct IOManagerFacet
{
  using FESpaceOrig_T = FESpaceOrig;
  using Mesh_T = typename FESpaceOrig_T::Mesh_T;
  using Elem_T = typename Mesh_T::Elem_T;
  using Facet_T = typename Elem_T::Facet_T;
  using MeshFacet_T = Mesh<Facet_T>;
  using FESpaceFacet_T = FESpace<
      MeshFacet_T,
      typename LagrangeFE<Facet_T, 0>::RefFE_T,
      typename LagrangeFE<Facet_T, 0>::RecommendedQR>;

  IOManagerFacet() = default;
  ~IOManagerFacet() = default;

  IOManagerFacet(
      FESpaceOrig_T const & fe, std::filesystem::path const fp, uint const it = 0):
      feSpaceOrig{&fe},
      meshFacet{new MeshFacet_T}
  {
    buildFacetMesh(*meshFacet, *fe.mesh);
    feSpaceFacet.init(*meshFacet);
    io.init(feSpaceFacet, fp, it);
  }

  void init(FESpaceOrig_T const & fe, std::filesystem::path const fp, uint const it = 0)
  {
    feSpaceOrig->init(&fe);
    meshFacet.reset(new MeshFacet_T);
    buildFacetMesh(*meshFacet, *fe.mesh);
    feSpaceFacet.init(*meshFacet);
    io.init(feSpaceFacet, fp, it);
  }

  template <typename... Vars>
  void print(std::tuple<Vars...> const && data, double const t = 0.0)
  {
    std::tuple<ToVar_T<Vars>...> dataFacet{};
    static_for(
        dataFacet,
        data,
        [this](uint const /*i*/, auto & vFacet, auto const & v)
        {
          using Var_T = std::decay_t<decltype(v)>;
          if constexpr (!std::is_same_v<Var_T, Var>)
          {
            // we need only the mesh and reffe to be the same, no need to check the
            // qr
            static_assert(std::is_same_v<
                          typename Var_T::FESpace_T::Mesh_T,
                          typename FESpaceOrig_T::Mesh_T>);
            static_assert(std::is_same_v<
                          typename Var_T::FESpace_T::RefFE_T,
                          typename FESpaceOrig_T::RefFE_T>);
          }
          vFacet.name = v.name + "Facet";
          vFacet.data = Vec::Zero(static_cast<uint>(meshFacet->elementList.size()));
          for (auto const & facet: v.feSpace->mesh->facetList)
          {
            auto const [insideElemPtr, side] = facet.facingElem[0];
            auto const dofId = v.feSpace->dof.getId(insideElemPtr->id, side);
            vFacet.data[facet.id] = v.data[dofId];
          }
        });

    io.print(std::move(dataFacet), t);
  }

  FESpaceOrig_T const * feSpaceOrig;
  std::unique_ptr<MeshFacet_T> meshFacet;
  FESpaceFacet_T feSpaceFacet;
  IOManager<FESpaceFacet_T> io;
};

} // namespace proxpde
