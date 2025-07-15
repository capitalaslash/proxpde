#pragma once

// defines
#include "def.hpp"

// stl
#include <filesystem>

// external libs
#include <hdf5.h>
#include <pugixml.hpp>

// local
#include "fe.hpp"
#include "hdf5.hpp"
#include "l2projection.hpp"
#include "mesh.hpp"
#include "var.hpp"
#include "xdmf_doc.hpp"
#include "xdmf_traits.hpp"

namespace proxpde
{

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

  void print(std::vector<Var> const & data, double const t = 0.0);

  void print(std::vector<FEVar<FESpace_T>> const & data, double const t = 0.0);

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
void IOManager<FESpace>::print(std::vector<Var> const & vars, double const t)
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
  doc.setVar(
      mp,
      {"ptId",
       XDMFType::SCALAR,
       XDMFCenter::NODE,
       XDMFNumberType::INT,
       feSpace->dof.mapSize,
       1u});

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
void IOManager<FESpace>::print(
    std::vector<FEVar<FESpace_T>> const & vars, double const t)
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
  doc.setVar(
      mp,
      {"ptId",
       XDMFType::SCALAR,
       XDMFCenter::NODE,
       XDMFNumberType::INT,
       feSpace->dof.mapSize,
       1u});

  std::filesystem::path mh = filePath.value();
  mh += fmt::format(".{}.h5", iter);
  HDF5 h5Iter{mh, HDF5FileMode::OVERWRITE};
  std::filesystem::path dp = filePath->filename();
  dp += fmt::format(".{}.h5", iter);

  for (auto const & v: vars)
  {
    assert(v.name != "");
    assert(v.data.size() == feSpace->dof.size * FESpace_T::dim);
    using Var_T = std::decay_t<decltype(v)>;
    if constexpr (!std::is_same_v<Var_T, Var>)
    {
      // no need to check the qr
      static_assert(
          std::is_same_v<typename Var_T::FESpace_T::Mesh_T, typename FESpace::Mesh_T>);
      static_assert(
          std::
              is_same_v<typename Var_T::FESpace_T::RefFE_T, typename FESpace::RefFE_T>);
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
  }
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

  Table<double, 3u> coords(feSpace->dof.mapSize, 3u);
  for (auto const & e: feSpace->mesh->elementList)
  {
    auto const elemMappingPts = RefFE_T::mappingPts(e);
    for (auto p = 0u; p < RefFE_T::numGeoDOFs; ++p)
    {
      auto const row = feSpace->dof.geoMap(e.id, p);
      coords.row(row) = elemMappingPts[p];
    }
  }
  h5Mesh.print(coords, "coords");
  Table<DOFid_T, 1u> pointIds(feSpace->dof.mapSize, 1u);
  for (auto k = 0u; k < feSpace->dof.mapSize; k++)
    pointIds[k] = uintNotSet;
  for (auto k = 0u; k < feSpace->mesh->pointList.size(); k++)
    pointIds[feSpace->dof.ptMap[k]] = k;
  h5Mesh.print(pointIds, "ptId");
  // elemIds are the same as IDs in mesh already, no need to print them
  // Table<id_T, 1u> elemIds(feSpace->mesh->elementList.size(), 1u);
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
  doc.setTopology<typename LagrangeFE<Facet_T, 1u>::RefFE_T>(
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
    for (auto p = 0u; p < Facet_T::numPts; ++p)
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

  void print(std::vector<Var> const & data, double const t = 0.0)
  {
    std::vector<Var> dataP0(data.size());
    for (auto i = 0u; i < data.size(); ++i)
    {
      dataP0[i].name = data[i].name + "P0";
      l2Projector.setRhs(data[i].data);
      dataP0[i].data = l2Projector.apply();
    }
    io.print(dataP0, t);
  }

  template <typename FEVar>
  void print(std::vector<FEVar> const & vars, double const t = 0.0)
  {
    std::vector<Var> varsP0;
    for (auto const & v: vars)
    {
      // no need to check the qr is the same
      using Var_T = std::decay_t<decltype(v)>;
      static_assert(std::is_same_v<
                    typename Var_T::FESpace_T::Mesh_T,
                    typename FESpaceOrig_T::Mesh_T>);
      static_assert(std::is_same_v<
                    typename Var_T::FESpace_T::RefFE_T,
                    typename FESpaceOrig_T::RefFE_T>);
      varsP0.emplace_back(v.name + "P0");
      l2Projector.setRhs(v.data);
      varsP0.back().data = l2Projector.apply();
    }

    io.print(varsP0, t);
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

  template <typename FEVar>
  void print(std::vector<FEVar> const & vars, double const t = 0.0)
  {
    std::vector<Var> varsFacet;
    for (auto const & v: vars)
    {
      // no need to check the qr is the same
      using Var_T = std::decay_t<decltype(v)>;
      static_assert(std::is_same_v<
                    typename Var_T::FESpace_T::Mesh_T,
                    typename FESpaceOrig_T::Mesh_T>);
      static_assert(std::is_same_v<
                    typename Var_T::FESpace_T::RefFE_T,
                    typename FESpaceOrig_T::RefFE_T>);
      varsFacet.emplace_back(v.name + "Facet");
      varsFacet.back().data = Vec::Zero(meshFacet->elementList.size());
      for (auto const & facet: v.feSpace->mesh->facetList)
      {
        auto const [insideElemPtr, side] = facet.facingElem[0];
        auto const dofId = v.feSpace->dof.getId(insideElemPtr->id, side);
        varsFacet.back().data[facet.id] = v.data[dofId];
      }
    }
    io.print(varsFacet, t);
  }

  FESpaceOrig_T const * feSpaceOrig;
  std::unique_ptr<MeshFacet_T> meshFacet;
  FESpaceFacet_T feSpaceFacet;
  IOManager<FESpaceFacet_T> io;
};

} // namespace proxpde
