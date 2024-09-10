#pragma once

#include "def.hpp"

#include "geo.hpp"

namespace proxpde
{

// =====================================================================
enum class MeshFlags : uint8_t
{
  NONE = 0U,
  BOUNDARY_FACETS = 0x1 << 0,
  INTERNAL_FACETS = 0x1 << 1,
  NORMALS = 0x1 << 2,
  FACET_PTRS = 0x1 << 3,
  CHECK_FACET_PLANAR = 0x1 << 4,
};

static constexpr std::array<MeshFlags, 5> MeshFlags_ALL = {
    MeshFlags::BOUNDARY_FACETS,
    MeshFlags::INTERNAL_FACETS,
    MeshFlags::NORMALS,
    MeshFlags::FACET_PTRS,
    MeshFlags::CHECK_FACET_PLANAR,
};

template <>
struct enable_bitmask_operators<MeshFlags>
{
  static const bool value = true;
};

static constexpr const char * to_string(MeshFlags const flag)
{
  switch (flag)
  {
  case MeshFlags::NONE:
    return "NONE";
  case MeshFlags::BOUNDARY_FACETS:
    return "BOUNDARY_FACETS";
  case MeshFlags::INTERNAL_FACETS:
    return "INTERNAL_FACETS";
  case MeshFlags::NORMALS:
    return "NORMALS";
  case MeshFlags::FACET_PTRS:
    return "FACET_PTRS";
  case MeshFlags::CHECK_FACET_PLANAR:
    return "CHECK_FACET_PLANAR";
  }
  abort();
  return "ERROR";
}

static constexpr MeshFlags to_MeshFlags(std::string_view str)
{
  if (str == "BOUNDARY_FACETS")
    return MeshFlags::BOUNDARY_FACETS;
  else if (str == "INTERNAL_FACETS")
    return MeshFlags::INTERNAL_FACETS;
  else if (str == "NORMALS")
    return MeshFlags::NORMALS;
  else if (str == "FACET_PTRS")
    return MeshFlags::FACET_PTRS;
  else if (str == "CHECK_FACET_PLANAR")
    return MeshFlags::CHECK_FACET_PLANAR;
  fmt::print(stderr, "mesh flag not recognized: {}\n", str);
  abort();
  return MeshFlags::NONE;
}

// =====================================================================
template <typename Elem>
struct Mesh
{
  using Elem_T = Elem;
  using Facet_T = typename Elem::Facet_T;
  using PointList_T = std::vector<Point>;
  using ElementList_T = std::vector<Elem>;
  using FacetList_T = std::vector<Facet_T>;
  using elemToPoint_T = std::vector<std::array<id_T, Elem::numPts>>;
  using elemToFacet_T = std::vector<std::array<id_T, Elem::numFacets>>;

  Mesh() = default;
  Mesh(Mesh<Elem> const &) = default;
  Mesh<Elem> & operator=(Mesh<Elem> &) = default;

  void buildConnectivity();
  uint countRidges();

  PointList_T pointList;
  ElementList_T elementList;
  FacetList_T facetList;
  uint numBdFacets;
  elemToPoint_T elemToPoint;
  elemToFacet_T elemToFacet;
  Bitmask<MeshFlags> flags;
};

template <typename Elem>
void Mesh<Elem>::buildConnectivity()
{
  elemToPoint.reserve(elementList.size());
  for (auto const & e: elementList)
  {
    std::array<id_T, Elem::numPts> elemConn;
    uint counter = 0;
    for (auto const & p: e.pts)
    {
      elemConn[counter] = p->id;
      counter++;
    }
    elemToPoint.push_back(elemConn);
  }
}

template <typename Elem>
uint Mesh<Elem>::countRidges()
{
  // ridges are relevant only in 3d
  if (Elem_T::dim < 3)
  {
    return 0;
  }

  // identify a ridge by its points (unordered)
  using Ridge_T = std::set<id_T>;
  using RidgeList_T = std::set<Ridge_T>;

  RidgeList_T ridges;
  auto ridgeCount = 0U;
  for (auto const & elem: elementList)
  {
    for (auto const ridgeIds: Elem_T::elemToRidge)
    {
      Ridge_T ridge;
      ridge.insert(elem.pts[ridgeIds[0]]->id);
      ridge.insert(elem.pts[ridgeIds[1]]->id);

      [[maybe_unused]] auto const [ptr, inserted] = ridges.insert(ridge);
      if (inserted)
      {
        ridgeCount++;
      }
    }
  }
  return ridgeCount;
}

template <typename Elem>
std::ostream & operator<<(std::ostream & out, Mesh<Elem> const & mesh)
{
  out << "point list\n----------" << std::endl;
  for (auto & p: mesh.pointList)
    out << p << std::endl;
  out << "\nelement list\n------------" << std::endl;
  for (auto & e: mesh.elementList)
    out << e << std::endl;
  out << "\nfacet list\n------------" << std::endl;
  for (auto & f: mesh.facetList)
    out << f << std::endl;
  out << "\nelemToFacet map\n------------" << std::endl;
  for (auto & row: mesh.elemToFacet)
  {
    for (auto & clm: row)
    {
      out << clm << " ";
    }
    out << std::endl;
  }
  out << "\n------------" << std::endl;
  return out;
}

// =====================================================================
enum side : marker_T
{
  BOTTOM = 1,
  RIGHT = 2,
  TOP = 3,
  LEFT = 4,
  BACK = 5,
  FRONT = 6,
  CIRCLE = 101,
};

template <uint dim>
struct AllSides
{};

template <>
struct AllSides<1U>
{
  static constexpr std::array<marker_T, 2> values = {LEFT, RIGHT};
};

template <>
struct AllSides<2U>
{
  static constexpr std::array<marker_T, 4> values = {BOTTOM, RIGHT, TOP, LEFT};
};

template <>
struct AllSides<3U>
{
  static constexpr std::array<marker_T, 6> values = {
      BOTTOM, RIGHT, TOP, LEFT, BACK, FRONT};
};

// =====================================================================
template <typename Mesh>
void buildFacets(Mesh & mesh, Bitmask<MeshFlags> flags = MeshFlags::NONE)
{
  using Elem_T = typename Mesh::Elem_T;
  using Facet_T = typename Mesh::Facet_T;
  // facets are identified by the ordered set of their point ids
  using FacetMapKey_T = std::set<id_T>;
  using FacetMap_T = std::map<FacetMapKey_T, Facet_T>;

  // store facets already present in the mesh
  FacetMap_T oldFacetMap;
  for (auto const & facet: mesh.facetList)
  {
    FacetMapKey_T facetIds;
    for (short_T p = 0; p < Mesh::Facet_T::numPts; ++p)
    {
      facetIds.insert(facet.pts[p]->id);
    }
    oldFacetMap.insert(std::pair{facetIds, facet});
  }

  uint facetCount = 0;
  uint iFacetCount = 0;
  uint oldFacetCount = 0;
  FacetMap_T facetMap;
  for (auto & e: mesh.elementList)
  {
    auto side = 0U;
    for (auto const & f: Elem_T::elemToFacet)
    {
      std::vector<Point *> facetPts(Mesh::Facet_T::numPts);
      FacetMapKey_T facetIds;
      for (short_T p = 0; p < Mesh::Facet_T::numPts; ++p)
      {
        facetPts[p] = e.pts[f[p]];
        facetIds.insert(e.pts[f[p]]->id);
      }
      Facet_T facet;
      // check if this facet was already present, create it otherwise
      if (auto it = oldFacetMap.find(facetIds); it != oldFacetMap.end())
      {
        facet = it->second;
        oldFacetCount++;
      }
      else
      {
        facet.pts = facetPts;
        facet.id = facetCount;
      }
      auto && [it, inserted] = facetMap.insert(std::pair(facetIds, facet));
      if (inserted)
      {
        // we are the first element to cross this facet
        it->second.facingElem[0].ptr = &e;
        it->second.facingElem[0].side = side;
        facetCount++;
      }
      else
      {
        // we are the second element crossing this facet, this is an internal facet
        it->second.facingElem[1].ptr = &e;
        it->second.facingElem[1].side = side;
        iFacetCount++;
      }
      side++;
    }
  }
  [[maybe_unused]] uint const oldBd = std::count_if(
      mesh.facetList.begin(),
      mesh.facetList.end(),
      [](Facet_T const & f) { return f.onBoundary(); });
  assert(oldFacetCount == oldBd + 2 * (mesh.facetList.size() - oldBd));

  uint bFacetSize = facetCount - iFacetCount;
  mesh.facetList.resize((flags & MeshFlags::INTERNAL_FACETS) ? facetCount : bFacetSize);
  mesh.elemToFacet.resize(
      mesh.elementList.size(), fillArray<Elem_T::numFacets>(dofIdNotSet));

  // store facets in mesh ordering boundary facets first
  iFacetCount = bFacetSize;
  uint bFacetCount = 0;
  // check
  // https://stackoverflow.com/questions/40673080/stdignore-with-structured-bindings
  for ([[maybe_unused]] auto const & [idSet, facet]: facetMap)
  {
    if (!facet.facingElem[1])
    {
      // this is a boundary facet
      mesh.facetList[bFacetCount] = facet;
      mesh.facetList[bFacetCount].id = bFacetCount;
      mesh.elemToFacet[facet.facingElem[0].ptr->id][facet.facingElem[0].side] =
          bFacetCount;
      bFacetCount++;
    }
    else
    {
      // this is an internal facet
      if (flags & MeshFlags::INTERNAL_FACETS)
      {
        mesh.facetList[iFacetCount] = facet;
        mesh.facetList[iFacetCount].id = iFacetCount;
        mesh.elemToFacet[facet.facingElem[0].ptr->id][facet.facingElem[0].side] =
            iFacetCount;
        mesh.elemToFacet[facet.facingElem[1].ptr->id][facet.facingElem[1].side] =
            iFacetCount;
      }
      iFacetCount++;
    }
  }
  assert(bFacetCount == bFacetSize);
  assert(iFacetCount == facetCount);
  // always set the BOUNDARY_FACETS flags
  mesh.flags |= MeshFlags::BOUNDARY_FACETS;
  mesh.numBdFacets = bFacetSize;
  // set INTERNAL_FACETS flag based on the input flags
  mesh.flags |= (flags & MeshFlags::INTERNAL_FACETS);
}

// =====================================================================
void buildLine(
    Mesh<Line> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    uint const numElems,
    Bitmask<MeshFlags> flags);

void buildSquare(
    Mesh<Triangle> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 2> const numElems,
    Bitmask<MeshFlags> flags);

void buildSquare(
    Mesh<Quad> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 2> const numElems,
    Bitmask<MeshFlags> flags);

void buildCircleMesh(
    Mesh<Quad> & mesh,
    Vec3 const & origin,
    double const & radius,
    std::array<uint, 3> const numElems);

void buildCube(
    Mesh<Tetrahedron> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 3> const numElems,
    Bitmask<MeshFlags> flags);

void buildCube(
    Mesh<Hexahedron> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 3> const numElems,
    Bitmask<MeshFlags> flags);

template <typename Elem>
void buildHyperCube(
    Mesh<Elem> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 3> const numElems,
    Bitmask<MeshFlags> flags = MeshFlags::NONE)
{
  if constexpr (std::is_same_v<Elem, Line>)
  {
    buildLine(mesh, origin, length, numElems[0], flags);
  }
  else if constexpr (std::is_same_v<Elem, Triangle> || std::is_same_v<Elem, Quad>)
  {
    buildSquare(mesh, origin, length, {{numElems[0], numElems[1]}}, flags);
  }
  else if constexpr (
      std::is_same_v<Elem, Tetrahedron> || std::is_same_v<Elem, Hexahedron>)
  {
    buildCube(mesh, origin, length, {{numElems[0], numElems[1], numElems[2]}}, flags);
  }
  else
  {
    // this should never happen
    fmt::print(stderr, "elememnt type {} not recognized\n", typeid(Elem{}).name());
    abort();
  }

  if (flags & MeshFlags::NORMALS)
  {
    buildNormals(mesh, flags);
  }
  if (flags & MeshFlags::FACET_PTRS)
  {
    addElemFacetList(mesh);
  }
}

template <typename Elem>
void buildHyperCube(Mesh<Elem> & mesh, ParameterDict const & config)
{
  auto flags = Bitmask{MeshFlags::NONE};
  if (config["flags"])
  {
    flags |= config["flags"].as<Bitmask<MeshFlags>>();
  }

  config.validate({"origin", "length", "n"});

  buildHyperCube(
      mesh,
      config["origin"].as<Vec3>(),
      config["length"].as<Vec3>(),
      config["n"].as<std::array<uint, 3>>(),
      flags);
}

template <typename Mesh>
void markFacetsCube(Mesh & mesh, Vec3 const & origin, Vec3 const & length)
{
  for (auto & f: mesh.facetList)
  {
    if (f.onBoundary())
    {
      auto const [min, max] = f.bbox();
      // BOTTOM -> ymin
      if (std::fabs(max[1] - origin[1]) < 1e-6 * length[1])
      {
        f.marker = side::BOTTOM;
      }
      // RIGHT -> xmax
      else if (std::fabs(min[0] - origin[0] - length[0]) < 1e-6 * length[0])
      {
        f.marker = side::RIGHT;
      }
      // TOP -> ymax
      else if (std::fabs(min[1] - origin[1] - length[1]) < 1e-6 * length[1])
      {
        f.marker = side::TOP;
      }
      // LEFT -> xmin
      else if (std::fabs(max[0] - origin[0]) < 1e-6 * length[0])
      {
        f.marker = side::LEFT;
      }
      // BACK -> zmin
      else if (std::fabs(max[2] - origin[2]) < 1e-6 * length[2])
      {
        f.marker = side::BACK;
      }
      // FRONT -> zmax
      else if (std::fabs(min[2] - origin[2] - length[2]) < 1e-6 * length[2])
      {
        f.marker = side::FRONT;
      }
    }
  }
}

void buildWedge(
    Mesh<Triangle> & mesh,
    Vec3 const & origin,
    Vec3 const & radius,
    Vec3 const & normal,
    uint const numLayers,
    double const angle = M_PI / 180.);

// =====================================================================
void extrude(
    Mesh<Triangle> const & mesh2d,
    Mesh<Tetrahedron> & mesh3d,
    uint const numLayers,
    Vec3 const & direction,
    double const distance);

// =====================================================================
void refTriangleMesh(Mesh<Triangle> & mesh);
void refQuadMesh(Mesh<Quad> & mesh);
void refTetrahedronMesh(Mesh<Tetrahedron> & mesh);
void refHexahedronMesh(Mesh<Hexahedron> & mesh);
void hexagonMesh(Mesh<Triangle> & mesh);
void hexagonSquare(Mesh<Triangle> & mesh, MeshFlags flags = MeshFlags::NONE);

template <typename Elem>
void referenceMesh(Mesh<Elem> & mesh)
{
  if constexpr (std::is_same_v<Elem, Line>)
  {
    buildHyperCube(
        mesh,
        {-1., 0., 0.},
        {2., 0., 0.},
        {1, 0, 0},
        MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS);
  }
  else if constexpr (std::is_same_v<Elem, Quad>)
  {
    refQuadMesh(mesh);
  }
  else if constexpr (std::is_same_v<Elem, Hexahedron>)
  {
    refHexahedronMesh(mesh);
  }
  else if constexpr (std::is_same_v<Elem, Triangle>)
  {
    refTriangleMesh(mesh);
  }
  else if constexpr (std::is_same_v<Elem, Tetrahedron>)
  {
    refTetrahedronMesh(mesh);
  }
  else
  {
    abort();
  }
}

// =====================================================================
template <typename Mesh>
void buildNormals(Mesh & mesh, Bitmask<MeshFlags> flags)
{
  if (!(mesh.flags & MeshFlags::NORMALS))
  {
    Bitmask<GeoElemFlags> facetFlags = GeoElemFlags::NONE;
    // CHECK_PLANAR makes sense only on quad facets
    if constexpr (std::is_same_v<typename Mesh::Elem_T, Hexahedron>)
    {
      facetFlags |= (flags & MeshFlags::CHECK_FACET_PLANAR) ? GeoElemFlags::CHECK_PLANAR
                                                            : GeoElemFlags::NONE;
    }
    for (auto & facet: mesh.facetList)
    {
      facet.buildNormal(facetFlags);

      // normals on boundary facets should all point outside
      if (facet.onBoundary())
      {
        if ((facet.midpoint() - facet.facingElem[0].ptr->midpoint())
                .dot(facet._normal) < 0.)
        {
          facet._normal *= -1.0;
        }
      }
      // all internal normals should point from facing elem 0 towards facing elem 1
      else if (
          (facet.facingElem[1].ptr->midpoint() - facet.facingElem[0].ptr->midpoint())
              .dot(facet._normal) < 0.)
      {
        facet._normal *= -1.0;
      }
    }
    mesh.flags |= MeshFlags::NORMALS;
  }
  else
  {
    fmt::print(stderr, "warning: mesh normals were already there\n");
  }
}

// =====================================================================
template <typename Mesh>
void addElemFacetList(Mesh & mesh)
{
  if (!(mesh.flags & MeshFlags::FACET_PTRS))
  {
    for (auto & elem: mesh.elementList)
    {
      elem.facets.resize(Mesh::Elem_T::numFacets);
    }
    for (auto & facet: mesh.facetList)
    {
      auto insideElem = facet.facingElem[0].ptr;
      auto const insidePos = facet.facingElem[0].side;
      insideElem->facets[insidePos] = &facet;
      auto outsideElem = facet.facingElem[1].ptr;
      if (outsideElem)
      {
        auto const outsidePos = facet.facingElem[1].side;
        outsideElem->facets[outsidePos] = &facet;
      }
    }
    mesh.flags |= MeshFlags::FACET_PTRS;
  }
  else
  {
    fmt::print(stderr, "warning: facet pointers were already there\n");
  }
}

// =====================================================================
enum class GMSHElemType : uint8_t
{
  NONE = 0,
  LINE = 1,
  TRIANGLE = 2,
  QUAD = 3,
  TET = 4,
  HEX = 5,
  QUADRATIC_LINE = 8,
  QUADRATIC_TRIANGLE = 9,
  QUADRATIC_QUAD = 10,
  QUADRATIC_TET = 11,
  BIQUADRATIC_HEX = 12,
  POINT = 15,
};

template <typename Elem, short_T Order = 1>
struct ElemToGmsh
{};

template <>
struct ElemToGmsh<PointElem, 1>
{
  static GMSHElemType constexpr value = GMSHElemType::POINT;
};

template <>
struct ElemToGmsh<PointElem, 2>
{
  static GMSHElemType constexpr value = GMSHElemType::POINT;
  static short_T constexpr connSize = 1U;
};

template <>
struct ElemToGmsh<Line, 1>
{
  static GMSHElemType constexpr value = GMSHElemType::LINE;
};

template <>
struct ElemToGmsh<Line, 2>
{
  static GMSHElemType constexpr value = GMSHElemType::QUADRATIC_LINE;
  static short_T constexpr connSize = 3U;
};

template <>
struct ElemToGmsh<Triangle, 1>
{
  static GMSHElemType constexpr value = GMSHElemType::TRIANGLE;
};

template <>
struct ElemToGmsh<Triangle, 2>
{
  static GMSHElemType constexpr value = GMSHElemType::QUADRATIC_TRIANGLE;
  static short_T constexpr connSize = 6U;
};

template <>
struct ElemToGmsh<Quad, 1>
{
  static GMSHElemType constexpr value = GMSHElemType::QUAD;
};

template <>
struct ElemToGmsh<Quad, 2>
{
  static GMSHElemType constexpr value = GMSHElemType::QUADRATIC_QUAD;
  static short_T constexpr connSize = 9U;
};

template <>
struct ElemToGmsh<Tetrahedron, 1>
{
  static GMSHElemType constexpr value = GMSHElemType::TET;
};

template <>
struct ElemToGmsh<Tetrahedron, 2>
{
  static GMSHElemType constexpr value = GMSHElemType::QUADRATIC_TET;
  static short_T constexpr connSize = 10U;
};

template <>
struct ElemToGmsh<Hexahedron, 1>
{
  static GMSHElemType constexpr value = GMSHElemType::HEX;
};

template <>
struct ElemToGmsh<Hexahedron, 2>
{
  static GMSHElemType constexpr value = GMSHElemType::BIQUADRATIC_HEX;
  static short_T constexpr connSize = 27U;
};

template <typename Elem>
struct ElemStub
{
  std::array<id_T, Elem::numPts> conn;
  id_T id = idNotSet;
  marker_T marker = markerNotSet;
};

template <typename Elem>
inline bool operator<(ElemStub<Elem> const & e1, ElemStub<Elem> const & e2)
{
  return e1.id < e2.id;
}

// =====================================================================
template <typename Elem>
void setFacetMarkers(
    Mesh<Elem> & mesh,
    std::set<ElemStub<typename Elem::Facet_T>> & facets,
    std::unordered_map<id_T, id_T> ptIdMap)
{
  using Facet_T = typename Elem::Facet_T;
  assert(
      std::count_if(
          mesh.facetList.begin(),
          mesh.facetList.end(),
          [](Facet_T const & f)
          { return f.onBoundary(); }) == static_cast<long int>(facets.size()));

  // use file facets to set boundary flags
  for (auto & meshFacet: mesh.facetList)
  {
    // iterator loop required to use std::set::erase()?
    for (auto it = facets.begin(); it != facets.end(); ++it)
    {
      std::vector<Point *> connPts;
      connPts.reserve(Facet_T::numPts);
      for (auto const c: it->conn)
      {
        connPts.push_back(&mesh.pointList[ptIdMap.at(c)]);
      };
      Facet_T tempFacet{connPts, it->id, it->marker};
      if (geoEqual(meshFacet, tempFacet))
      {
        meshFacet.marker = it->marker;
        facets.erase(it);
        break;
      }
    }
  }
  // check that all facets have been found in the mesh
  assert(facets.size() == 0);
}

// =====================================================================
template <typename Elem>
void readGMSH(
    Mesh<Elem> & mesh,
    std::string_view const filename,
    Bitmask<MeshFlags> flags = MeshFlags::NONE)
{
  auto in = std::ifstream(filename.data());
  if (!in.is_open())
  {
    fmt::print(stderr, "mesh file {} not found\n", filename);
    std::exit(PROXPDE_GMSH_ERROR);
  }

  // header
  std::string buf;
  in >> buf;
  if (buf != "$MeshFormat")
  {
    fmt::print(stderr, "file format not recognized\n");
    std::exit(PROXPDE_GMSH_ERROR);
  }

  int filetype, datasize; // filetype = 0 -> ascii
  in >> buf >> filetype >> datasize;
  if (buf != "2.2" || filetype != 0 || datasize != 8)
  {
    fmt::print(stderr, "file format not recognized\n");
    std::exit(PROXPDE_GMSH_ERROR);
  }

  in >> buf;
  if (buf != "$EndMeshFormat")
  {
    fmt::print(stderr, "file format not recognized\n");
    std::exit(PROXPDE_GMSH_ERROR);
  }

  std::vector<Point> points;
  std::vector<ElemStub<Elem>> elems;
  using Facet_T = typename Elem::Facet_T;
  using Ridge_T = typename Facet_T::Facet_T;
  std::set<ElemStub<Facet_T>> facets;
  std::set<marker_T> volumeMarkers;
  std::set<marker_T> facetMarkers;
  std::map<marker_T, std::string> physicalNames;

  auto ignoredElements = 0U;
  in >> buf;
  while (!in.eof())
  {
    if (buf == "$PhysicalNames")
    {
      int numNames;
      in >> numNames;
      for (short_T n = 0; n < numNames; ++n)
      {
        int dim;
        marker_T marker;
        std::string name;
        in >> dim >> marker >> name;
        name = name.substr(1, name.size() - 2);
        assert(dim == Elem::dim || dim == Elem::dim - 1);
        physicalNames[marker] = name;
      }
      in >> buf;
      if (buf != "$EndPhysicalNames")
      {
        fmt::print(stderr, "error reading physical names\n");
        std::exit(PROXPDE_GMSH_ERROR);
      }
    }
    else if (buf == "$Nodes")
    {
      uint numNodes;
      in >> numNodes;
      points.resize(numNodes);

      // format: node-number(one-based) x-coord y-coord z-coord
      for (uint n = 0; n < numNodes; n++)
      {
        uint id;
        double x, y, z;
        in >> id >> x >> y >> z;
        // currently only consecutive ids are supported
        assert(n == id - 1);
        points[n] = Point{Vec3{x, y, z}, n};
      }
      in >> buf;
      if (buf != "$EndNodes")
      {
        fmt::print(stderr, "error reading nodes\n");
        std::exit(PROXPDE_GMSH_ERROR);
      }
    }
    else if (buf == "$Elements")
    {
      uint numElements;
      in >> numElements;
      elems.reserve(numElements);
      // format: elm-number(one-based) elm-type number-of-tags < tag > ...
      // node-number-list
      uint eVol = 0, eBd = 0;
      for (uint e = 0; e < numElements; e++)
      {
        uint id;
        int tmp;
        uint numTags;
        in >> id >> tmp >> numTags;
        GMSHElemType elType = static_cast<GMSHElemType>(tmp);
        assert(numTags > 0);
        std::vector<uint> tags(numTags);
        for (uint t = 0; t < numTags; t++)
        {
          in >> tags[t];
        }

        // check element type (volume or boundary, linear or quadratic)
        if (ElemToGmsh<Elem, 1>::value == elType ||
            ElemToGmsh<Elem, 2>::value == elType)
        {
          volumeMarkers.insert(tags[0]);

          // read connectivity from file
          std::array<uint, Elem::numPts> conn;
          for (uint c = 0; c < Elem::numPts; c++)
          {
            in >> conn[c];
            // gmsh starts from 1
            conn[c]--;
            // register element is point
            points[conn[c]].neighboringElemSize++;
          }
          // read eventual additional connectivity from quadratic elements
          if (ElemToGmsh<Elem, 2>::value == elType)
          {
            for (uint c = Elem::numPts; c < ElemToGmsh<Elem, 2>::connSize; ++c)
            {
              in >> buf;
            }
          }

          // create mesh element
          // elems.emplace_back(connPts, eVol, tags[0]);
          elems.emplace_back(conn, eVol, tags[0]);
          eVol++;
        }
        else if (
            ElemToGmsh<Facet_T, 1>::value == elType ||
            ElemToGmsh<Facet_T, 2>::value == elType)
        {
          facetMarkers.insert(tags[0]);

          // read connectivity from file
          std::array<uint, Facet_T::numPts> conn;
          for (uint c = 0; c < Facet_T::numPts; c++)
          {
            in >> conn[c];
            // gmsh starts from 1
            conn[c]--;
          }
          // read possible additional connectivity from quadratic elements
          if (ElemToGmsh<Facet_T, 2>::value == elType)
          {
            for (uint c = Facet_T::numPts; c < ElemToGmsh<Facet_T, 2>::connSize; ++c)
            {
              in >> buf;
            }
          }

          // create mesh boundary element
          facets.insert(ElemStub<Facet_T>{conn, eBd, static_cast<marker_T>(tags[0])});
          eBd++;
        }
        else if (
            ElemToGmsh<Ridge_T, 1>::value == elType ||
            ElemToGmsh<Ridge_T, 2>::value == elType)
        {
          // ignore ridges
          for (uint c = 0; c < Ridge_T::numPts; c++)
          {
            in >> buf;
          }
          // read possible additional connectivity from quadratic ridges
          if (ElemToGmsh<Ridge_T, 2>::value == elType)
          {
            for (uint c = Ridge_T::numPts; c < ElemToGmsh<Ridge_T, 2>::connSize; ++c)
            {
              in >> buf;
            }
          }
          ignoredElements++;
        }
        else
        {
          fmt::print(stderr, "error reading elements\n");
          std::exit(PROXPDE_GMSH_ERROR);
        }
      }
      in >> buf;
      if (buf != "$EndElements")
      {
        fmt::print(stderr, "error detecting the end of the elements section\n");
        std::exit(PROXPDE_GMSH_ERROR);
      }
    }
    else
    {
      fmt::print(stderr, "file section not recognized: {}\n", buf);
      // discard the whole unrecognized section
      auto sectionEnd = "$End" + buf.substr(1, buf.length() - 1);
      while (buf != sectionEnd)
      {
        in >> buf;
      }
      // std::exit(PROXPDE_GMSH_ERROR);
    }
    // get next section
    in >> buf;
  }
  fmt::print("number of ignored elements: {}\n", ignoredElements);
  fmt::print("available volume markers: [");
  for (auto const m: volumeMarkers)
    fmt::print("{}, ", m);
  fmt::print("]\n");
  fmt::print("available facet markers: [");
  for (auto const m: facetMarkers)
    fmt::print("{}, ", m);
  fmt::print("]\n");
  fmt::print("physical names: [");
  for (auto const & n: physicalNames)
    fmt::print("{} -> {}, ", n.first, n.second);
  fmt::print("]\n");

  // TODO: store physical names in the mesh

  id_T ptCount = 0;
  std::unordered_map<id_T, id_T> ptIdMap;
  mesh.pointList.reserve(points.size());
  for (auto const & p: points)
  {
    if (p.neighboringElemSize > 0)
    {
      mesh.pointList.emplace_back(p.coord, ptCount, p.marker, p.neighboringElemSize);
      ptIdMap[p.id] = ptCount;
      ptCount++;
    }
  };

  mesh.elementList.reserve(elems.size());
  for (auto const & e: elems)
  {
    std::vector<Point *> connPts;
    connPts.reserve(Elem::numPts);
    for (auto const c: e.conn)
    {
      connPts.push_back(&mesh.pointList[ptIdMap.at(c)]);
    };
    mesh.elementList.emplace_back(connPts, e.id, e.marker);
  };

  mesh.buildConnectivity();

  // the file should contain all the boundary facets
  // TODO: this does not work if the mesh has internal boundaries
  buildFacets(mesh, flags);

  setFacetMarkers(mesh, facets, ptIdMap);

  if (flags & MeshFlags::NORMALS)
  {
    buildNormals(mesh, flags);
  }
  if (flags & MeshFlags::FACET_PTRS)
  {
    addElemFacetList(mesh);
  }

  fmt::print("mesh file {} successfully read\n", filename);
}

// =====================================================================
template <typename Elem>
void writeGMSH(Mesh<Elem> const & mesh, std::string const filename)
{
  std::ofstream out{filename};

  out << "$MeshFormat" << std::endl;
  out << "2.2 0 8" << std::endl;
  out << "$EndMeshFormat" << std::endl;

  out << "$Nodes" << std::endl;
  out << mesh.pointList.size() << std::endl;
  for (auto const & p: mesh.pointList)
  {
    out << p.id + 1 << " " << p.coord[0] << " " << p.coord[1] << " " << p.coord[2]
        << std::endl;
  }
  out << "$EndNodes" << std::endl;

  out << "$Elements" << std::endl;
  out << mesh.facetList.size() + mesh.elementList.size() << std::endl;
  for (auto const & f: mesh.facetList)
  {
    out << f.id + 1 << " "
        << static_cast<uint>(ElemToGmsh<typename Elem::Facet_T>::value) << " 2 "
        << f.marker << " " << f.marker << " ";
    for (uint i = 0; i < Elem::Facet_T::numPts; ++i)
    {
      out << f.pts[i]->id + 1 << " ";
    }
    out << std::endl;
  }
  for (auto const & e: mesh.elementList)
  {
    out << e.id + 1 + mesh.facetList.size() << " "
        << static_cast<uint>(ElemToGmsh<Elem>::value) << " 2 " << e.marker << " "
        << e.marker << " ";
    for (uint i = 0; i < Elem::numPts; ++i)
    {
      out << e.pts[i]->id + 1 << " ";
    }
    out << std::endl;
  }
  out << "$EndElements" << std::endl;
}

// =====================================================================
enum class MeshType : uint8_t
{
  STRUCTURED,
  GMSH,
};

static constexpr const char * to_string(MeshType const type)
{
  switch (type)
  {
  case MeshType::STRUCTURED:
    return "STRUCTURED";
  case MeshType::GMSH:
    return "GMSH";
  default:
    abort();
  }
}

static constexpr MeshType to_MeshType(std::string_view str)
{
  if (str == "STRUCTURED")
    return MeshType::STRUCTURED;
  else if (str == "GMSH")
    return MeshType::GMSH;
  fmt::print(stderr, "mesh type not recognized; {}\n", str);
  abort();
}

} // namespace proxpde

namespace fmt
{

template <>
struct formatter<proxpde::MeshType>: formatter<std::string_view>
{
  auto format(proxpde::MeshType t, format_context & ctx) const
  {
    return formatter<std::string_view>::format(proxpde::to_string(t), ctx);
  }
};

} // namespace fmt

namespace proxpde
{

// =====================================================================
template <typename Mesh>
void readMesh(Mesh & mesh, ParameterDict const & config)
{
  config.validate({"type"});
  auto const meshType = config["type"].as<MeshType>();

  switch (meshType)
  {
  case MeshType::STRUCTURED:
  {
    buildHyperCube(mesh, config);
    break;
  }
  case MeshType::GMSH:
  {
    Bitmask flags{MeshFlags::NONE};
    if (config["flags"])
    {
      flags |= config["flags"].as<Bitmask<MeshFlags>>();
    }
    readGMSH(mesh, config["filename"].as<std::string>(), flags);
    break;
  }
  default:
  {
    fmt::print(stderr, "mesh type not recognized: {}\n", meshType);
    abort();
  }
  }
}

// =====================================================================
template <typename Elem>
void buildFacetMesh(Mesh<typename Elem::Facet_T> & facetMesh, Mesh<Elem> const & mesh)
{
  facetMesh.pointList = mesh.pointList;
  facetMesh.elementList = mesh.facetList;
  facetMesh.buildConnectivity();
}

// =====================================================================
uint checkPlanarFacets(Mesh<Hexahedron> const & mesh);

// =====================================================================
template <typename Elem, typename Selector>
void extractBlock(
    proxpde::Mesh<Elem> const & mesh,
    proxpde::Mesh<Elem> & block,
    Selector const & selector)
{
  using namespace proxpde;
  std::vector<id_T> blockElemIds;
  blockElemIds.reserve(mesh.elementList.size());
  std::map<id_T, id_T> blockPtIds;
  id_T blockPtCounter = 0U;
  for (auto const & elem: mesh.elementList)
  {
    if (selector(elem))
    {
      blockElemIds.push_back(elem.id);
      for (auto const & pt: elem.pts)
      {
        [[maybe_unused]] auto const [ptr, inserted] =
            blockPtIds.insert({pt->id, blockPtCounter});
        if (inserted)
        {
          blockPtCounter++;
        }
      }
    }
  }

  fmt::print("number of elements in block: {}\n", blockElemIds.size());
  fmt::print("number of points in block: {}\n", blockPtCounter);

  block.pointList.resize(blockPtCounter);
  block.elementList.resize(blockElemIds.size());

  for (auto const [meshId, blockId]: blockPtIds)
  {
    block.pointList[blockId] = mesh.pointList[meshId];
    block.pointList[blockId].id = blockId;
  }

  id_T elemCounter = 0U;
  for (auto const meshId: blockElemIds)
  {
    Elem & blockElem = block.elementList[elemCounter];
    Elem const & meshElem = mesh.elementList[meshId];
    blockElem = meshElem;
    blockElem.id = elemCounter;
    for (uint p = 0U; p < Elem::numPts; p++)
    {
      id_T const blockId = blockPtIds[meshElem.pts[p]->id];
      blockElem.pts[p] = &block.pointList[blockId];
    }
    elemCounter++;
  }

  block.buildConnectivity();

  buildFacets(block, mesh.flags);
  // TODO: copy original markers

  if (mesh.flags & MeshFlags::NORMALS)
  {
    buildNormals(block, mesh.flags);
  }
  if (mesh.flags & MeshFlags::FACET_PTRS)
  {
    addElemFacetList(block);
  }

  block.flags = mesh.flags;

  // fmt::print("block:\n{}\n", block);
}

} // namespace proxpde

namespace YAML
{

template <>
struct convert<proxpde::Bitmask<proxpde::MeshFlags>>
{
  using T = proxpde::Bitmask<proxpde::MeshFlags>;
  static Node encode(T const & rhs)
  {
    // string array implementation
    Node node;
    for (auto const flag: proxpde::MeshFlags_ALL)
    {
      if (rhs & flag)
        node.push_back(proxpde::to_string(flag));
    }
    return node;

    // // ulong implementation
    // return Node{rhs.value};
  }

  static bool decode(Node const & node, T & rhs)
  {
    // string array implementation
    if (!node.IsSequence())
      return false;

    for (auto const & flag: node)
    {
      rhs |= proxpde::to_MeshFlags(flag.as<std::string>());
    }
    return true;

    // // ulong implementation
    // if (!node.IsScalar())
    //   return false;
    // rhs = node.as<T>();
    // return true;
  }
};

template <>
struct convert<proxpde::MeshType>
{
  using T = proxpde::MeshType;
  static Node encode(T const & rhs)
  {
    // string implementation
    return Node(proxpde::to_string(rhs));

    // // ulong implementation
    // return Node{static_cast<std::underlying_t<T>>(rhs)};
  }

  static bool decode(Node const & node, T & rhs)
  {
    // string implementation
    if (!node.IsScalar())
      return false;

    rhs = proxpde::to_MeshType(node.as<std::string>());
    return true;

    // // ulong implementation
    // if (!node.IsScalar())
    //   return false;
    // rhs = node.as<T>();
    // return true;
  }
};

} // namespace YAML
