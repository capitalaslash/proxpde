#pragma once

#include "def.hpp"

#include "geo.hpp"

struct MeshFlags
{
  using T = std::bitset<4>;
  static constexpr T NONE = 0b0000;
  static constexpr T BOUNDARY_FACETS = 0b0001;
  static constexpr T INTERNAL_FACETS = 0b0010;
  static constexpr T NORMALS = 0b0100;
  static constexpr T FACET_PTRS = 0b1000;
};

template <typename Elem>
class Mesh
{
public:
  using Elem_T = Elem;
  using Facet_T = typename Elem::Facet_T;
  using PointList_T = std::vector<Point>;
  using ElementList_T = std::vector<Elem>;
  using FacetList_T = std::vector<Facet_T>;
  using elemToPoint_T = std::vector<array<id_T, Elem::numPts>>;
  using elemToFacet_T = std::vector<array<id_T, Elem::numFacets>>;

  void buildConnectivity()
  {
    elemToPoint.reserve(elementList.size());
    for (auto & l: elementList)
    {
      array<id_T, Elem::numPts> elemConn;
      uint counter = 0;
      for (auto & p: l.pointList)
      {
        elemConn[counter] = p->id;
        counter++;
      }
      elemToPoint.push_back(elemConn);
    }
  }

  PointList_T pointList;
  ElementList_T elementList;
  FacetList_T facetList;
  elemToPoint_T elemToPoint;
  elemToFacet_T elemToFacet;
  MeshFlags::T flags;
};

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

enum side : marker_T
{
  BOTTOM = 1,
  RIGHT = 2,
  TOP = 3,
  LEFT = 4,
  BACK = 5,
  FRONT = 6,
  CIRCLE = 101
};

template <typename Mesh>
void buildFacets(Mesh & mesh, MeshFlags::T flags = MeshFlags::NONE)
{
  using Facet_T = typename Mesh::Facet_T;
  std::map<std::set<id_T>, Facet_T> facetMap;

  uint facetCount = 0;
  uint iFacetCount = 0;
  for (auto & e: mesh.elementList)
  {
    uint side = 0;
    for (auto const & row: Mesh::Elem_T::elemToFacet)
    {
      std::vector<Point *> facetPts(Mesh::Facet_T::numPts);
      std::set<id_T> facetIds;
      uint i = 0;
      for (auto const & clm: row)
      {
        facetPts[i] = e.pointList[clm];
        facetIds.insert(e.pointList[clm]->id);
        i++;
      }
      Facet_T facet{facetPts, facetCount};
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

  uint bFacetSize = facetCount - iFacetCount;
  if ((flags & MeshFlags::INTERNAL_FACETS).any())
  {
    mesh.facetList.resize(facetCount);
  }
  else
  {
    mesh.facetList.resize(bFacetSize);
  }
  mesh.elemToFacet.resize(mesh.elementList.size());
  for (auto & row: mesh.elemToFacet)
  {
    row.fill(dofIdNotSet);
  }

  // store facets in mesh ordering boundary facets first
  iFacetCount = bFacetSize;
  uint bFacetCount = 0;
  // check
  // https://stackoverflow.com/questions/40673080/stdignore-with-structured-bindings
  for ([[maybe_unused]] auto const & [idSet, facet]: facetMap)
  {
    if (facet.facingElem[1].ptr == nullptr)
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
      if ((flags & MeshFlags::INTERNAL_FACETS).any())
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
  // set INTERNAL_FACETS flag based on the input flags
  mesh.flags |= flags & MeshFlags::INTERNAL_FACETS;
}

void buildLine(
    Mesh<Line> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    uint const numElems,
    MeshFlags::T flags);

void buildSquare(
    Mesh<Triangle> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    array<uint, 2> const numElems,
    MeshFlags::T flags);

void buildSquare(
    Mesh<Quad> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    array<uint, 2> const numElems,
    MeshFlags::T flags);

void buildCircleMesh(
    Mesh<Quad> & mesh,
    Vec3 const & origin,
    double const & radius,
    array<uint, 3> const numElems);

void buildCube(
    Mesh<Tetrahedron> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    array<uint, 3> const numElems,
    MeshFlags::T flags);

void buildCube(
    Mesh<Hexahedron> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    array<uint, 3> const numElems,
    MeshFlags::T flags);

template <typename Elem>
void buildHyperCube(
    Mesh<Elem> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    array<uint, 3> const numElems,
    MeshFlags::T flags = MeshFlags::NONE)
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
    std::cerr << "element type " << typeid(Elem{}).name() << " not recognized."
              << std::endl;
    abort();
  }

  if ((flags & MeshFlags::NORMALS).any())
  {
    buildNormals(mesh);
  }
  if ((flags & MeshFlags::FACET_PTRS).any())
  {
    addElemFacetList(mesh);
  }
}

void refTriangleMesh(Mesh<Triangle> & mesh);
void refQuadMesh(Mesh<Quad> & mesh);
void refTetrahedronMesh(Mesh<Tetrahedron> & mesh);
void refHexahedronMesh(Mesh<Hexahedron> & mesh);
void hexagonMesh(Mesh<Triangle> & mesh);
void hexagonSquare(Mesh<Triangle> & mesh, MeshFlags::T flags = MeshFlags::NONE);

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

template <typename Mesh>
void buildNormals(Mesh & mesh)
{
  if ((mesh.flags & MeshFlags::NORMALS).none())
  {
    for (auto & facet: mesh.facetList)
    {
      facet.buildNormal();

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
    std::cerr << "warning: mesh normals were already there" << std::endl;
  }
}

template <typename Mesh>
void addElemFacetList(Mesh & mesh)
{
  if ((mesh.flags & MeshFlags::FACET_PTRS).none())
  {
    for (auto & elem: mesh.elementList)
    {
      elem.facetList.resize(Mesh::Elem_T::numFacets);
    }
    for (auto & facet: mesh.facetList)
    {
      auto insideElem = facet.facingElem[0].ptr;
      auto const insidePos = facet.facingElem[0].side;
      insideElem->facetList[insidePos] = &facet;
      auto outsideElem = facet.facingElem[1].ptr;
      if (outsideElem)
      {
        auto const outsidePos = facet.facingElem[1].side;
        outsideElem->facetList[outsidePos] = &facet;
      }
    }
    mesh.flags |= MeshFlags::FACET_PTRS;
  }
  else
  {
    std::cerr << "warning: facet pointers were already there" << std::endl;
  }
}

enum GMSHElemType : int8_t
{
  GMSHNull = 0,
  GMSHLine = 1,
  GMSHTriangle = 2,
  GMSHQuad = 3,
  GMSHTet = 4,
  GMSHHexa = 5,
  // GMSHQuadraticLine = 8,
  // GMSHQuadraticTriangle = 9,
  // GMSHQuadraticQuad = 10,
  // GMSHQuadraticTet = 11,
  // GMSHQuadraticHexa = 12
};

template <typename Elem>
struct ElemToGmsh
{};

template <>
struct ElemToGmsh<Line>
{
  static GMSHElemType constexpr value = GMSHLine;
};

template <>
struct ElemToGmsh<Triangle>
{
  static GMSHElemType constexpr value = GMSHTriangle;
};

template <>
struct ElemToGmsh<Quad>
{
  static GMSHElemType constexpr value = GMSHQuad;
};

template <>
struct ElemToGmsh<Tetrahedron>
{
  static GMSHElemType constexpr value = GMSHTet;
};

template <>
struct ElemToGmsh<Hexahedron>
{
  static GMSHElemType constexpr value = GMSHHexa;
};

template <typename Elem>
void readGMSH(
    Mesh<Elem> & mesh,
    std::string_view const filename,
    MeshFlags::T flags = MeshFlags::NONE)
{
  auto in = std::ifstream(filename.data());
  if (!in.is_open())
  {
    std::cerr << "mesh file " << filename << " not found" << std::endl;
    std::exit(ERROR_GMSH);
  }

  // header
  std::string buf;
  in >> buf;
  if (buf != "$MeshFormat")
  {
    std::cerr << "file format not recognized" << std::endl;
    std::exit(ERROR_GMSH);
  }

  int filetype, datasize; // filetype = 0 -> ascii
  in >> buf >> filetype >> datasize;
  if (buf != "2.2" || filetype != 0 || datasize != 8)
  {
    std::cerr << "file format not recognized" << std::endl;
    std::exit(ERROR_GMSH);
  }

  in >> buf;
  if (buf != "$EndMeshFormat")
  {
    std::cerr << "file format not recognized" << std::endl;
    std::exit(ERROR_GMSH);
  }

  std::set<typename Elem::Facet_T> facets;

  in >> buf;
  while (!in.eof())
  {
    if (buf == "$PhysicalNames")
    {
      // TODO: stop discarding physical names
      while (buf != "$EndPhysicalNames")
      {
        in >> buf;
      }
    }
    else if (buf == "$Nodes")
    {
      uint numNodes;
      in >> numNodes;
      mesh.pointList.reserve(numNodes);

      // format: node-number(one-based) x-coord y-coord z-coord
      for (uint n = 0; n < numNodes; n++)
      {
        uint id;
        double x, y, z;
        in >> id >> x >> y >> z;
        // currently only sonsecutive ids are supported
        assert(n == id - 1);
        mesh.pointList.emplace_back(Vec3{x, y, z}, n);
      }
      in >> buf;
      if (buf != "$EndNodes")
      {
        std::cerr << "error reading nodes" << std::endl;
        std::exit(ERROR_GMSH);
      }
    }
    else if (buf == "$Elements")
    {
      uint numElements;
      in >> numElements;
      mesh.elementList.reserve(numElements);
      // format: elm-number(one-based) elm-type number-of-tags < tag > ...
      // node-number-list
      uint eVol = 0, eBd = 0;
      for (uint e = 0; e < numElements; e++)
      {
        uint id, elType, numTags;
        in >> id >> elType >> numTags;
        assert(numTags > 0);
        std::vector<uint> tags(numTags);
        for (uint t = 0; t < numTags; t++)
        {
          in >> tags[t];
        }

        // check if volume or boundary element
        if (ElemToGmsh<Elem>::value == elType)
        {
          // read connectivity from file
          array<uint, Elem::numPts> conn;
          for (uint c = 0; c < Elem::numPts; c++)
          {
            in >> conn[c];
          }

          // get points pointers from connectivity
          // TODO: use array or pre-fix size
          std::vector<Point *> connPts;
          std::for_each(
              conn.begin(),
              conn.end(),
              [&connPts, &mesh](uint const c)
              { connPts.push_back(&mesh.pointList[c - 1]); });

          // create mesh element
          mesh.elementList.emplace_back(Elem(connPts, eVol));
          eVol++;
        }
        else if (ElemToGmsh<typename Elem::Facet_T>::value == elType)
        {
          // read connectivity from file
          array<uint, Elem::Facet_T::numPts> conn;
          for (uint c = 0; c < Elem::Facet_T::numPts; c++)
          {
            in >> conn[c];
          }

          // get points pointers from connectivity
          // TODO: use array or pre-fix size
          std::vector<Point *> connPts;
          std::for_each(
              conn.begin(),
              conn.end(),
              [&connPts, &mesh](uint const c)
              { connPts.push_back(&mesh.pointList[c - 1]); });

          // create mesh boundary element
          facets.insert(typename Elem::Facet_T(connPts, eBd, tags[0]));
          eBd++;
        }
        else
        {
          std::cerr << "error reading elements" << std::endl;
          std::exit(ERROR_GMSH);
        }
      }
      in >> buf;
      if (buf != "$EndElements")
      {
        std::cerr << "error reading elements" << std::endl;
        std::exit(ERROR_GMSH);
      }
    }
    else
    {
      std::cerr << "file section not recognized: " << buf << std::endl;
      // discard the whole unrecognized section
      auto sectionEnd = "$End" + buf.substr(1, buf.length() - 1);
      while (buf != sectionEnd)
      {
        in >> buf;
      }
      // std::exit(ERROR_GMSH);
    }
    // get next section
    in >> buf;
  }
  mesh.buildConnectivity();
  buildFacets(mesh, flags);

  // the file should contain all the boundary facets
  // TODO: this does not work if the mesh has internal boundaries
  assert(
      std::count_if(
          mesh.facetList.begin(),
          mesh.facetList.end(),
          [](typename Elem::Facet_T const & f)
          { return f.onBoundary(); }) == static_cast<long int>(facets.size()));

  // use file facets to set boundary flags
  for (auto & meshFacet: mesh.facetList)
  {
    for (auto it = facets.begin(); it != facets.end(); ++it)
    {
      if (geoEqual(meshFacet, *it))
      {
        meshFacet.marker = it->marker;
        facets.erase(it);
        break;
      }
    }
  }
  // check that all facets have been found in the mesh
  assert(facets.size() == 0);

  if ((flags & MeshFlags::NORMALS).any())
  {
    buildNormals(mesh);
  }
  if ((flags & MeshFlags::FACET_PTRS).any())
  {
    addElemFacetList(mesh);
  }
  std::cout << "mesh file " << filename << " successfully read" << std::endl;
}

template <typename Mesh>
void readMesh(
    Mesh & mesh, YAML::Node const & config, MeshFlags::T flags = MeshFlags::NONE)
{
  auto const mesh_type = config["mesh_type"].as<std::string>();
  if (mesh_type == "structured")
  {
    array<uint, 3> numElems = {
        config["nx"].as<uint>(), config["ny"].as<uint>(), config["nz"].as<uint>()};
    // TODO: auto const origin = config["origin"].as<Vec3>();
    Vec3 const origin{0., 0., 0.};
    Vec3 const length{1., 1., 1.};
    buildHyperCube(mesh, origin, length, numElems, flags);
  }
  else if (mesh_type == "msh")
  {
    readGMSH(mesh, config["mesh_file"].as<std::string>(), flags);
  }
}

template <typename Elem>
void buildFacetMesh(Mesh<typename Elem::Facet_T> & facetMesh, Mesh<Elem> const & mesh)
{
  facetMesh.pointList = mesh.pointList;
  facetMesh.elementList = mesh.facetList;
  facetMesh.buildConnectivity();
}
