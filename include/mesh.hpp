#pragma once

#include "def.hpp"
#include "geo.hpp"

#include <set>
#include <unordered_set>
#include <map>
#include <fstream>

template <typename Elem>
class Mesh
{
public:
  using Elem_T = Elem;
  using Facet_T = typename Elem::Facet_T;
  using PointList_T = std::vector<Point>;
  using ElementList_T = std::vector<Elem>;
  using FacetList_T = std::vector<Facet_T>;
  using elemToPoint_T = std::vector<array<id_T,Elem::numPts>>;
  using elemToFacet_T = std::vector<array<id_T,Elem::numFacets>>;

  void buildConnectivity()
  {
    elemToPoint.reserve(elementList.size());
    for(auto& l: elementList)
    {
      array<id_T,Elem::numPts> elemConn;
      uint counter = 0;
      for(auto& p: l.pointList)
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
};

template <typename Elem>
std::ostream& operator<<(std::ostream& out, Mesh<Elem> const & mesh)
{
  out << "point list\n----------" << std::endl;
  for(auto &p: mesh.pointList)
    out << p << std::endl;
  out << "\nelement list\n------------" << std::endl;
  for(auto &e: mesh.elementList)
    out << e << std::endl;
  out << "\nfacet list\n------------" << std::endl;
  for(auto &f: mesh.facetList)
    out << f << std::endl;
  out << "\nelemToFacet map\n------------" << std::endl;
  for(auto &row: mesh.elemToFacet)
  {
    for(auto &clm: row)
    {
      out << clm << " ";
    }
    out << std::endl;
  }
  out << "\n------------" << std::endl;
  return out;
}

enum side
{
  BOTTOM = 1,
  RIGHT = 2,
  TOP = 3,
  LEFT = 4,
  CIRCLE = 101
};

template <typename Mesh>
void buildFacets(std::shared_ptr<Mesh> meshPtr, bool keep_internal = false)
{
  using Facet_T = typename Mesh::Facet_T;
  std::map<std::set<id_T>, Facet_T> facetMap;

  uint facetCount = 0;
  uint iFacetCount = 0;
  for(auto const & e: meshPtr->elementList)
  {
    uint side = 0;
    for(auto const & row: Mesh::Elem_T::elemToFacet)
    {
      std::vector<Point*> facetPts(Mesh::Facet_T::numPts);
      std::set<id_T> facetIds;
      uint i=0;
      for(auto const & clm: row)
      {
        facetPts[i] = e.pointList[clm];
        facetIds.insert(e.pointList[clm]->id);
        i++;
      }
      Facet_T facet(facetPts);
      auto inserted = facetMap.insert(std::make_pair(facetIds, facet));
      if(inserted.second)
      {
        // we are the first element to cross this facet
        inserted.first->second.facingElem[0].first = &e;
        inserted.first->second.facingElem[0].second = side;
        facetCount++;
      }
      else
      {
        // we are the second element crossing this facet, this is internal
        inserted.first->second.facingElem[1].first = &e;
        inserted.first->second.facingElem[1].second = side;
        iFacetCount++;
      }
      side++;
    }
  }

  uint bFacetTotal = facetCount - iFacetCount;
  if(keep_internal)
  {
    meshPtr->facetList.resize(facetCount);
  }
  else
  {
    meshPtr->facetList.resize(bFacetTotal);
  }
  meshPtr->elemToFacet.resize(meshPtr->elementList.size());
  for(uint i=0; i<meshPtr->elemToFacet.size(); ++i)
  {
    meshPtr->elemToFacet[i].fill(DOFidNotSet);
  }

  iFacetCount = bFacetTotal;
  uint bFacetCount = 0;
  for(auto const & facet: facetMap)
  {
    if(facet.second.facingElem[1].first == nullptr)
    {
      // this is a boundary facet
      meshPtr->facetList[bFacetCount] = facet.second;
      meshPtr->facetList[bFacetCount].id = bFacetCount;
      meshPtr->elemToFacet
        [facet.second.facingElem[0].first->id]
        [facet.second.facingElem[0].second] = bFacetCount;
      bFacetCount++;
    }
    else
    {
      // this is an internal facet
      if(keep_internal)
      {
        meshPtr->facetList[iFacetCount] = facet.second;
        meshPtr->facetList[iFacetCount].id = iFacetCount;
        meshPtr->elemToFacet
          [facet.second.facingElem[0].first->id]
          [facet.second.facingElem[0].second] = iFacetCount;
        meshPtr->elemToFacet
          [facet.second.facingElem[1].first->id]
          [facet.second.facingElem[1].second] = iFacetCount;
      }
      iFacetCount++;
    }
  }
  assert(bFacetCount == bFacetTotal);
  assert(iFacetCount == facetCount);
}

void buildMesh1D(std::shared_ptr<Mesh<Line>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 uint const numPts);

void buildMesh2D(std::shared_ptr<Mesh<Triangle>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 2> const numPts);

void buildMesh2D(std::shared_ptr<Mesh<Quad>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 2> const numPts);

void buildCircleMesh(std::shared_ptr<Mesh<Quad>> meshPtr,
                     Vec3 const& origin,
                     double const& radius,
                     array<uint, 3> const numPts);

template <class Elem>
struct MeshBuilder
{
  void build(std::shared_ptr<Mesh<Line>> meshPtr,
             Vec3 const& origin,
             Vec3 const& length,
             array<uint, 3> const numPts);
};

template <>
struct MeshBuilder<Line>
{
  void build(std::shared_ptr<Mesh<Line>> meshPtr,
             Vec3 const& origin,
             Vec3 const& length,
             array<uint, 3> const numPts)
  {
    buildMesh1D(meshPtr, origin, length, numPts[0]);
  }
};

template <>
struct MeshBuilder<Triangle>
{
  void build(std::shared_ptr<Mesh<Triangle>> meshPtr,
             Vec3 const& origin,
             Vec3 const& length,
             array<uint, 3> const numPts)
  {
    buildMesh2D(meshPtr, origin, length, {{numPts[0], numPts[1]}});
  }
};

template <>
struct MeshBuilder<Quad>
{
  void build(std::shared_ptr<Mesh<Quad>> meshPtr,
             Vec3 const& origin,
             Vec3 const& length,
             array<uint, 3> const numPts)
  {
    buildMesh2D(meshPtr, origin, length, {{numPts[0], numPts[1]}});
  }
};

template <typename Mesh>
void buildNormals(std::shared_ptr<Mesh> meshPtr)
{
  for (auto & facet: meshPtr->facetList)
  {
      if (facet.onBoundary())
      {
        facet.buildNormal();

        // check orientation wrt. midpoint of internal element
        if ((facet.midpoint() - facet.facingElem[0].first->midpoint()).dot(facet.normal) < 0.)
        {
          facet.normal = -1.0 * facet.normal;
        }

      }
  }
}

enum GMSHElemType
{
  GMSHNull = 0,
  GMSHLine = 1,
  GMSHTriangle = 2,
  GMSHQuad = 3,
//  GMSHTet = 4,
//  GMSHHexa = 5,
//  GMSHQuadraticLine = 8,
//  GMSHQuadraticTriangle = 9,
//  GMSHQuadraticQuad = 10,
//  GMSHQuadraticTet = 11,
//  GMSHQuadraticHexa = 12
};

template <typename Elem>
struct ElemToGmsh { static GMSHElemType constexpr type = GMSHNull; };

template <>
struct ElemToGmsh<Line> { static GMSHElemType constexpr type = GMSHLine; };
template <>
struct ElemToGmsh<Triangle> { static GMSHElemType constexpr type = GMSHTriangle; };
template <>
struct ElemToGmsh<Quad> { static GMSHElemType constexpr type = GMSHQuad; };

template <typename Elem>
void readGMSH(std::shared_ptr<Mesh<Elem>> meshPtr,
              std::string filename)
{
  auto in = std::ifstream(filename);
  if (!in.is_open())
  {
    std::cerr << "mesh file " << filename << " not found" << std::endl;
    std::exit(1);
  }

  // header
  std::string buf;
  in >> buf;
  if (buf != "$MeshFormat")
  {
    std::cerr << "file format not recognized" << std::endl;
    std::exit(1);
  }

  int filetype, datasize; // filetype = 0 -> ascii
  in >> buf >> filetype >> datasize;
  if (buf != "2.2" || filetype != 0 || datasize != 8)
  {
    std::cerr << "file format not recognized" << std::endl;
    std::exit(1);
  }

  in >> buf;
  if (buf != "$EndMeshFormat")
  {
    std::cerr << "file format not recognized" << std::endl;
    std::exit(1);
  }

  std::set<typename Elem::Facet_T> facets;

  in >> buf;
  while(!in.eof())
  {
    if (buf == "$PhysicalNames")
    {
      // TODO
      while(buf != "$EndPhysicalNames")
      {
        in >> buf;
      }
    }
    else if (buf == "$Nodes")
    {
      uint numNodes;
      in >> numNodes;
      meshPtr->pointList.reserve(numNodes);

      // format: node-number(one-based) x-coord y-coord z-coord
      for (uint n=0; n<numNodes; n++)
      {
        uint id;
        double x,y,z;
        in >> id >> x >> y >> z;
        // currently only sonsecutive ids are supported
        assert(n == id-1);
        meshPtr->pointList.emplace_back(Vec3{x, y, z}, n);
      }
      in >> buf;
      if (buf != "$EndNodes")
      {
        std::cerr << "error reading nodes" << std::endl;
        std::exit(1);
      }
    }
    else if (buf == "$Elements")
    {
      uint numElements;
      in >> numElements;
      meshPtr->elementList.reserve(numElements);
      // format: elm-number(one-based) elm-type number-of-tags < tag > ... node-number-list
      uint eVol = 0, eBd = 0;
      for (uint e=0; e<numElements; e++)
      {
        uint id, elType, numTags;
        in >> id >> elType >> numTags;
        std::vector<uint> tags(numTags);
        for (uint t=0; t<numTags; t++)
        {
          in >> tags[t];
        }

        // check if volume or boundary element
        if (ElemToGmsh<Elem>::type == elType)
        {
          // read connectivity from file
          array<uint, Elem::numPts> conn;
          for (uint c=0; c<Elem::numPts; c++)
          {
            in >> conn[c];
          }

          // get points pointers from connectivity
          std::vector<Point*> connPts;
          std::for_each(conn.begin(), conn.end(), [&connPts, &meshPtr](uint const c){
            connPts.push_back(&meshPtr->pointList[c-1]);
          });

          // create mesh element
          meshPtr->elementList.emplace_back(Elem(connPts, eVol));
          eVol++;
        }
        else if (ElemToGmsh<typename Elem::Facet_T>::type == elType)
        {
          // read connectivity from file
          array<uint, Elem::Facet_T::numPts> conn;
          for (uint c=0; c<Elem::Facet_T::numPts; c++)
          {
            in >> conn[c];
          }

          // get points pointers from connectivity
          std::vector<Point*> connPts;
          std::for_each(conn.begin(), conn.end(), [&connPts, &meshPtr](uint const c){
            connPts.push_back(&meshPtr->pointList[c-1]);
          });

          // create mesh boundary element
          facets.insert(typename Elem::Facet_T(connPts, eBd, tags[0]));
          eBd++;
        }
        else
        {
          std::cerr << "error reading elements" << std::endl;
          std::exit(1);
        }
      }
      in >> buf;
      if (buf != "$EndElements")
      {
        std::cerr << "error reading elements" << std::endl;
        std::exit(1);
      }
    }
    else
    {
      std::cerr << "file section not recognized" << std::endl;
      std::exit(1);
    }
    // get next section
    in >> buf;
  }
  meshPtr->buildConnectivity();
  buildFacets(meshPtr);

  // the file should contain all the boundary facets
  assert(meshPtr->facetList.size() == facets.size());

  // use file facets to set boundary flags
  for (auto & meshFacet: meshPtr->facetList)
  {
    for (auto it = facets.begin(); it != facets.end(); ++it)
    {
      if(geoEqual(meshFacet, *it))
      {
        meshFacet.marker = it->marker;
        facets.erase(it);
        break;
      }
    }
  }
  // check that all facets have been found in the mesh
  assert(facets.size() == 0);

  std::cout << "mesh file " << filename << " successfully read" << std::endl;
}
