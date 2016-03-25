#pragma once

#include "def.hpp"
#include "geo.hpp"

#include <set>
#include <map>

template <typename Elem>
class Mesh
{
public:
  typedef Elem Elem_T;
  typedef typename Elem::Facet_T Facet_T;
  typedef std::vector<Point> PointList_T;
  typedef std::vector<Elem> ElementList_T;
  typedef std::vector<Facet_T> FacetList_T;
  typedef std::array<id_T,Elem::numPts> ElementConn_T;
  typedef std::vector<std::array<id_T,Elem::numPts>> ConnList_T;

  void buildConnectivity()
  {
    _connList.reserve(elementList.size());
    for(auto& l: elementList)
    {
      ElementConn_T elemConn;
      uint counter = 0;
      for(auto& p: l.pointList)
      {
        elemConn[counter] = p->id;
        counter++;
      }
      _connList.push_back(elemConn);
    }
  }

  PointList_T pointList;
  ElementList_T elementList;
  FacetList_T facetList;
  ConnList_T _connList;
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
  out << "\n------------" << std::endl;
  return out;
}

enum side
{
  BOTTOM,
  RIGHT,
  TOP,
  LEFT
};

template <typename Mesh>
void buildFacets(std::shared_ptr<Mesh> meshPtr, bool keep_internal = false)
{
  typedef typename Mesh::Facet_T Facet_T;
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

  iFacetCount = bFacetTotal;
  uint bFacetCount = 0;
  for(auto const & facet: facetMap)
  {
    if(facet.second.facingElem[1].first == nullptr)
    {
      // this is a boundary facet
      meshPtr->facetList[bFacetCount] = facet.second;
      meshPtr->facetList[bFacetCount].id = bFacetCount;
      bFacetCount++;
    }
    else
    {
      // this is an internal facet
      if(keep_internal)
      {
        meshPtr->facetList[iFacetCount] = facet.second;
        meshPtr->facetList[iFacetCount].id = iFacetCount;
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
                 std::array<uint, 2> const numPts);

void buildMesh2D(std::shared_ptr<Mesh<Quad>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 std::array<uint, 2> const numPts);

template <class Elem>
struct MeshBuilder
{
  void build(std::shared_ptr<Mesh<Line>> meshPtr,
             Vec3 const& origin,
             Vec3 const& length,
             std::array<uint, 3> const numPts);
};

template <>
struct MeshBuilder<Line>
{
  void build(std::shared_ptr<Mesh<Line>> meshPtr,
             Vec3 const& origin,
             Vec3 const& length,
             std::array<uint, 3> const numPts)
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
             std::array<uint, 3> const numPts)
  {
    buildMesh2D(meshPtr, origin, length, {numPts[0], numPts[1]});
  }
};

template <>
struct MeshBuilder<Quad>
{
  void build(std::shared_ptr<Mesh<Quad>> meshPtr,
             Vec3 const& origin,
             Vec3 const& length,
             std::array<uint, 3> const numPts)
  {
    buildMesh2D(meshPtr, origin, length, {numPts[0], numPts[1]});
  }
};
