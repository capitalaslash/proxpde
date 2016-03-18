#pragma once

#include <array>
#include <vector>
#include <memory>

#include "def.hpp"
#include "geo.hpp"

template <typename Elem>
class Mesh
{
public:
  typedef Elem Elem_T;
  typedef std::vector<Point> PointList_T;
  typedef std::vector<Elem> ElementList_T;
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
