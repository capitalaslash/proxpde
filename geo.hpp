#pragma once

#include <array>
#include <vector>
#include <memory>

#include "def.hpp"

class Point
{
public:
  explicit Point(Vec3 const & c = Vec3::Zero(),
                 id_T const i = -1,
                 marker_T const m = -1):
    coord(c),
    id(i),
    marker(m)
  {}
  explicit Point(double const x,
                 double const y,
                 double const z,
                 id_T const i = -1,
                 marker_T const m = -1):
    Point(Vec3(x, y, z), i, m)
  {}

  double operator()(uint const i) const
  {
    return this->coord(i);
  }

  Vec3 coord;
  id_T id;
  marker_T marker;
};

inline std::ostream& operator<<(std::ostream& out, Point const & p)
{
  out << "(" << p.coord[0] << "," << p.coord[1] << "," << p.coord[2] << "), id: "
      << p.id << ", m: " << p.marker;
  return out;
}

struct GeoElem
{
  typedef std::vector<Point*> PointList_T;

  explicit GeoElem(std::initializer_list<Point*> list, id_T i):
    pointList(list),
    id(i)
  {}

  virtual Vec3 midpoint() const = 0;
  virtual Vec3 origin() const = 0;
  virtual double volume() const = 0;

  PointList_T pointList;
  id_T id;
};

inline std::ostream& operator<<(std::ostream& out, GeoElem const & e)
{
  out << "(";
  for(auto & p: e.pointList)
  {
    out << p->id << ", ";
  }
  out << "\b\b), " << "id: " << e.id;
  return out;
}

class Line: public GeoElem
{
public:
  static uint const numPts = 2U;
  typedef Eigen::Vector2d LocalVec_T;
  typedef Eigen::Matrix2d LocalMat_T;
  typedef std::vector<Point*> PointList_T;

  explicit Line(std::initializer_list<Point*> list, id_T id):
    GeoElem(list, id)
  {}

  Vec3 midpoint() const
  {
    return Vec3(0.5*(pointList[1]->coord+pointList[0]->coord));
  }

  Vec3 origin() const
  {
    return this->midpoint();
  }

  double volume() const
  {
    return (pointList[1]->coord-pointList[0]->coord).norm();
  }
};

class Triangle: public GeoElem
{
public:
  static uint const numPts = 3U;
  typedef std::array<Point*,numPts> PointList_T;

  explicit Triangle(std::initializer_list<Point*> list, id_T id):
    GeoElem(list, id)
  {}

  Vec3 midpoint() const
  {
    return Vec3((pointList[0]->coord
                +pointList[1]->coord
                +pointList[2]->coord)/3.);
  }

  Vec3 origin() const
  {
    return pointList[0]->coord;
  }

  double volume() const
  {
    return ((pointList[1]->coord-pointList[0]->coord).cross(
            (pointList[2]->coord-pointList[0]->coord))).norm();
  }
};

template <typename Elem>
class Mesh
{
public:
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
