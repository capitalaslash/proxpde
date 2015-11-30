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

struct GeoElem
{
  virtual Vec3 midpoint() const = 0;
  virtual double volume() const = 0;
};

class Line: public GeoElem
{
public:
  static uint const numPts = 2U;
  typedef Eigen::Vector2d LocalVec_T;
  typedef Eigen::Matrix2d LocalMat_T;
  typedef std::array<Point*,numPts> PointList_T;

  explicit Line(PointList_T const&& pl = {nullptr, nullptr}):
    pointList(pl)
    {}

  Vec3 midpoint() const
  {
    return Vec3(0.5*(pointList[1]->coord+pointList[0]->coord));
  }

  double volume() const
  {
    return (pointList[1]->coord-pointList[0]->coord).norm();
  }

  PointList_T pointList;
};

template <class Elem>
class Mesh
{
public:
  typedef std::vector<Point> PointList_T;
  typedef std::vector<Line> ElementList_T;
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

enum side
{
  LEFT,
  RIGHT
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
