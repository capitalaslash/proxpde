#pragma once

#include <array>
#include <vector>
#include <memory>

#include "def.hpp"
#include "dof_object.hpp"

class Point: public DOFobject
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

struct GeoElem: public DOFobject
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

class Quad: public GeoElem
{
public:
  static uint const numPts = 4U;
  typedef std::array<Point*,numPts> PointList_T;

  explicit Quad(std::initializer_list<Point*> list, id_T id):
    GeoElem(list, id)
  {}

  Vec3 midpoint() const
  {
    return Vec3(0.25*(pointList[0]->coord
                     +pointList[1]->coord
                     +pointList[2]->coord
                     +pointList[3]->coord));
  }

  Vec3 origin() const
  {
    return midpoint();
  }

  double volume() const
  {
    return ((pointList[2]->coord-pointList[0]->coord).cross(
            (pointList[3]->coord-pointList[1]->coord))).norm();
  }
};
