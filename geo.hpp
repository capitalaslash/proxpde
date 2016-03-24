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

  explicit GeoElem(std::initializer_list<Point*> const & list = {nullptr},
                   id_T i = -1,
                   marker_T m = -1):
    pointList(list),
    id(i),
    marker(m),
    facingElem{nullptr, nullptr}
  {}

  explicit GeoElem(std::vector<Point*> const & list,
                   id_T i = -1,
                   marker_T m = -1):
    pointList(list),
    id(i),
    marker(m),
    facingElem{nullptr, nullptr}
  {}

  virtual Vec3 midpoint() const = 0;
  virtual Vec3 origin() const = 0;
  virtual double volume() const = 0;

  PointList_T pointList;
  id_T id;
  marker_T marker;
  std::array<GeoElem const*, 2> facingElem;
};

inline std::ostream& operator<<(std::ostream& out, GeoElem const & e)
{
  out << "(";
  for(auto & p: e.pointList)
  {
    out << p->id << ", ";
  }
  out << "\b\b), " << "id: " << e.id << ", m: " << e.marker;
  if(e.facingElem[0])
  {
    out << ", fe: " << e.facingElem[0]->id;
    if(e.facingElem[1])
    {
      out << ", " << e.facingElem[1]->id;
    }
    else
    {
      out << ", on boundary";
    }
  }
  return out;
}

struct PointElem: public GeoElem
{
  static uint const numPts = 1U;
  static uint const numEdges = 0U;
  static uint const numFaces = 0U;

  explicit PointElem(std::initializer_list<Point*> list = {nullptr},
                 id_T i = -1,
                 marker_T m = -1):
    GeoElem(list, i, m)
  {}

  explicit PointElem(std::vector<Point*> const & list,
                id_T i = -1,
                marker_T m = -1):
    GeoElem(list, i, m)
  {}

  virtual Vec3 midpoint() const {return pointList[0]->coord;}
  virtual Vec3 origin() const {return pointList[0]->coord;}
  virtual double volume() const {return 1.;}
};

class Line: public GeoElem
{
public:
  static uint const numPts = 2U;
  typedef Eigen::Vector2d LocalVec_T;
  typedef Eigen::Matrix2d LocalMat_T;
  typedef std::vector<Point*> PointList_T;

  explicit Line(std::initializer_list<Point*> list, id_T id = -1):
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

  explicit Triangle(std::initializer_list<Point*> list, id_T id = -1):
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

  explicit Quad(std::initializer_list<Point*> list, id_T id = -1):
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
