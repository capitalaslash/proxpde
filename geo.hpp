#pragma once

#include "def.hpp"

class Point
{
public:
  explicit Point(Vec3 const & c = Vec3::Zero(),
                 id_T const i = -1,
                 marker_T const m = -1):
    coord(c),
    id(i),
    dof_id(DOFidNotSet),
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
  DOFid_T dof_id;
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
  typedef std::pair<GeoElem const*, uint> FacingElem_T;

  explicit GeoElem(std::initializer_list<Point*> const & list = {nullptr},
                   id_T i = -1,
                   marker_T m = -1):
    pointList(list),
    id(i),
    marker(m),
    facingElem{std::make_pair(nullptr, -1), std::make_pair(nullptr, -1)}
  {}

  explicit GeoElem(std::vector<Point*> const & list,
                   id_T i = -1,
                   marker_T m = -1):
    pointList(list),
    id(i),
    marker(m),
    facingElem{std::make_pair(nullptr, -1), std::make_pair(nullptr, -1)}
  {}

  virtual Vec3 midpoint() const = 0;
  virtual Vec3 origin() const = 0;
  virtual double volume() const = 0;

  PointList_T pointList;
  id_T id;
  marker_T marker;
  std::array<FacingElem_T, 2> facingElem;
};

inline std::ostream& operator<<(std::ostream& out, GeoElem const & e)
{
  out << "(";
  for(auto & p: e.pointList)
  {
    out << p->id << ", ";
  }
  out << "\b\b), " << "id: " << e.id << ", m: " << e.marker;
  if(e.facingElem[0].first)
  {
    out << ", fe: (" << e.facingElem[0].first->id << ", " << e.facingElem[0].second << ")";
    if(e.facingElem[1].first)
    {
      out << ", " << e.facingElem[1].first->id << ", " << e.facingElem[1].second << ")";
    }
    else
    {
      out << ", on boundary";
    }
  }
  return out;
}

struct NullElem: public GeoElem
{
  static uint const numPts = 0U;
  static uint const numEdges = 0U;
  static uint const numFaces = 0U;
};

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
  typedef PointElem Facet_T;
  typedef NullElem Face_T;
  typedef Line Edge_T;
  static uint const numPts = 2U;
  static uint const numEdges = 1U;
  static uint const numFaces = 0U;
  static std::array<uint,4> constexpr dof_place{0,0,0,1};
  static std::array<std::array<id_T,1>,2> constexpr elemToFacet{
    {{0}, {1}}
  };
  static std::array<std::array<id_T,0>,0> constexpr elemToFace = {{}};
  static std::array<std::array<id_T,2>,1> constexpr elemToEdge{
    {{0,1}}
  };

  explicit Line(std::initializer_list<Point*> const & list = {nullptr},
                id_T i = -1,
                marker_T m = -1):
    GeoElem(list, i, m)
  {}

  explicit Line(std::vector<Point*> const & list,
                id_T i = -1,
                marker_T m = -1):
    GeoElem(list, i, m)
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
  typedef Line Facet_T;
  typedef Triangle Face_T;
  typedef Line Edge_T;
  static uint const numPts = 3U;
  static uint const numEdges = 3U;
  static uint const numFaces = 1U;
  static std::array<uint,4> constexpr dof_place{0,0,0,1};
  static std::array<std::array<id_T,2>,3> constexpr elemToFacet{
    {{0,1}, {1,2}, {2,0}}
  };
  static std::array<std::array<id_T,2>,3> constexpr elemToEdge = elemToFacet;
  static std::array<std::array<id_T,3>,1> constexpr elemToFace{
    {{0,1,2}}
  };

  explicit Triangle(std::initializer_list<Point*> const & list = {nullptr},
                    id_T i = -1,
                    marker_T m = -1):
        GeoElem(list, i, m)
  {}

  explicit Triangle(std::vector<Point*> const & list,
                id_T i = -1,
                marker_T m = -1):
    GeoElem(list, i, m)
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
  typedef Line Facet_T;
  typedef Quad Face_T;
  typedef Line Edge_T;
  static uint const numPts = 4U;
  static uint const numEdges = 4U;
  static uint const numFaces = 1U;
  static std::array<std::array<id_T,2>,4> constexpr elemToFacet{
    {{0,1}, {1,2}, {2,3}, {3,0}}
  };
  static std::array<std::array<id_T,2>,4> constexpr elemToEdge = elemToFacet;
  static std::array<std::array<id_T,4>,1> constexpr elemToFace{
    {{0,1,2,3}}
  };

  explicit Quad(std::initializer_list<Point*> list = {nullptr},
                id_T i = -1,
                marker_T m = -1):
    GeoElem(list, i, m)
  {}

  explicit Quad(std::vector<Point*> const & list,
                id_T i = -1,
                marker_T m = -1):
    GeoElem(list, i, m)
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
