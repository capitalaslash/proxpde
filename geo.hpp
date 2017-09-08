#pragma once

#include "def.hpp"

#include <algorithm>

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
  using PointList_T = std::vector<Point*>;
  using FacingElem_T = std::pair<GeoElem const*, uint>;

  explicit GeoElem(std::initializer_list<Point*> const & list = {nullptr},
                   id_T i = -1,
                   marker_T m = -1):
    pointList(list),
    id(i),
    marker(m),
    facingElem{{std::make_pair(nullptr, -1), std::make_pair(nullptr, -1)}}
  {}

  // TODO: use array instead of vector, the number of points is fixed by the elem type
  explicit GeoElem(std::vector<Point*> const & list,
                   id_T i = -1,
                   marker_T m = -1):
    pointList(list),
    id(i),
    marker(m),
    facingElem{{std::make_pair(nullptr, -1), std::make_pair(nullptr, -1)}}
  {}

  virtual Vec3 midpoint() const = 0;
  virtual Vec3 origin() const = 0;
  virtual double volume() const = 0;
  virtual void buildNormal() = 0;

  // check if geoelem is on boundary
  bool onBoundary()
  {
    return facingElem[1].first == nullptr && facingElem[0].first != nullptr;
  }

  PointList_T pointList;
  id_T id;
  marker_T marker;
  array<FacingElem_T, 2> facingElem;
  Vec3 normal;
};

inline bool operator< (GeoElem const & e1, GeoElem const & e2)
{
  return e1.id < e2.id;
}

inline std::ostream& operator<<(std::ostream& out, GeoElem const & e)
{
  out << "pts: ";
  for(auto & p: e.pointList)
  {
    out << p->id << " ";
  }
  out << "id: " << e.id << ", m: " << e.marker;
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
  static int const dim = -1;
  static uint const numPts = 0U;
  static uint const numEdges = 0U;
  static uint const numFaces = 0U;
  static uint const numFacets = 0U;
};

struct PointElem: public GeoElem
{
  static int const dim = 0;
  static uint const numPts = 1U;
  static uint const numEdges = 0U;
  static uint const numFaces = 0U;
  static uint const numFacets = 0U;

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

  virtual Vec3 midpoint() const final {return pointList[0]->coord;}
  virtual Vec3 origin() const {return pointList[0]->coord;}
  virtual double volume() const {return 1.;}
  virtual void buildNormal() final
  {
    normal = Vec3{0.0, 0.0, 0.0};
  }
};

class Line: public GeoElem
{
public:
  using Facet_T = PointElem;
  using Face_T = NullElem;
  using Edge_T = Line;
  static int const dim = 1;
  static uint const numPts = 2U;
  static uint const numEdges = 1U;
  static uint const numFaces = 0U;
  static uint const numFacets = 2U;
  static array<uint,4> constexpr dof_place{{0,0,0,1}};
  static array<array<id_T,1>,2> constexpr elemToFacet{{
    {{0}}, {{1}}
  }};
  static array<array<id_T,0>,0> constexpr elemToFace = {{}};
  static array<array<id_T,2>,1> constexpr elemToEdge{{
    {{0,1}}
  }};

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

  virtual void buildNormal() final
  {
    normal = Vec3{
        pointList[1]->coord[1]-pointList[0]->coord[1],
        pointList[0]->coord[0]-pointList[1]->coord[0],
        0.0};
    normal.normalize();
  }
};

class Triangle: public GeoElem
{
public:
  using Facet_T = Line;
  using Face_T = Triangle;
  using Edge_T = Line;
  static int const dim = 2;
  static uint const numPts = 3U;
  static uint const numEdges = 3U;
  static uint const numFaces = 1U;
  static uint const numFacets = 3U;
  static array<uint,4> constexpr dof_place{{0,0,0,1}};
  static array<array<id_T,2>,3> constexpr elemToFacet{
    {{{0,1}}, {{1,2}}, {{2,0}}}
  };
  static array<array<id_T,2>,3> constexpr elemToEdge = elemToFacet;
  static array<array<id_T,3>,1> constexpr elemToFace{
    {{{0,1,2}}}
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

  virtual void buildNormal() final
  {
    normal = (pointList[1]->coord-pointList[0]->coord).cross(
          pointList[2]->coord-pointList[0]->coord);
    normal.normalize();
  }
};

class Quad: public GeoElem
{
public:
  using Facet_T = Line;
  using Face_T = Quad;
  using Edge_T = Line;
  static int const dim = 2;
  static uint const numPts = 4U;
  static uint const numEdges = 4U;
  static uint const numFaces = 1U;
  static uint const numFacets = 4U;
  static array<array<id_T,2>,4> constexpr elemToFacet{{
      {{0,1}}, {{1,2}}, {{2,3}}, {{3,0}}
  }};
  static array<array<id_T,2>,4> constexpr elemToEdge = elemToFacet;
  static array<array<id_T,4>,1> constexpr elemToFace{{
    {{0,1,2,3}}
  }};

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

  virtual void buildNormal() final
  {
    normal = (pointList[1]->coord-pointList[0]->coord).cross(
          pointList[2]->coord-pointList[0]->coord);
    normal.normalize();
  }
};

// this method checks only to see if the 2 elems are equivalent from the geometric pov, i.e. they
// have the same point ids (or a permutation of it)
template <typename Elem>
bool geoEqual(Elem const & e1, Elem const & e2)
{
  std::array<id_T, Elem::numPts> ids1, ids2;
  uint counter= 0;
  std::for_each(e1.pointList.begin(), e1.pointList.end(), [&ids1, &counter](Point* const p)
  {
    ids1[counter++] = p->id;
  });
  counter = 0;
  std::for_each(e2.pointList.begin(), e2.pointList.end(), [&ids2, &counter](Point* const p)
  {
    ids2[counter++] = p->id;
  });

  return std::is_permutation(ids1.begin(), ids1.end(), ids2.begin());
}

// rotation matrix from axis and angle
// https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
class RotationMatrix
{
  using Mat_T = Eigen::Matrix3d;
public:
  RotationMatrix(Vec3 const & axis, double const angle)
  {
    Mat_T skew;
    skew <<       0., -axis[2],  axis[1],
             axis[2],       0., -axis[0],
            -axis[1],  axis[0],       0.;
    _m = std::cos(angle)*Mat_T::Identity() + std::sin(angle) * skew + (1.-std::cos(angle)) * (axis * axis.transpose());
    _mi = _m.inverse();
  }

  Vec3 apply(Vec3 const & v) const
  {
    return _m * v;
  }

  Vec3 applyInverse(Vec3 const & v) const
  {
    return _mi * v;
  }

private:
  Mat_T _m;
  Mat_T _mi;
};
