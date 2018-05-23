#pragma once

#include "def.hpp"

#include <algorithm>

class Point
{
public:
  explicit Point(Vec3 const & c = Vec3::Zero(),
                 id_T const i = DOFidNotSet,
                 marker_T const m = MarkerNotSet):
    coord(c),
    id(i),
    marker(m)
  {}
  explicit Point(double const x,
                 double const y,
                 double const z,
                 id_T const i = DOFidNotSet,
                 marker_T const m = MarkerNotSet):
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
                   id_T const i = DOFidNotSet,
                   marker_T const m = MarkerNotSet):
    pointList(list),
    id(i),
    marker(m),
    facingElem{{std::pair(nullptr, -1), std::pair(nullptr, -1)}}
  {}

  // TODO: use array instead of vector, the number of points is fixed by the elem type
  explicit GeoElem(std::vector<Point*> const & list,
                   id_T const i = DOFidNotSet,
                   marker_T const m = MarkerNotSet):
    pointList(list),
    id(i),
    marker(m),
    facingElem{{std::pair(nullptr, -1), std::pair(nullptr, -1)}}
  {}

  virtual ~GeoElem() {}

  virtual Vec3 midpoint() const = 0;
  virtual Vec3 origin() const = 0;
  virtual double volume() const = 0;
  virtual void buildNormal() = 0;

  virtual std::tuple<Vec3, Vec3> bbox() const final
  {
    Vec3 min = Vec3::Constant(std::numeric_limits<double>::max());
    Vec3 max = Vec3::Constant(std::numeric_limits<double>::min());
    for (auto const & p: pointList)
    {
      for(uint c=0; c<3; ++c)
      {
        min[c] = std::min(min[c], p->coord[c]);
        max[c] = std::max(max[c], p->coord[c]);
      }
    }
    return std::tie(min, max);
  }

  // check if geoelem is on boundary
  bool onBoundary()
  {
    // the facet is on the boundary iff there is an inside element and no outside element
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
  auto const & [insideFacePtr, insideFaceLoc] = e.facingElem[0];
  if (insideFacePtr)
  {
    out << ", fe: (" << insideFacePtr->id << ", " << insideFaceLoc << ")";
    auto const [outsideFacePtr, outsideFaceLoc] = e.facingElem[1];
    if(outsideFacePtr)
    {
      out << ", " << outsideFacePtr->id << ", " << outsideFaceLoc << ")";
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
  static int constexpr dim = -1;
  static uint constexpr numPts = 0U;
  static uint constexpr numEdges = 0U;
  static uint constexpr numFaces = 0U;
  static uint constexpr numFacets = 0U;
};

struct PointElem: public GeoElem
{
  static int constexpr dim = 0;
  static uint constexpr numPts = 1U;
  static uint constexpr numEdges = 0U;
  static uint constexpr numFaces = 0U;
  static uint constexpr numFacets = 0U;

  explicit PointElem(std::initializer_list<Point*> const & list = {nullptr},
                 id_T const i = DOFidNotSet,
                 marker_T const m = MarkerNotSet):
    GeoElem(list, i, m)
  {}

  explicit PointElem(std::vector<Point*> const & list,
                id_T const i = DOFidNotSet,
                marker_T const m = MarkerNotSet):
    GeoElem(list, i, m)
  {}

  virtual ~PointElem() {}

  virtual Vec3 midpoint() const final {return pointList[0]->coord;}
  virtual Vec3 origin() const final {return pointList[0]->coord;}
  virtual double volume() const final {return 1.;}
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
  static int constexpr dim = 1;
  static uint constexpr numPts = 2U;
  static uint constexpr numEdges = 1U;
  static uint constexpr numFaces = 0U;
  static uint constexpr numFacets = 2U;
  static array<array<id_T,1>,2> constexpr elemToFacet = {{
    {{0}}, {{1}}
  }};
  static array<array<id_T,0>,0> constexpr elemToFace = {{}};
  static array<array<id_T,2>,1> constexpr elemToEdge = {{
    {{0,1}}
  }};

  explicit Line(std::initializer_list<Point*> const & list = {nullptr},
                id_T const i = DOFidNotSet,
                marker_T const m = MarkerNotSet):
    GeoElem(list, i, m)
  {}

  explicit Line(std::vector<Point*> const & list,
                id_T const i = DOFidNotSet,
                marker_T const m = MarkerNotSet):
    GeoElem(list, i, m)
  {}

  virtual ~Line() {}

  Vec3 midpoint() const final
  {
    return Vec3(0.5*(pointList[1]->coord+pointList[0]->coord));
  }

  Vec3 origin() const final
  {
    return this->midpoint();
  }

  double volume() const final
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
  static int constexpr dim = 2;
  static uint constexpr numPts = 3U;
  static uint constexpr numEdges = 3U;
  static uint constexpr numFaces = 1U;
  static uint constexpr numFacets = 3U;
  static array<array<id_T,2>,3> constexpr elemToFacet = {{
    {{0,1}}, {{1,2}}, {{2,0}}
  }};
  static array<array<id_T,2>,3> constexpr elemToEdge = elemToFacet;
  static array<array<id_T,3>,1> constexpr elemToFace = {{
    {{0,1,2}}
  }};

  explicit Triangle(std::initializer_list<Point*> const & list = {nullptr},
                    id_T const i = DOFidNotSet,
                    marker_T const m = MarkerNotSet):
        GeoElem(list, i, m)
  {}

  explicit Triangle(std::vector<Point*> const & list,
                id_T const i = DOFidNotSet,
                marker_T const m = MarkerNotSet):
    GeoElem(list, i, m)
  {}

  virtual ~Triangle() {}

  Vec3 midpoint() const final
  {
    return Vec3{
      (pointList[0]->coord + pointList[1]->coord + pointList[2]->coord)/3.
    };
  }

  Vec3 origin() const final
  {
    return pointList[0]->coord;
  }

  double volume() const final
  {
    auto const v1 = pointList[1]->coord-pointList[0]->coord;
    auto const v2 = pointList[2]->coord-pointList[0]->coord;
    return .5 * v1.cross(v2).norm();
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
  static int constexpr dim = 2;
  static uint constexpr numPts = 4U;
  static uint constexpr numEdges = 4U;
  static uint constexpr numFaces = 1U;
  static uint constexpr numFacets = 4U;
  static array<array<id_T,2>,4> constexpr elemToFacet = {{
      {{0,1}}, {{1,2}}, {{2,3}}, {{3,0}}
  }};
  static array<array<id_T,2>,4> constexpr elemToEdge = elemToFacet;
  static array<array<id_T,4>,1> constexpr elemToFace = {{
    {{0,1,2,3}}
  }};

  explicit Quad(std::initializer_list<Point*> const & list = {nullptr},
                id_T const i = DOFidNotSet,
                marker_T const m = MarkerNotSet):
    GeoElem(list, i, m)
  {}

  explicit Quad(std::vector<Point*> const & list,
                id_T const i = DOFidNotSet,
                marker_T const m = MarkerNotSet):
    GeoElem(list, i, m)
  {}

  virtual ~Quad() {}

  Vec3 midpoint() const final
  {
    return Vec3{
      0.25*(pointList[0]->coord + pointList[1]->coord +
            pointList[2]->coord + pointList[3]->coord)
    };
  }

  Vec3 origin() const final
  {
    return midpoint();
  }

  double volume() const final
  {
    auto const v1 = pointList[1]->coord-pointList[0]->coord;
    auto const v2 = pointList[2]->coord-pointList[0]->coord;
    auto const v3 = pointList[3]->coord-pointList[0]->coord;
    return .5*(v2.cross(v1).norm() + v2.cross(v3).norm());
  }

  virtual void buildNormal() final
  {
    // TODO: we consider only planar quads
    auto const v1 = pointList[1]->coord-pointList[0]->coord;
    auto const v2 = pointList[2]->coord-pointList[0]->coord;
    normal = v1.cross(v2);
    normal.normalize();
  }
};

// this method checks only to see if the 2 elems are equivalent from the geometric pov, i.e. they
// have the same point ids (or a permutation of it)
template <typename Elem>
bool geoEqual(Elem const & e1, Elem const & e2)
{
  array<id_T, Elem::numPts> ids1, ids2;
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
