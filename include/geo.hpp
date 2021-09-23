#pragma once

#include "def.hpp"

class Point
{
public:
  Point(Vec3 const c, id_T const i, marker_T const m = markerNotSet):
      coord(std::move(c)),
      id(i),
      marker(m)
  {}

  Point() = default;

  double operator[](short_T const i) const { return this->coord[i]; }

  Vec3 coord = Vec3::Zero();
  id_T id = dofIdNotSet;
  marker_T marker = markerNotSet;
};

inline std::ostream & operator<<(std::ostream & out, Point const & p)
{
  out << "(" << p[0] << "," << p[1] << "," << p[2] << "), id: " << p.id
      << ", m: " << p.marker;
  return out;
}

struct GeoElem;

struct FacingElem
{
  GeoElem * ptr;
  short_T side;
};

struct ChildElem
{
  GeoElem * ptr;
  short_T corner;
};

struct GeoElem
{
  using PointList_T = std::vector<Point *>;
  using FacetList_T = std::vector<GeoElem *>;
  GeoElem(std::initializer_list<Point *> const & pList, id_T const i, marker_T const m):
      pointList{pList},
      id{i},
      marker{m}
  {}

  GeoElem(PointList_T const pList, id_T const i, marker_T const m):
      pointList{std::move(pList)},
      id{i},
      marker{m}
  {}

  GeoElem() = default;

  virtual ~GeoElem();

  virtual Vec3 midpoint() const = 0;
  virtual Vec3 origin() const = 0;
  virtual double volume() const = 0;
  virtual void buildNormal() = 0;
  virtual Vec3 normal() const = 0;
  virtual double hMin() const = 0;
  virtual double hMax() const = 0;

  virtual std::tuple<Vec3, Vec3> bbox() const final
  {
    Vec3 min = Vec3::Constant(std::numeric_limits<double>::max());
    Vec3 max = Vec3::Constant(std::numeric_limits<double>::min());
    for (auto const & p: pointList)
    {
      for (uint c = 0; c < 3; ++c)
      {
        min[c] = std::min(min[c], p->coord[c]);
        max[c] = std::max(max[c], p->coord[c]);
      }
    }
    return std::tie(min, max);
  }

  // check if geoelem is on boundary
  bool onBoundary() const
  {
    // check that facingElems have been initialized
    assert(facingElem[0].ptr != nullptr);
    // the facet is on the boundary iff there is no outside element
    return facingElem[1].ptr == nullptr;
  }

  PointList_T pointList = {};
  FacetList_T facetList = {};
  id_T id = idNotSet;
  marker_T marker = markerNotSet;
  ChildElem parent = ChildElem{nullptr, shortNotSet};
  std::vector<ChildElem> children = {};
  // stores internal and external (element, facet side) pairs
  std::array<FacingElem, 2> facingElem = {
      FacingElem{nullptr, shortNotSet}, FacingElem{nullptr, shortNotSet}};
  Vec3 _normal;
};

inline bool operator<(GeoElem const & e1, GeoElem const & e2) { return e1.id < e2.id; }

inline std::ostream & operator<<(std::ostream & out, GeoElem const & e)
{
  out << "pts: ";
  for (auto & p: e.pointList)
  {
    out << p->id << " ";
  }
  out << "id: " << e.id << ", m: " << e.marker;
  auto const & [insideFacePtr, insideFaceLoc] = e.facingElem[0];
  if (insideFacePtr)
  {
    out << ", fe: (" << insideFacePtr->id << ", " << insideFaceLoc << ")";
    auto const [outsideFacePtr, outsideFaceLoc] = e.facingElem[1];
    if (outsideFacePtr)
    {
      out << ", (" << outsideFacePtr->id << ", " << outsideFaceLoc << ")";
    }
    else
    {
      out << ", on boundary";
    }
  }
  if (e.parent.ptr)
  {
    out << ", parent id: " << e.parent.ptr->id;
  }
  else
  {
    out << ", no parent";
  }
  out << ", children ids (" << e.children.size() << "): ";
  for (auto const ch: e.children)
  {
    out << ch.ptr->id << " ";
  }
  return out;
}

struct NullElem: public GeoElem
{
  static short_T constexpr dim = shortNotSet;
  static short_T constexpr numPts = 0U;
  static short_T constexpr numEdges = 0U;
  static short_T constexpr numFaces = 0U;
  static short_T constexpr numFacets = 0U;
};

struct PointElem: public GeoElem
{
  static short_T constexpr dim = 0;
  static short_T constexpr numPts = 1U;
  static short_T constexpr numEdges = 0U;
  static short_T constexpr numFaces = 0U;
  static short_T constexpr numFacets = 0U;

  PointElem(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  PointElem(PointList_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  PointElem() = default;

  ~PointElem() override = default;

  Vec3 midpoint() const final { return pointList[0]->coord; }

  Vec3 origin() const final { return pointList[0]->coord; }

  double volume() const final { return 1.; }

  void buildNormal() final { _normal = Vec3{1.0, 0.0, 0.0}; }

  Vec3 normal() const final { return Vec3{1.0, 0.0, 0.0}; }

  double hMin() const final { return 1.; }

  double hMax() const final { return 1.; }
};

class Line: public GeoElem
{
public:
  using Facet_T = PointElem;
  using Face_T = NullElem;
  using Edge_T = Line;
  static short_T constexpr dim = 1;
  static short_T constexpr numPts = 2U;
  static short_T constexpr numEdges = 1U;
  static short_T constexpr numFaces = 0U;
  static short_T constexpr numFacets = 2U;
  static std::array<std::array<short_T, 1>, 2> constexpr elemToFacet = {{{{0}}, {{1}}}};
  static std::array<std::array<short_T, 0>, 0> constexpr elemToFace = {{}};
  static std::array<std::array<short_T, 2>, 1> constexpr elemToEdge = {{{{0, 1}}}};
  static double constexpr refVolume = 2.;

  explicit Line(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Line(PointList_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Line() = default;

  ~Line() override = default;

  Vec3 midpoint() const final
  {
    return Vec3(0.5 * (pointList[1]->coord + pointList[0]->coord));
  }

  Vec3 origin() const final { return this->midpoint(); }

  double volume() const final
  {
    return (pointList[1]->coord - pointList[0]->coord).norm();
  }

  void buildNormal() final
  {
    Vec3 const length = pointList[1]->coord - pointList[0]->coord;
    _normal = Vec3{length[1], -length[0], 0.0};
    _normal.normalize();
  }

  Vec3 normal() const final
  {
    Vec3 const length = pointList[1]->coord - pointList[0]->coord;
    auto n = Vec3{length[1], -length[0], 0.0};
    n.normalize();
    return n;
  }

  double hMin() const final { return this->volume(); }

  double hMax() const final { return this->volume(); }
};

class Triangle: public GeoElem
{
public:
  using Facet_T = Line;
  using Face_T = Triangle;
  using Edge_T = Line;
  static short_T constexpr dim = 2;
  static short_T constexpr numPts = 3U;
  static short_T constexpr numEdges = 3U;
  static short_T constexpr numFaces = 1U;
  static short_T constexpr numFacets = 3U;
  static std::array<std::array<short_T, 2>, 3> constexpr elemToFacet = {
      {{{0, 1}}, {{1, 2}}, {{2, 0}}}};
  static std::array<std::array<short_T, 2>, 3> constexpr elemToEdge = elemToFacet;
  static std::array<std::array<short_T, 3>, 1> constexpr elemToFace = {{{{0, 1, 2}}}};
  static double constexpr refVolume = 0.5;

  Triangle(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Triangle(PointList_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Triangle() = default;

  ~Triangle() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{(pointList[0]->coord + pointList[1]->coord + pointList[2]->coord) / 3.};
  }

  Vec3 origin() const final { return pointList[0]->coord; }

  double volume() const final
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[0]->coord;
    return .5 * v1.cross(v2).norm();
  }

  void buildNormal() final
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[0]->coord;
    _normal = v1.cross(v2);
    _normal.normalize();
  }

  Vec3 normal() const final
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[0]->coord;
    auto n = v1.cross(v2);
    n.normalize();
    return n;
  }

  double perimeter() const
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[1]->coord;
    auto const v3 = pointList[0]->coord - pointList[2]->coord;

    return v1.norm() + v2.norm() + v3.norm();
  }

  double hMin() const final
  {
    // diameter of the inscribed circle
    return 4. * volume() / perimeter();
  }

  double hMax() const final
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[1]->coord;
    auto const v3 = pointList[0]->coord - pointList[2]->coord;

    // diameter of the circumscribed circle
    return v1.norm() * v2.norm() * v3.norm() / (v2.cross(v3).norm());
  }
};

class Quad: public GeoElem
{
public:
  using Facet_T = Line;
  using Face_T = Quad;
  using Edge_T = Line;
  static short_T constexpr dim = 2;
  static short_T constexpr numPts = 4U;
  static short_T constexpr numEdges = 4U;
  static short_T constexpr numFaces = 1U;
  static short_T constexpr numFacets = 4U;
  static std::array<std::array<short_T, 2>, 4> constexpr elemToFacet = {
      {{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}}};
  static std::array<std::array<short_T, 2>, 4> constexpr elemToEdge = elemToFacet;
  static std::array<std::array<short_T, 4>, 1> constexpr elemToFace = {
      {{{0, 1, 2, 3}}}};
  static double constexpr refVolume = 4.;

  Quad(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Quad(PointList_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Quad() = default;

  ~Quad() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{
        0.25 * (pointList[0]->coord + pointList[1]->coord + pointList[2]->coord +
                pointList[3]->coord)};
  }

  Vec3 origin() const final { return midpoint(); }

  double volume() const final
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[0]->coord;
    auto const v3 = pointList[3]->coord - pointList[0]->coord;
    return .5 * (v2.cross(v1).norm() + v2.cross(v3).norm());
  }

  void buildNormal() final
  {
    // TODO: we consider only planar quads
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[0]->coord;
    _normal = v1.cross(v2);
    _normal.normalize();
  }

  Vec3 normal() const final
  {
    // TODO: we consider only planar quads
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[0]->coord;
    auto n = v1.cross(v2);
    n.normalize();
    return n;
  }

  double perimeter() const
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[1]->coord;
    auto const v3 = pointList[3]->coord - pointList[2]->coord;
    auto const v4 = pointList[0]->coord - pointList[3]->coord;

    return v1.norm() + v2.norm() + v3.norm() + v4.norm();
  }

  double hMin() const final
  {
    return 4. * volume() / perimeter();
    // auto const d1 = pointList[2]->coord - pointList[0]->coord;
    // auto const d2 = pointList[3]->coord - pointList[1]->coord;
    // return std::min(d1.norm(), d2.norm());
  }

  double hMax() const final
  {
    // this holds only for cyclic quads
    auto const d1 = pointList[2]->coord - pointList[0]->coord;
    auto const d2 = pointList[3]->coord - pointList[1]->coord;
    auto const m1 = 0.5 * (pointList[2]->coord + pointList[0]->coord);
    auto const m2 = 0.5 * (pointList[3]->coord + pointList[1]->coord);
    auto const x = m1 - m2;
    return std::sqrt(0.5 * (d1.squaredNorm() + d2.squaredNorm() + x.squaredNorm()));
  }
};

class Tetrahedron: public GeoElem
{
public:
  using Facet_T = Triangle;
  using Face_T = Triangle;
  using Edge_T = Line;
  static short_T constexpr dim = 3;
  static short_T constexpr numPts = 4U;
  static short_T constexpr numEdges = 6U;
  static short_T constexpr numFaces = 4U;
  static short_T constexpr numFacets = 4U;
  static std::array<std::array<short_T, 3>, 4> constexpr elemToFacet = {
      {{{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}}, {{1, 2, 3}}}};
  static std::array<std::array<short_T, 2>, 6> constexpr elemToEdge = {
      {{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
  static std::array<std::array<short_T, 3>, 4> constexpr elemToFace = elemToFacet;
  static double constexpr refVolume = 1. / 6;

  Tetrahedron(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Tetrahedron(PointList_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Tetrahedron() = default;

  ~Tetrahedron() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{
        .25 * (pointList[0]->coord + pointList[1]->coord + pointList[2]->coord +
               pointList[3]->coord)};
  }

  Vec3 origin() const final { return pointList[0]->coord; }

  double volume() const final
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[0]->coord;
    auto const v3 = pointList[3]->coord - pointList[0]->coord;
    return std::fabs((v1.cross(v2)).dot(v3)) / 6.;
  }

  void buildNormal() final { std::abort(); }

  Vec3 normal() const final { std::abort(); }

  double hMin() const final
  {
    std::abort();
    return 1.;
  }

  double hMax() const final
  {
    std::abort();
    return 0;
  }
};

class Hexahedron: public GeoElem
{
public:
  using Facet_T = Quad;
  using Face_T = Quad;
  using Edge_T = Line;
  static short_T constexpr dim = 3;
  static short_T constexpr numPts = 8U;
  static short_T constexpr numEdges = 12U;
  static short_T constexpr numFaces = 6U;
  static short_T constexpr numFacets = 6U;
  static std::array<std::array<short_T, 4>, 6> constexpr elemToFacet = {
      {{{0, 3, 2, 1}},
       {{0, 1, 5, 4}},
       {{1, 2, 6, 5}},
       {{2, 3, 7, 6}},
       {{3, 0, 4, 7}},
       {{4, 5, 6, 7}}}};
  static std::array<std::array<short_T, 2>, 12> constexpr elemToEdge = {
      {{{0, 1}},
       {{1, 2}},
       {{2, 3}},
       {{3, 0}},
       {{0, 4}},
       {{1, 5}},
       {{2, 6}},
       {{3, 7}},
       {{4, 5}},
       {{5, 6}},
       {{6, 7}},
       {{7, 4}}}};
  static std::array<std::array<short_T, 4>, 6> constexpr elemToFace = elemToFacet;
  static double constexpr refVolume = 8.;

  Hexahedron(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Hexahedron(PointList_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Hexahedron() = default;

  ~Hexahedron() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{
        .125 * (pointList[0]->coord + pointList[1]->coord + pointList[2]->coord +
                pointList[3]->coord + pointList[4]->coord + pointList[5]->coord +
                pointList[6]->coord + pointList[7]->coord)};
  }

  Vec3 origin() const final { return this->midpoint(); }

  double volume() const final
  {
    // TODO: this is not correct for general hexahedrons, just for parallelepids
    // a better way is to split the hexahedron in tetrahedra and sum their volumes
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[3]->coord - pointList[0]->coord;
    auto const v3 = pointList[4]->coord - pointList[0]->coord;
    return std::fabs((v1.cross(v2)).dot(v3));
  }

  void buildNormal() final { std::abort(); }

  Vec3 normal() const final { std::abort(); }

  double hMin() const final
  {
    std::abort();
    return 1.;
  }

  double hMax() const final
  {
    std::abort();
    return 0;
  }
};

// this method checks only to see if the 2 elems are equivalent from the geometric pov,
// i.e. they have the same point ids (or a permutation of it)
template <typename Elem>
bool geoEqual(Elem const & e1, Elem const & e2)
{
  std::array<id_T, Elem::numPts> ids1, ids2;
  uint counter = 0;
  std::for_each(
      e1.pointList.begin(),
      e1.pointList.end(),
      [&ids1, &counter](Point const * p) { ids1[counter++] = p->id; });
  counter = 0;
  std::for_each(
      e2.pointList.begin(),
      e2.pointList.end(),
      [&ids2, &counter](Point const * p) { ids2[counter++] = p->id; });

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
    skew << 0., -axis[2], axis[1], axis[2], 0., -axis[0], -axis[1], axis[0], 0.;
    _m = std::cos(angle) * Mat_T::Identity() + std::sin(angle) * skew +
         (1. - std::cos(angle)) * (axis * axis.transpose());
    _mi = _m.inverse();
  }

  Vec3 apply(Vec3 const & v) const { return _m * v; }

  Vec3 applyInverse(Vec3 const & v) const { return _mi * v; }

private:
  Mat_T _m;
  Mat_T _mi;
};

// http://blackpawn.com/texts/pointinpoly/default.html
inline bool
sameSide2d(Vec3 const & p1, Vec3 const & p2, std::array<Vec3, 2> const & line)
{
  // the two points p1 and p2 are on the same side of the line through a and b
  // iff their cross product with the line itself point in the same direction
  Vec3 const diff = line[1] - line[0];
  return ((p1 - line[0]).cross(diff)).dot((p2 - line[0]).cross(diff)) >= 0.;
}

inline bool
sameSide3d(Vec3 const & p1, Vec3 const & p2, std::array<Vec3, 3> const & plane)
{
  // the two points p1 and p2 are on the same side of the plane through a, b and c
  // iff the dot products of of them with the normal to the plane has the same sign.
  Vec3 const normal = (plane[1] - plane[0]).cross(plane[2] - plane[0]);
  return (p1 - plane[0]).dot(normal) * (p2 - plane[0]).dot(normal) >= 0.;
}

template <typename Elem>
bool inside(Elem const & e, Vec3 const & pt)
{
  if constexpr (std::is_same_v<Elem, Line>)
  {
    Vec3 const & p0 = e.pointList[0]->coord;
    Vec3 const & p1 = e.pointList[1]->coord;
    if (pt[0] >= p0[0] && pt[0] <= p1[0])
      return true;
  }
  else if constexpr (std::is_same_v<Elem, Triangle>)
  {
    Vec3 const & p0 = e.pointList[0]->coord;
    Vec3 const & p1 = e.pointList[1]->coord;
    Vec3 const & p2 = e.pointList[2]->coord;
    // check that the point is on the same side of each edge wrt to the third point
    if (sameSide2d(pt, p0, {p1, p2}) && sameSide2d(pt, p1, {p2, p0}) &&
        sameSide2d(pt, p2, {p0, p1}))
      return true;
  }
  else if constexpr (std::is_same_v<Elem, Quad>)
  {
    // this approach requires 3 + 0.5 * 3 checks asymptotically
    // // split in two triangles and check for them
    // if (inside(Triangle{{e.pointList[0], e.pointList[1], e.pointList[2]}, idNotSet},
    // pt))
    //   return true;
    // else if (inside(Triangle{{e.pointList[2], e.pointList[3], e.pointList[0]},
    // idNotSet}, pt))
    //   return true;

    // this approach always requires 4 checks
    Vec3 const & p0 = e.pointList[0]->coord;
    Vec3 const & p1 = e.pointList[1]->coord;
    Vec3 const & p2 = e.pointList[2]->coord;
    Vec3 const & p3 = e.pointList[3]->coord;
    if (sameSide2d(pt, p0, {p2, p3}) && sameSide2d(pt, p2, {p0, p1}) &&
        sameSide2d(pt, p3, {p1, p2}) && sameSide2d(pt, p1, {p3, p0}))
      return true;
  }
  else if constexpr (std::is_same_v<Elem, Tetrahedron>)
  {
    Vec3 const & p0 = e.pointList[0]->coord;
    Vec3 const & p1 = e.pointList[1]->coord;
    Vec3 const & p2 = e.pointList[2]->coord;
    Vec3 const & p3 = e.pointList[3]->coord;
    // check that the point is on the same side of each face wrt to the forth point
    if (sameSide3d(pt, p3, {p0, p1, p2}) && sameSide3d(pt, p0, {p1, p2, p3}) &&
        sameSide3d(pt, p1, {p2, p3, p0}) && sameSide3d(pt, p2, {p3, p0, p1}))
      return true;
  }
  else if constexpr (std::is_same_v<Elem, Hexahedron>)
  {
    // the faces of an hexahedron may not lay on a plane.
    // in these cases all methods below are approximations that
    // relay on splitting each face in a number of planes

    // // split in 5 (* 8) tetrahedrons and check for them
    // asymptotic cost:
    //   (largest tet first)
    //   4 * (1 + 2/3 * 1/6 + 1/2 * 1/6 + 1/3 * 1/6 + 1/6 * 1/6) = 4 * 23 / 18 = 5.11111
    if (inside(
            Tetrahedron{
                {e.pointList[0], e.pointList[2], e.pointList[5], e.pointList[7]},
                idNotSet},
            pt))
      return true;
    else if (inside(
                 Tetrahedron{
                     {e.pointList[1], e.pointList[2], e.pointList[0], e.pointList[5]},
                     idNotSet},
                 pt))
      return true;
    else if (inside(
                 Tetrahedron{
                     {e.pointList[3], e.pointList[0], e.pointList[2], e.pointList[7]},
                     idNotSet},
                 pt))
      return true;
    else if (inside(
                 Tetrahedron{
                     {e.pointList[4], e.pointList[5], e.pointList[7], e.pointList[0]},
                     idNotSet},
                 pt))
      return true;
    else if (inside(
                 Tetrahedron{
                     {e.pointList[6], e.pointList[7], e.pointList[5], e.pointList[2]},
                     idNotSet},
                 pt))
      return true;

    // // minimum cost: 6 checks (1 per face)
    // // split each face in 2: 12 checks
    // // split each face in 8: 48 checks
    // Vec3 const & p0 = e.pointList[0]->coord;
    // Vec3 const & p1 = e.pointList[1]->coord;
    // Vec3 const & p2 = e.pointList[2]->coord;
    // Vec3 const & p3 = e.pointList[3]->coord;
    // Vec3 const & p4 = e.pointList[0]->coord;
    // Vec3 const & p5 = e.pointList[1]->coord;
    // Vec3 const & p6 = e.pointList[2]->coord;
    // Vec3 const & p7 = e.pointList[3]->coord;
    // if (sameSide3d(pt, p0, {p1, p2, p5}) &&
    //     sameSide3d(pt, p1, {p0, p3, p4}) &&
    //     sameSide3d(pt, p0, {p2, p3, p7}) &&
    //     sameSide3d(pt, p3, {p0, p1, p5}) &&
    //     sameSide3d(pt, p0, {p4, p5, p6}) &&
    //     sameSide3d(pt, p4, {p0, p1, p2}))
    //   return true;
  }
  else
  {
    abort();
  }
  return false;
}
