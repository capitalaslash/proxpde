#pragma once

#include "def.hpp"

struct Point
{
  double operator[](short_T const i) const { return this->coord[i]; }

  Vec3 coord = Vec3::Zero();
  id_T id = idNotSet;
  marker_T marker = markerNotSet;
  short_T neighboringElemSize = 0;
};

std::ostream & operator<<(std::ostream & out, Point const & p);

struct GeoElem;

struct FacingElem
{
  GeoElem * ptr;
  short_T side;

  explicit operator bool() const { return ptr; }
};

inline bool operator==(FacingElem const & e1, FacingElem const & e2)
{
  return e1.ptr == e2.ptr && e1.side == e2.side;
}

struct ChildElem
{
  GeoElem * ptr;
  short_T corner;

  explicit operator bool() const { return ptr; }
};

inline bool operator==(ChildElem const & e1, ChildElem const & e2)
{
  return e1.ptr == e2.ptr && e1.corner == e2.corner;
}

struct GeoElem
{
  using Pts_T = std::vector<Point *>;
  using Facets_T = std::vector<GeoElem *>;
  GeoElem(std::initializer_list<Point *> const & pList, id_T const i, marker_T const m):
      pts{pList},
      id{i},
      marker{m}
  {}

  GeoElem(Pts_T const pList, id_T const i, marker_T const m):
      pts{std::move(pList)},
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
    Vec3 max = -Vec3::Constant(std::numeric_limits<double>::max());
    for (auto const & p: pts)
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

  Pts_T pts = {};
  Facets_T facets = {};
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

// this method checks only to see if the 2 elems are equivalent from the geometric pov,
// i.e. they have the same point ids (or a permutation of it)
bool geoEqual(GeoElem const & e1, GeoElem const & e2);

bool operator==(GeoElem const & e1, GeoElem const & e2);

std::ostream & operator<<(std::ostream & out, GeoElem const & e);

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
  using Facet_T = NullElem;
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

  PointElem(Pts_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  PointElem() = default;

  ~PointElem() override = default;

  Vec3 midpoint() const final { return pts[0]->coord; }

  Vec3 origin() const final { return pts[0]->coord; }

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
  static array2d<short_T, numFacets, Facet_T::numPts> constexpr elemToFacet = {
      {{{0}}, {{1}}}};
  static array2d<short_T, 0, 0> constexpr elemToFace = {{}};
  static array2d<short_T, numEdges, Edge_T::numPts> constexpr elemToEdge = {{{{0, 1}}}};
  static short_T constexpr numChildren = 2U;
  static array2d<short_T, numChildren, numPts> constexpr elemToChild = {
      {{{0, 2}}, {{2, 1}}}};
  static std::array<FMat<numPts, numPts>, numChildren> const embeddingMatrix;
  static double constexpr refVolume = 2.;

  explicit Line(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Line(Pts_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Line() = default;

  ~Line() override = default;

  Vec3 midpoint() const final { return Vec3(0.5 * (pts[1]->coord + pts[0]->coord)); }

  Vec3 origin() const final { return this->midpoint(); }

  double volume() const final { return (pts[1]->coord - pts[0]->coord).norm(); }

  void buildNormal() final
  {
    Vec3 const length = pts[1]->coord - pts[0]->coord;
    _normal = Vec3{length[1], -length[0], 0.0};
    _normal.normalize();
  }

  Vec3 normal() const final
  {
    Vec3 const length = pts[1]->coord - pts[0]->coord;
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
  static array2d<short_T, numFacets, Facet_T::numPts> constexpr elemToFacet = {
      {{{0, 1}}, {{1, 2}}, {{2, 0}}}};
  static array2d<short_T, numEdges, Edge_T::numPts> constexpr elemToEdge = elemToFacet;
  static array2d<short_T, numFaces, Face_T::numPts> constexpr elemToFace = {
      {{{0, 1, 2}}}};
  static short_T constexpr numChildren = 4U;
  static array2d<short_T, numChildren, numPts> constexpr elemToChild = {
      {{{0, 3, 5}}, {{1, 4, 3}}, {{2, 5, 4}}, {{3, 4, 5}}}};
  static std::array<FMat<numPts, numPts>, numChildren> const embeddingMatrix;
  static array3d<
      short_T,
      numFacets,
      Facet_T::numChildren,
      Facet_T::numPts> constexpr elemToFacetChild = {{
      {{{{0, 3}}, {{3, 1}}}},
      {{{{1, 4}}, {{4, 2}}}},
      {{{{2, 5}}, {{5, 0}}}},
  }};
  static array3d<
      short_T,
      numFacets,
      Facet_T::numChildren,
      2> constexpr elemToFacetChildFacing = {{
      {{{{0, 0}}, {{1, 2}}}},
      {{{{1, 0}}, {{2, 2}}}},
      {{{{2, 0}}, {{0, 2}}}},
  }};
  static double constexpr refVolume = 0.5;

  Triangle(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Triangle(Pts_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Triangle() = default;

  ~Triangle() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{(pts[0]->coord + pts[1]->coord + pts[2]->coord) / 3.};
  }

  Vec3 origin() const final { return pts[0]->coord; }

  double volume() const final
  {
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[0]->coord;
    return .5 * v1.cross(v2).norm();
  }

  void buildNormal() final
  {
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[0]->coord;
    _normal = v1.cross(v2);
    _normal.normalize();
  }

  Vec3 normal() const final
  {
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[0]->coord;
    auto n = v1.cross(v2);
    n.normalize();
    return n;
  }

  double perimeter() const
  {
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[1]->coord;
    auto const v3 = pts[0]->coord - pts[2]->coord;

    return v1.norm() + v2.norm() + v3.norm();
  }

  double hMin() const final
  {
    // diameter of the inscribed circle
    return 4. * volume() / perimeter();
  }

  double hMax() const final
  {
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[1]->coord;
    auto const v3 = pts[0]->coord - pts[2]->coord;

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
  static array2d<short_T, numFacets, Facet_T::numPts> constexpr elemToFacet = {
      {{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}}};
  static array2d<short_T, numEdges, Edge_T::numPts> constexpr elemToEdge = elemToFacet;
  static array2d<short_T, numFaces, Face_T::numPts> constexpr elemToFace = {
      {{{0, 1, 2, 3}}}};
  static short_T constexpr numChildren = 4U;
  static array2d<short_T, numPts, numChildren> constexpr elemToChild = {
      {{{0, 4, 8, 7}}, {{4, 1, 5, 8}}, {{8, 5, 2, 6}}, {{7, 8, 6, 3}}}};
  static std::array<FMat<numPts, numPts>, numChildren> const embeddingMatrix;
  static array3d<
      short_T,
      numFacets,
      Facet_T::numChildren,
      Facet_T::numPts> constexpr elemToFacetChild = {{
      {{{{0, 4}}, {{4, 1}}}},
      {{{{1, 5}}, {{5, 2}}}},
      {{{{2, 6}}, {{6, 3}}}},
      {{{{3, 7}}, {{7, 0}}}},
  }};
  static array3d<
      short_T,
      numFacets,
      Facet_T::numChildren,
      Facet_T::numPts> constexpr elemToFacetChildFacing = {{
      {{{{0, 0}}, {{1, 0}}}},
      {{{{1, 1}}, {{2, 1}}}},
      {{{{2, 2}}, {{3, 2}}}},
      {{{{3, 3}}, {{0, 3}}}},
  }};
  static double constexpr refVolume = 4.;

  Quad(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Quad(Pts_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Quad() = default;

  ~Quad() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{0.25 * (pts[0]->coord + pts[1]->coord + pts[2]->coord + pts[3]->coord)};
  }

  Vec3 origin() const final { return midpoint(); }

  double volume() const final
  {
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[0]->coord;
    auto const v3 = pts[3]->coord - pts[0]->coord;
    return .5 * (v2.cross(v1).norm() + v2.cross(v3).norm());
  }

  void buildNormal() final
  {
    // TODO: we consider only planar quads
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[0]->coord;
    _normal = v1.cross(v2);
    _normal.normalize();
  }

  Vec3 normal() const final
  {
    // TODO: we consider only planar quads
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[0]->coord;
    auto n = v1.cross(v2);
    n.normalize();
    return n;
  }

  double perimeter() const
  {
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[1]->coord;
    auto const v3 = pts[3]->coord - pts[2]->coord;
    auto const v4 = pts[0]->coord - pts[3]->coord;

    return v1.norm() + v2.norm() + v3.norm() + v4.norm();
  }

  double hMin() const final
  {
    return 4. * volume() / perimeter();
    // auto const d1 = pts[2]->coord - pts[0]->coord;
    // auto const d2 = pts[3]->coord - pts[1]->coord;
    // return std::min(d1.norm(), d2.norm());
  }

  double hMax() const final
  {
    // this holds only for cyclic quads
    auto const d1 = pts[2]->coord - pts[0]->coord;
    auto const d2 = pts[3]->coord - pts[1]->coord;
    auto const m1 = 0.5 * (pts[2]->coord + pts[0]->coord);
    auto const m2 = 0.5 * (pts[3]->coord + pts[1]->coord);
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
  static array2d<short_T, numFacets, Facet_T::numPts> constexpr elemToFacet = {
      {{{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}}, {{1, 2, 3}}}};
  static array2d<short_T, numEdges, Edge_T::numPts> constexpr elemToEdge = {
      {{{0, 1}}, {{1, 2}}, {{2, 0}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};
  static array2d<short_T, numFaces, Face_T::numPts> constexpr elemToFace = elemToFacet;
  static double constexpr refVolume = 1. / 6;

  Tetrahedron(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Tetrahedron(Pts_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Tetrahedron() = default;

  ~Tetrahedron() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{.25 * (pts[0]->coord + pts[1]->coord + pts[2]->coord + pts[3]->coord)};
  }

  Vec3 origin() const final { return pts[0]->coord; }

  double volume() const final
  {
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[2]->coord - pts[0]->coord;
    auto const v3 = pts[3]->coord - pts[0]->coord;
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
  static array2d<short_T, numFacets, Facet_T::numPts> constexpr elemToFacet = {
      {{{0, 3, 2, 1}},
       {{0, 1, 5, 4}},
       {{1, 2, 6, 5}},
       {{2, 3, 7, 6}},
       {{3, 0, 4, 7}},
       {{4, 5, 6, 7}}}};
  static array2d<short_T, numEdges, Edge_T::numPts> constexpr elemToEdge = {
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
  static array2d<short_T, numFaces, Face_T::numPts> constexpr elemToFace = elemToFacet;
  static double constexpr refVolume = 8.;

  Hexahedron(
      std::initializer_list<Point *> const & pList,
      id_T const i,
      marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Hexahedron(Pts_T const & pList, id_T const i, marker_T const m = markerNotSet):
      GeoElem{pList, i, m}
  {}

  Hexahedron() = default;

  ~Hexahedron() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{
        .125 * (pts[0]->coord + pts[1]->coord + pts[2]->coord + pts[3]->coord +
                pts[4]->coord + pts[5]->coord + pts[6]->coord + pts[7]->coord)};
  }

  Vec3 origin() const final { return this->midpoint(); }

  double volume() const final
  {
    // TODO: this is not correct for general hexahedrons, just for parallelepids
    // a better way is to split the hexahedron in tetrahedra and sum their volumes
    auto const v1 = pts[1]->coord - pts[0]->coord;
    auto const v2 = pts[3]->coord - pts[0]->coord;
    auto const v3 = pts[4]->coord - pts[0]->coord;
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
    Vec3 const & p0 = e.pts[0]->coord;
    Vec3 const & p1 = e.pts[1]->coord;
    if (pt[0] >= p0[0] && pt[0] <= p1[0])
      return true;
  }
  else if constexpr (std::is_same_v<Elem, Triangle>)
  {
    Vec3 const & p0 = e.pts[0]->coord;
    Vec3 const & p1 = e.pts[1]->coord;
    Vec3 const & p2 = e.pts[2]->coord;
    // check that the point is on the same side of each edge wrt to the third point
    if (sameSide2d(pt, p0, {p1, p2}) && sameSide2d(pt, p1, {p2, p0}) &&
        sameSide2d(pt, p2, {p0, p1}))
      return true;
  }
  else if constexpr (std::is_same_v<Elem, Quad>)
  {
    // this approach requires 3 + 0.5 * 3 checks asymptotically
    // // split in two triangles and check for them
    // if (inside(Triangle{{e.pts[0], e.pts[1], e.pts[2]}, idNotSet},
    // pt))
    //   return true;
    // else if (inside(Triangle{{e.pts[2], e.pts[3], e.pts[0]},
    // idNotSet}, pt))
    //   return true;

    // this approach always requires 4 checks
    Vec3 const & p0 = e.pts[0]->coord;
    Vec3 const & p1 = e.pts[1]->coord;
    Vec3 const & p2 = e.pts[2]->coord;
    Vec3 const & p3 = e.pts[3]->coord;
    if (sameSide2d(pt, p0, {p2, p3}) && sameSide2d(pt, p2, {p0, p1}) &&
        sameSide2d(pt, p3, {p1, p2}) && sameSide2d(pt, p1, {p3, p0}))
      return true;
  }
  else if constexpr (std::is_same_v<Elem, Tetrahedron>)
  {
    Vec3 const & p0 = e.pts[0]->coord;
    Vec3 const & p1 = e.pts[1]->coord;
    Vec3 const & p2 = e.pts[2]->coord;
    Vec3 const & p3 = e.pts[3]->coord;
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
    if (inside(Tetrahedron{{e.pts[0], e.pts[2], e.pts[5], e.pts[7]}, idNotSet}, pt))
      return true;
    else if (inside(
                 Tetrahedron{{e.pts[1], e.pts[2], e.pts[0], e.pts[5]}, idNotSet}, pt))
      return true;
    else if (inside(
                 Tetrahedron{{e.pts[3], e.pts[0], e.pts[2], e.pts[7]}, idNotSet}, pt))
      return true;
    else if (inside(
                 Tetrahedron{{e.pts[4], e.pts[5], e.pts[7], e.pts[0]}, idNotSet}, pt))
      return true;
    else if (inside(
                 Tetrahedron{{e.pts[6], e.pts[7], e.pts[5], e.pts[2]}, idNotSet}, pt))
      return true;

    // // minimum cost: 6 checks (1 per face)
    // // split each face in 2: 12 checks
    // // split each face in 8: 48 checks
    // Vec3 const & p0 = e.pts[0]->coord;
    // Vec3 const & p1 = e.pts[1]->coord;
    // Vec3 const & p2 = e.pts[2]->coord;
    // Vec3 const & p3 = e.pts[3]->coord;
    // Vec3 const & p4 = e.pts[0]->coord;
    // Vec3 const & p5 = e.pts[1]->coord;
    // Vec3 const & p6 = e.pts[2]->coord;
    // Vec3 const & p7 = e.pts[3]->coord;
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
