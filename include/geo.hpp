#pragma once

#include "def.hpp"

#include <algorithm>

class Point
{
public:
  Point(Vec3 const c,
        id_T const i,
        marker_T const m = markerNotSet):
    coord(std::move(c)),
    id(i),
    marker(m)
  {}

  Point() = default;

  double operator[](uint const i) const
  {
    return this->coord[i];
  }

  Vec3 coord  = Vec3::Zero();
  id_T id = dofIdNotSet;
  marker_T marker = markerNotSet;
};

inline std::ostream& operator<<(std::ostream& out, Point const & p)
{
  out << "(" << p[0] << "," << p[1] << "," << p[2] << "), id: "
      << p.id << ", m: " << p.marker;
  return out;
}

struct GeoElem
{
  using PointList_T = std::vector<Point *>;
  struct FacingElem
  {
    GeoElem const * ptr;
    id_T side;
  };

  GeoElem(std::initializer_list<Point *> const & pList,
          id_T const i,
          marker_T const m):
    pointList{pList},
    id{i},
    marker{m}
  {}

  GeoElem(PointList_T const pList,
          id_T const i,
          marker_T const m):
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
  bool onBoundary() const
  {
    // check that facingElems have been initialized
    assert(facingElem[0].ptr != nullptr);
    // the facet is on the boundary iff there is no outside element
    return facingElem[1].ptr == nullptr;
  }

  PointList_T pointList = {};
  id_T id = idNotSet;
  marker_T marker = markerNotSet;
  // stores internal and external (element, facet side) pairs
  array<FacingElem, 2> facingElem = {FacingElem{nullptr, idNotSet}, FacingElem{nullptr, idNotSet}};
  Vec3 _normal;
};

inline bool operator< (GeoElem const & e1, GeoElem const & e2)
{
  return e1.id < e2.id;
}

inline std::ostream & operator<<(std::ostream & out, GeoElem const & e)
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

  PointElem(std::initializer_list<Point *> const & pList,
            id_T const i,
            marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  PointElem(PointList_T const & pList,
            id_T const i,
            marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  PointElem() = default;

  ~PointElem() override = default;

  Vec3 midpoint() const final {return pointList[0]->coord;}
  Vec3 origin() const final {return pointList[0]->coord;}
  double volume() const final {return 1.;}
  void buildNormal() final
  {
    _normal = Vec3{1.0, 0.0, 0.0};
  }

  Vec3 normal() const final
  {
    return Vec3{1.0, 0.0, 0.0};
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

  explicit Line(std::initializer_list<Point *> const & pList,
                id_T const i,
                marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Line(PointList_T const & pList,
       id_T const i,
       marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Line() = default;

  ~Line() override = default;

  Vec3 midpoint() const final
  {
    return Vec3(0.5*(pointList[1]->coord + pointList[0]->coord));
  }

  Vec3 origin() const final
  {
    return this->midpoint();
  }

  double volume() const final
  {
    return (pointList[1]->coord - pointList[0]->coord).norm();
  }

  void buildNormal() final
  {
    Vec3 const length =  pointList[1]->coord - pointList[0]->coord;
    _normal = Vec3{length[1], -length[0], 0.0};
    _normal.normalize();
  }

  Vec3 normal() const final
  {
    Vec3 const length =  pointList[1]->coord - pointList[0]->coord;
    auto n = Vec3{length[1], -length[0], 0.0};
    n.normalize();
    return n;
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

  Triangle(std::initializer_list<Point *> const & pList,
           id_T const i,
           marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Triangle(PointList_T const & pList,
           id_T const i,
           marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Triangle() = default;

  ~Triangle() override = default;

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

  double h_min() const
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[1]->coord;
    auto const v3 = pointList[0]->coord - pointList[2]->coord;

    // triangle size based on diameter of the inscribed circle
    return 4 * volume() / (v1.norm() + v2.norm() + v3.norm());
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

  Quad(std::initializer_list<Point *> const & pList,
       id_T const i,
       marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Quad(PointList_T const & pList,
       id_T const i,
       marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Quad() = default;

  ~Quad() override = default;

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

  void buildNormal() final
  {
    // TODO: we consider only planar quads
    auto const v1 = pointList[1]->coord-pointList[0]->coord;
    auto const v2 = pointList[2]->coord-pointList[0]->coord;
    _normal = v1.cross(v2);
    _normal.normalize();
  }

  Vec3 normal() const final
  {
    // TODO: we consider only planar quads
    auto const v1 = pointList[1]->coord-pointList[0]->coord;
    auto const v2 = pointList[2]->coord-pointList[0]->coord;
    auto n = v1.cross(v2);
    n.normalize();
    return n;
  }
};

class Tetrahedron: public GeoElem
{
public:
  using Facet_T = Triangle;
  using Face_T = Triangle;
  using Edge_T = Line;
  static int constexpr dim = 3;
  static uint constexpr numPts = 4U;
  static uint constexpr numEdges = 6U;
  static uint constexpr numFaces = 4U;
  static uint constexpr numFacets = 4U;
  static array<array<id_T,3>,4> constexpr elemToFacet = {{
    {{0,2,1}}, {{0,1,3}}, {{0,3,2}}, {{1,2,3}}
  }};
  static array<array<id_T,2>,6> constexpr elemToEdge = {{
    {{0,1}}, {{1,2}}, {{2,0}}, {{0,3}}, {{1,3}}, {{2,3}}
  }};
  static array<array<id_T,3>,4> constexpr elemToFace = elemToFacet;

  Tetrahedron(std::initializer_list<Point *> const & pList,
              id_T const i,
              marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Tetrahedron(PointList_T const & pList,
              id_T const i,
              marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Tetrahedron() = default;

  ~Tetrahedron() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{.25*(
      pointList[0]->coord + pointList[1]->coord +
      pointList[2]->coord + pointList[3]->coord)
    };
  }

  Vec3 origin() const final
  {
    return pointList[0]->coord;
  }

  double volume() const final
  {
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[2]->coord - pointList[0]->coord;
    auto const v3 = pointList[3]->coord - pointList[0]->coord;
    return (v1.cross(v2)).dot(v3) / 6.;
  }

  void buildNormal() final
  {
    std::abort();
  }

  Vec3 normal() const final
  {
    std::abort();
  }
};

class Hexahedron: public GeoElem
{
public:
  using Facet_T = Quad;
  using Face_T = Quad;
  using Edge_T = Line;
  static int constexpr dim = 3;
  static uint constexpr numPts = 8U;
  static uint constexpr numEdges = 12U;
  static uint constexpr numFaces = 6U;
  static uint constexpr numFacets = 6U;
  static array<array<id_T,4>,6> constexpr elemToFacet = {{
    {{0,3,2,1}}, {{0,1,5,4}}, {{1,2,6,5}}, {{2,3,7,6}}, {{3,0,4,7}}, {{4,5,6,7}}
  }};
  static array<array<id_T,2>,12> constexpr elemToEdge = {{
    {{0,1}}, {{1,2}}, {{2,3}}, {{3,0}}, {{0,4}}, {{1,5}}, {{2,6}}, {{3,7}}, {{4,5}}, {{5,6}}, {{6,7}}, {{7,4}}
  }};
  static array<array<id_T,4>,6> constexpr elemToFace = elemToFacet;

  Hexahedron(std::initializer_list<Point *> const & pList,
             id_T const i,
             marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Hexahedron(PointList_T const & pList,
             id_T const i,
             marker_T const m = markerNotSet):
    GeoElem{pList, i, m}
  {}

  Hexahedron() = default;

  ~Hexahedron() override = default;

  Vec3 midpoint() const final
  {
    return Vec3{.125*(
      pointList[0]->coord + pointList[1]->coord +
      pointList[2]->coord + pointList[3]->coord +
      pointList[4]->coord + pointList[5]->coord +
      pointList[6]->coord + pointList[7]->coord)
    };
  }

  Vec3 origin() const final
  {
    return this->midpoint();
  }

  double volume() const final
  {
    // TODO: this is not correct for general hexahedrons, just for parallelepids
    // a betetr way is to split the hexahedron in ettrahedra and sum their volumes
    auto const v1 = pointList[1]->coord - pointList[0]->coord;
    auto const v2 = pointList[3]->coord - pointList[0]->coord;
    auto const v3 = pointList[4]->coord - pointList[0]->coord;
    return (v1.cross(v2)).dot(v3);
  }

  void buildNormal() final
  {
    std::abort();
  }

  Vec3 normal() const final
  {
    std::abort();
  }
};


// this method checks only to see if the 2 elems are equivalent from the geometric pov, i.e. they
// have the same point ids (or a permutation of it)
template <typename Elem>
bool geoEqual(Elem const & e1, Elem const & e2)
{
  array<id_T, Elem::numPts> ids1, ids2;
  uint counter= 0;
  std::for_each(e1.pointList.begin(), e1.pointList.end(), [&ids1, &counter](Point const * p)
  {
    ids1[counter++] = p->id;
  });
  counter = 0;
  std::for_each(e2.pointList.begin(), e2.pointList.end(), [&ids2, &counter](Point const * p)
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
