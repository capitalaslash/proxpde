#pragma once

#include <array>
#include <vector>
#include <functional>

#include "def.hpp"

class Point
{
public:
  explicit Point(Vec3 const & c = Vec3::Zero(3),
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

class Line
{
public:
  static uint const numPts = 2;
  typedef Eigen::Vector2d localVec_T;
  typedef Eigen::Matrix2d localMat_T;
  typedef std::array<Point*,numPts> pointList_T;

  explicit Line(pointList_T const&& pl = {nullptr, nullptr}):
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

  pointList_T pointList;
  static std::array<scalarFun_T,numPts> const shapeFuns;
  static localMat_T const reactMat;
  static localMat_T const gradMat;
  static constexpr double _refVolume = 2.L;
};

std::array<scalarFun_T,Line::numPts> const Line::shapeFuns =
{
  [] (Vec3 const & p) { return 0.5*(1-p(0)); },
  [] (Vec3 const & p) { return 0.5*(1+p(0)); }
};

Line::localMat_T const Line::reactMat =
  (Eigen::Matrix2d() << 2.L/3, 1.L/3,
                        1.L/3, 2.L/3 ).finished();
Line::localMat_T const Line::gradMat =
  (Eigen::Matrix2d() <<  0.5L, -0.5L,
                        -0.5L,  0.5L ).finished();

template <class Elem>
class Mesh
{
public:
  typedef std::vector<Point> pointList_T;
  typedef std::vector<Line> elementList_T;
  typedef std::array<id_T,Elem::numPts> elementConn_T;
  typedef std::vector<std::array<id_T,Elem::numPts>> connList_T;

  void buildConnectivity()
  {
    _connList.reserve(elementList.size());
    for(auto& l: elementList)
    {
      elementConn_T elemConn;
      uint counter = 0;
      for(auto& p: l.pointList)
      {
        elemConn[counter] = p->id;
        counter++;
      }
      _connList.push_back(elemConn);
    }
  }

  pointList_T pointList;
  elementList_T elementList;
  connList_T _connList;
};

enum side
{
  LEFT,
  RIGHT
};

void buildMesh1D(std::shared_ptr<Mesh<Line>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 uint const numPts)
{
  Vec3 const h = length / (numPts-1);
  meshPtr->pointList.reserve(numPts);
  for(uint p=0; p<numPts; ++p)
  {
    meshPtr->pointList.emplace_back(origin + p * h, p);
  }
  meshPtr->pointList[0].marker = side::LEFT;
  meshPtr->pointList[numPts-1].marker = side::RIGHT;

  uint const numElems = numPts-1;
  meshPtr->elementList.reserve(numElems);
  for(uint e=0; e<numElems; ++e)
  {
    meshPtr->elementList.emplace_back(
      Line::pointList_T{&meshPtr->pointList[e], &meshPtr->pointList[e+1]});
  }

  meshPtr->buildConnectivity();
}
