#pragma once

#include "def.hpp"
#include "geo.hpp"

struct RefLineP1
{
  static uint const numPts = 2U;
  typedef Line GeoElem_T;
  typedef Eigen::Vector2d LocalVec_T;
  typedef Eigen::Matrix2d LocalMat_T;
  typedef std::array<Point*,numPts> PointList_T;

  static std::array<Vec3,numPts> const points;
  static std::array<scalarFun_T,numPts> const phiFun;
  static std::array<vectorFun_T,numPts> const dphiFun;
  static LocalMat_T const massMat;
  static LocalMat_T const gradMat;
  static double constexpr volume = 2.L;
};

struct RefLineP2
{
  static uint const numPts = 3U;
  typedef Line GeoElem_T;
  typedef Eigen::Vector3d LocalVec_T;
  typedef Eigen::Matrix3d LocalMat_T;
  typedef std::array<Point*,numPts> PointList_T;

  static std::array<Vec3,numPts> const points;
  static std::array<scalarFun_T,numPts> const phiFun;
  static std::array<vectorFun_T,numPts> const dphiFun;
  static LocalMat_T const massMat;
  static LocalMat_T const gradMat;
  static double constexpr volume = 2.L;
};
