#pragma once

#include "def.hpp"

struct RefLineP1
{
  static uint const numPts = 2U;
  typedef Line GeoElem_T;
  typedef Eigen::Vector2d localVec_T;
  typedef Eigen::Matrix2d localMat_T;
  typedef std::array<Point*,numPts> pointList_T;

  static std::array<scalarFun_T,numPts> const phiFun;
  static std::array<vectorFun_T,numPts> const dphiFun;
  static localMat_T const massMat;
  static localMat_T const gradMat;
  static double constexpr volume = 2.L;
};

std::array<scalarFun_T,RefLineP1::numPts> const RefLineP1::phiFun =
{
  [] (Vec3 const & p) { return 0.5*(1-p(0)); },
  [] (Vec3 const & p) { return 0.5*(1+p(0)); }
};

std::array<vectorFun_T,RefLineP1::numPts> const RefLineP1::dphiFun =
{
  [] (Vec3 const & p) { return Vec3( 0.5L, 0.0, 0.0); },
  [] (Vec3 const & p) { return Vec3(-0.5L, 0.0, 0.0); }
};

RefLineP1::localMat_T const RefLineP1::massMat =
  (Eigen::Matrix2d() << 2.L/3, 1.L/3,
                        1.L/3, 2.L/3 ).finished();
RefLineP1::localMat_T const RefLineP1::gradMat =
  (Eigen::Matrix2d() <<  0.5L, -0.5L,
                        -0.5L,  0.5L ).finished();
