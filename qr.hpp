#pragma once

#include "def.hpp"
#include "geo.hpp"
#include "reffe.hpp"

template <typename GeoElem, uint N>
struct GaussQR
{
  static uint const numPts = N;

  static Eigen::Array<double,N,1> const w;
  static std::array<Vec3,N> const n;
};

template<> Eigen::Array<double,1,1> const GaussQR<Line,1>::w =
  Eigen::Array<double,1,1>::Constant(2.L);
template<> std::array<Vec3,1> const GaussQR<Line,1>::n =
{
  Vec3( 0., 0., 0.),
};

static double constexpr sqrt13rd = std::sqrt(1.L/3);
template<> Eigen::Array<double,2,1> const GaussQR<Line,2>::w = {1.L, 1.L};
template<> std::array<Vec3,2> const GaussQR<Line,2>::n =
{
  Vec3(-sqrt13rd, 0., 0.),
  Vec3( sqrt13rd, 0., 0.)
};

static double constexpr sqrt35th = std::sqrt(3.L/5);
template<> Eigen::Array<double,3,1> const GaussQR<Line,3>::w = {5.L/9, 8.L/9, 5.L/9};
template<> std::array<Vec3,3> const GaussQR<Line,3>::n =
{
  Vec3(-sqrt35th, 0., 0.),
  Vec3(       0., 0., 0.),
  Vec3( sqrt35th, 0., 0.)
};
