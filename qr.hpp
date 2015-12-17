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

template<> Eigen::Array<double,3,1> const GaussQR<Triangle,3>::w =
    Eigen::Array<double,3,1>::Constant(1.L/6);
template<> std::array<Vec3,3> const GaussQR<Triangle,3>::n =
{
  Vec3(0.5, 0.0, 0.),
  Vec3(0.5, 0.5, 0.),
  Vec3(0.0, 0.5, 0.)
};

template<> Eigen::Array<double,9,1> const GaussQR<Quad,9>::w = (Eigen::Array<double,9,1>() <<
  25.L/81, 40.L/81, 25.L/81, 40.L/81, 64.L/81, 40.L/81, 25.L/81, 40.L/81, 25.L/81).finished();

template<> std::array<Vec3,9> const GaussQR<Quad,9>::n =
{
  Vec3(-sqrt35th, -sqrt35th, 0.),
  Vec3(       0., -sqrt35th, 0.),
  Vec3( sqrt35th, -sqrt35th, 0.),
  Vec3(-sqrt35th,        0., 0.),
  Vec3(       0.,        0., 0.),
  Vec3( sqrt35th,        0., 0.),
  Vec3(-sqrt35th,  sqrt35th, 0.),
  Vec3(       0.,  sqrt35th, 0.),
  Vec3( sqrt35th,  sqrt35th, 0.)
};

template <typename GeoElem>
struct TrapQR
{
  static uint const numPts = GeoElem::numPts;

  static Eigen::Array<double,GeoElem::numPts,1> const w;
  static std::array<Vec3,GeoElem::numPts> const n;
};

template<> Eigen::Array<double,2,1> const TrapQR<Line>::w =
  Eigen::Array<double,2,1>::Constant(1.L);
template<> std::array<Vec3,2> const TrapQR<Line>::n =
{
  Vec3(-1., 0., 0.),
  Vec3( 1., 0., 0.)
};

