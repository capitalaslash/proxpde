#pragma once

#include "def.hpp"
#include "geo.hpp"
#include "reffe.hpp"

template <typename GeoElem, int N>
struct GaussQR
{
  using Vec_T = FVec<GeoElem::dim>;
  static int const numPts = N;

  static FMat<N,1> const weight;
  static array<Vec_T,N> const node;
};

template<> FVec<0> const GaussQR<NullElem,0>::weight = {};
template<> array<GaussQR<NullElem,0>::Vec_T,0> const GaussQR<NullElem,0>::node = {};

template<> FVec<1> const GaussQR<Line,1>::weight =
  FVec<1>::Constant(2.L);
template<> array<GaussQR<Line,1>::Vec_T,1> const GaussQR<Line,1>::node =
{{
  GaussQR<Line,1>::Vec_T::Constant(0.0L)
}};

static double constexpr sqrt13rd = 0.5773502691896258L;
template<> FVec<2> const GaussQR<Line,2>::weight = {1.L, 1.L};
template<> array<GaussQR<Line,2>::Vec_T,2> const GaussQR<Line,2>::node =
{{
  GaussQR<Line,2>::Vec_T::Constant(-sqrt13rd),
  GaussQR<Line,2>::Vec_T::Constant( sqrt13rd),
}};

static long double constexpr sqrt35th = 0.774596669241483L;
template<> FVec<3> const GaussQR<Line,3>::weight = {5.L/9, 8.L/9, 5.L/9};
template<> array<GaussQR<Line,3>::Vec_T,3> const GaussQR<Line,3>::node =
{{
  GaussQR<Line,2>::Vec_T::Constant(-sqrt35th),
  GaussQR<Line,2>::Vec_T::Constant( 0.0L),
  GaussQR<Line,2>::Vec_T::Constant( sqrt35th)
}};

template<> FVec<3> const GaussQR<Triangle,3>::weight =
    FVec<3>::Constant(1.L/6);
template<> array<GaussQR<Triangle,3>::Vec_T,3> const GaussQR<Triangle,3>::node =
{{
  GaussQR<Triangle,3>::Vec_T(0.5L, 0.0L),
  GaussQR<Triangle,3>::Vec_T(0.5L, 0.5L),
  GaussQR<Triangle,3>::Vec_T(0.0L, 0.5L)
}};

template<> FVec<4> const GaussQR<Triangle,4>::weight =
    (FVec<4>() <<
     25.L/96, 25.L/96, 25.L/96, -27.L/96).finished();
template<> array<GaussQR<Triangle,4>::Vec_T,4> const GaussQR<Triangle,4>::node =
{{
  GaussQR<Triangle,4>::Vec_T(0.2L, 0.2L),
  GaussQR<Triangle,4>::Vec_T(0.6L, 0.2L),
  GaussQR<Triangle,4>::Vec_T(0.2L, 0.6L),
  GaussQR<Triangle,4>::Vec_T(1.L/3., 1.L/3.)
}};

template<> FVec<9> const GaussQR<Quad,9>::weight =
    (FVec<9>() <<
     25.L/81, 40.L/81, 25.L/81,
     40.L/81, 64.L/81, 40.L/81,
     25.L/81, 40.L/81, 25.L/81).finished();
template<> array<GaussQR<Quad,9>::Vec_T,9> const GaussQR<Quad,9>::node =
{{
  GaussQR<Quad,9>::Vec_T(-sqrt35th, -sqrt35th),
  GaussQR<Quad,9>::Vec_T(       0., -sqrt35th),
  GaussQR<Quad,9>::Vec_T( sqrt35th, -sqrt35th),
  GaussQR<Quad,9>::Vec_T(-sqrt35th,        0.),
  GaussQR<Quad,9>::Vec_T(       0.,        0.),
  GaussQR<Quad,9>::Vec_T( sqrt35th,        0.),
  GaussQR<Quad,9>::Vec_T(-sqrt35th,  sqrt35th),
  GaussQR<Quad,9>::Vec_T(       0.,  sqrt35th),
  GaussQR<Quad,9>::Vec_T( sqrt35th,  sqrt35th)
}};

template <typename GeoElem>
struct TrapQR
{
  using Vec_T = FVec<GeoElem::dim>;
  static uint const numPts = GeoElem::numPts;

  static FVec<GeoElem::numPts> const weight;
  static array<Vec_T,GeoElem::numPts> const node;
};

template<> FVec<2> const TrapQR<Line>::weight =
  FVec<2>::Constant(1.L);
template<> array<TrapQR<Line>::Vec_T,2> const TrapQR<Line>::node =
{{
  TrapQR<Line>::Vec_T::Constant(-1.L),
  TrapQR<Line>::Vec_T::Constant( 1.L)
}};

template<> FVec<3> const TrapQR<Triangle>::weight =
  FVec<3>::Constant(1.L/6);
template<> array<TrapQR<Triangle>::Vec_T,3> const TrapQR<Triangle>::node =
{{
  TrapQR<Triangle>::Vec_T( 0., 0.),
  TrapQR<Triangle>::Vec_T( 1., 0.),
  TrapQR<Triangle>::Vec_T( 0., 1.)
}};

template<> FVec<4> const TrapQR<Quad>::weight =
  FVec<4>::Constant(1.L);
template<> array<TrapQR<Quad>::Vec_T,4> const TrapQR<Quad>::node =
{{
  TrapQR<Quad>::Vec_T(-1.,-1.),
  TrapQR<Quad>::Vec_T(-1., 1.),
  TrapQR<Quad>::Vec_T( 1.,-1.),
  TrapQR<Quad>::Vec_T( 1., 1.)
}};

template <typename QR>
struct SideQR
{
  using Type = GaussQR<NullElem,0>;
};

template <>
struct SideQR<GaussQR<Line,3>>
{
  using Type = GaussQR<PointElem,1>;
};

template <>
struct SideQR<GaussQR<Triangle,4>>
{
  using Type = GaussQR<Line,3>;
};

template <>
struct SideQR<GaussQR<Quad,9>>
{
  using Type = GaussQR<Line,3>;
};
