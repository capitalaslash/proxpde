#pragma once

#include "def.hpp"
#include "geo.hpp"
#include "reffe.hpp"

template <typename GeoElem, int N>
struct GaussQR
{
  using Scalar_T = double;
  using GeoElem_T = GeoElem;
  using Vec_T = FVec<GeoElem_T::dim>;
  using Weights_T = FVec<N>;
  // using Vec_T = Eigen::Matrix<long double, GeoElem::dim, 1>;
  // using Weights_T = Eigen::Matrix<long double, N, 1>;
  static int const numPts = N;

  static Weights_T const weight;
  static array<Vec_T,N> const node;
};

// ----------------------------------------------------------------------------
template<> GaussQR<NullElem, 0>::Weights_T const GaussQR<NullElem, 0>::weight = {};
template<> array<GaussQR<NullElem, 0>::Vec_T, 0> const GaussQR<NullElem, 0>::node = {};

// ----------------------------------------------------------------------------
template<> GaussQR<Line, 1>::Weights_T const GaussQR<Line, 1>::weight =
    GaussQR<Line, 1>::Weights_T::Constant(2.L);
template<> array<GaussQR<Line, 1>::Vec_T,1> const GaussQR<Line, 1>::node =
{{
   GaussQR<Line, 1>::Vec_T::Constant(0.0L)
}};

static GaussQR<NullElem, 0>::Scalar_T constexpr sqrt13rd = 0.5773502691896258L;
template<> GaussQR<Line, 2>::Weights_T const GaussQR<Line, 2>::weight = {1.L, 1.L};
template<> array<GaussQR<Line, 2>::Vec_T, 2> const GaussQR<Line, 2>::node =
{{
   GaussQR<Line,2>::Vec_T::Constant(-sqrt13rd),
   GaussQR<Line,2>::Vec_T::Constant( sqrt13rd),
}};

static GaussQR<NullElem, 0>::Scalar_T constexpr sqrt35th = 0.774596669241483L;
template<> GaussQR<Line, 3>::Weights_T const GaussQR<Line, 3>::weight = {5.L/9, 8.L/9, 5.L/9};
template<> array<GaussQR<Line,3>::Vec_T, 3> const GaussQR<Line, 3>::node =
{{
   GaussQR<Line, 3>::Vec_T::Constant(-sqrt35th),
   GaussQR<Line, 3>::Vec_T::Constant( 0.0L),
   GaussQR<Line, 3>::Vec_T::Constant( sqrt35th)
}};

static GaussQR<NullElem, 0>::Scalar_T constexpr sqrt37thm = 0.3399810435848563L; // sqrt((3.L - 2.L * sqrt(1.2L)) / 7);
static GaussQR<NullElem, 0>::Scalar_T constexpr sqrt37thp = 0.8611363115940526L; // sqrt((3.L + 2.L * sqrt(1.2L)) / 7);
template<> GaussQR<Line, 4>::Weights_T const GaussQR<Line, 4>::weight = {
  (18.L - sqrt(30.L))/36, (18.L + sqrt(30.L))/36, (18.L + sqrt(30.L))/36, (18.L - sqrt(30.L))/36
};
template<> array<GaussQR<Line, 4>::Vec_T, 4> const GaussQR<Line, 4>::node =
{{
   GaussQR<Line, 4>::Vec_T::Constant(-sqrt37thp),
   GaussQR<Line, 4>::Vec_T::Constant(-sqrt37thm),
   GaussQR<Line, 4>::Vec_T::Constant( sqrt37thm),
   GaussQR<Line, 4>::Vec_T::Constant( sqrt37thp),
}};

// ----------------------------------------------------------------------------
template<> GaussQR<Triangle, 1>::Weights_T const GaussQR<Triangle, 1>::weight =
    GaussQR<Triangle, 1>::Weights_T::Constant(0.5L);
template<> array<GaussQR<Triangle, 1>::Vec_T,1> const GaussQR<Triangle, 1>::node =
{{
   GaussQR<Triangle, 1>::Vec_T(1.L/3, 1.L/3)
}};

template<> GaussQR<Triangle, 3>::Weights_T const GaussQR<Triangle, 3>::weight =
    GaussQR<Triangle, 3>::Weights_T::Constant(1.L/6);
template<> array<GaussQR<Triangle, 3>::Vec_T, 3> const GaussQR<Triangle,3>::node =
{{
   GaussQR<Triangle, 3>::Vec_T(1./6, 1./6),
   GaussQR<Triangle, 3>::Vec_T(2./3, 1./6),
   GaussQR<Triangle, 3>::Vec_T(1./6, 2./3)
 }};
// this node positioning seems to give smaller errors in some cases but it is
// not ideal since they coincide with the location of the dofs of the P2 element
// template<> array<GaussQR<Triangle,3>::Vec_T,3> const GaussQR<Triangle,3>::node =
// {{
//     GaussQR<Triangle, 3>::Vec_T(0.5L, 0.0L),
//     GaussQR<Triangle, 3>::Vec_T(0.5L, 0.5L),
//     GaussQR<Triangle, 3>::Vec_T(0.0L, 0.5L)
// }};

template<> GaussQR<Triangle, 4>::Weights_T const GaussQR<Triangle, 4>::weight =
    (GaussQR<Triangle, 4>::Weights_T() <<
     25.L, 25.L, 25.L, -27.L).finished() / 96.L;
template<> array<GaussQR<Triangle, 4>::Vec_T, 4> const GaussQR<Triangle, 4>::node =
{{
   GaussQR<Triangle, 4>::Vec_T(0.2L, 0.2L),
   GaussQR<Triangle, 4>::Vec_T(0.6L, 0.2L),
   GaussQR<Triangle, 4>::Vec_T(0.2L, 0.6L),
   GaussQR<Triangle, 4>::Vec_T(1.L/3, 1.L/3)
 }};

template<> GaussQR<Triangle, 6>::Weights_T const GaussQR<Triangle, 6>::weight = .5*
    (GaussQR<Triangle, 6>::Weights_T() <<
     0.109951743655322L, 0.109951743655322L, 0.109951743655322L,
     0.223381589678011L, 0.223381589678011L, 0.223381589678011L).finished();
template<> array<GaussQR<Triangle, 6>::Vec_T, 6> const GaussQR<Triangle, 6>::node =
{{
   GaussQR<Triangle, 6>::Vec_T(0.091576213509771L, 0.091576213509771L),
   GaussQR<Triangle, 6>::Vec_T(0.816847572980459L, 0.091576213509771L),
   GaussQR<Triangle, 6>::Vec_T(0.091576213509771L, 0.816847572980459L),
   GaussQR<Triangle, 6>::Vec_T(0.445948490915965L, 0.445948490915965L),
   GaussQR<Triangle, 6>::Vec_T(0.108103018168070L, 0.445948490915965L),
   GaussQR<Triangle, 6>::Vec_T(0.445948490915965L, 0.108103018168070L),
 }};

template<> GaussQR<Triangle, 7>::Weights_T const GaussQR<Triangle, 7>::weight = .5*
    (GaussQR<Triangle, 7>::Weights_T() <<
     0.22500000000000L,
     0.13239415278851L, 0.13239415278851L, 0.13239415278851L,
     0.12593918054483L, 0.12593918054483L, 0.12593918054483L).finished();
template<> array<GaussQR<Triangle, 7>::Vec_T, 7> const GaussQR<Triangle, 7>::node =
{{
   GaussQR<Triangle,7>::Vec_T(1.L/3, 1.L/3),
   GaussQR<Triangle,7>::Vec_T(0.47014206410511L, 0.47014206410511L),
   GaussQR<Triangle,7>::Vec_T(0.47014206410511L, 0.05971587178977L),
   GaussQR<Triangle,7>::Vec_T(0.05971587178977L, 0.47014206410511L),
   GaussQR<Triangle,7>::Vec_T(0.10128650732346L, 0.10128650732346L),
   GaussQR<Triangle,7>::Vec_T(0.10128650732346L, 0.79742698535309L),
   GaussQR<Triangle,7>::Vec_T(0.79742698535309L, 0.10128650732346L),
 }};

// ----------------------------------------------------------------------------
template<> GaussQR<Quad, 1>::Weights_T const GaussQR<Quad, 1>::weight =
    GaussQR<Quad, 1>::Weights_T::Constant(4.L);
template<> array<GaussQR<Quad, 1>::Vec_T, 1> const GaussQR<Quad, 1>::node =
{{
   GaussQR<Quad, 1>::Vec_T(0.L, 0.L)
}};

template<> GaussQR<Quad, 4>::Weights_T const GaussQR<Quad, 4>::weight =
    GaussQR<Quad, 4>::Weights_T::Constant(1.L);
template<> array<GaussQR<Quad, 4>::Vec_T, 4> const GaussQR<Quad, 4>::node =
{{
   GaussQR<Quad, 4>::Vec_T(-sqrt13rd, -sqrt13rd),
   GaussQR<Quad, 4>::Vec_T( sqrt13rd, -sqrt13rd),
   GaussQR<Quad, 4>::Vec_T(-sqrt13rd,  sqrt13rd),
   GaussQR<Quad, 4>::Vec_T( sqrt13rd,  sqrt13rd),
}};

template<> GaussQR<Quad, 9>::Weights_T const GaussQR<Quad, 9>::weight =
    (GaussQR<Quad, 9>::Weights_T() <<
     25.L, 40.L, 25.L,
     40.L, 64.L, 40.L,
     25.L, 40.L, 25.L).finished() / 81.L;
template<> array<GaussQR<Quad, 9>::Vec_T, 9> const GaussQR<Quad, 9>::node =
{{
   GaussQR<Quad, 9>::Vec_T(-sqrt35th, -sqrt35th),
   GaussQR<Quad, 9>::Vec_T(       0., -sqrt35th),
   GaussQR<Quad, 9>::Vec_T( sqrt35th, -sqrt35th),
   GaussQR<Quad, 9>::Vec_T(-sqrt35th,        0.),
   GaussQR<Quad, 9>::Vec_T(       0.,        0.),
   GaussQR<Quad, 9>::Vec_T( sqrt35th,        0.),
   GaussQR<Quad, 9>::Vec_T(-sqrt35th,  sqrt35th),
   GaussQR<Quad, 9>::Vec_T(       0.,  sqrt35th),
   GaussQR<Quad, 9>::Vec_T( sqrt35th,  sqrt35th)
}};

// ----------------------------------------------------------------------------
template<> GaussQR<Tetrahedron, 1>::Weights_T const GaussQR<Tetrahedron,1>::weight =
    GaussQR<Tetrahedron, 1>::Weights_T::Constant(1.L/6);
template<> array<GaussQR<Tetrahedron, 1>::Vec_T, 1> const GaussQR<Tetrahedron, 1>::node =
{{
  GaussQR<Tetrahedron, 1>::Vec_T(.25L, .25L, .25L)
}};

template<> GaussQR<Tetrahedron, 4>::Weights_T const GaussQR<Tetrahedron, 4>::weight =
    GaussQR<Tetrahedron, 4>::Weights_T::Constant(.25L/6);
template<> array<GaussQR<Tetrahedron, 4>::Vec_T, 4> const GaussQR<Tetrahedron, 4>::node =
{{
   GaussQR<Tetrahedron, 4>::Vec_T(0.1381966011250105L, 0.1381966011250105L, 0.1381966011250105L),
   GaussQR<Tetrahedron, 4>::Vec_T(0.5854101966249685L, 0.1381966011250105L, 0.1381966011250105L),
   GaussQR<Tetrahedron, 4>::Vec_T(0.1381966011250105L, 0.5854101966249685L, 0.1381966011250105L),
   GaussQR<Tetrahedron, 4>::Vec_T(0.1381966011250105L, 0.1381966011250105L, 0.5854101966249685L)
}};

template<> GaussQR<Tetrahedron, 11>::Weights_T const GaussQR<Tetrahedron, 11>::weight =
    (GaussQR<Tetrahedron, 11>::Weights_T() <<
     -0.0131555555555556L,
     0.0076222222222222L, 0.0076222222222222L, 0.0076222222222222L, 0.0076222222222222L,
     0.024888888888889L, 0.024888888888889L, 0.024888888888889L, 0.024888888888889L, 0.024888888888889L, 0.024888888888889L
     ).finished();
template<> array<GaussQR<Tetrahedron, 11>::Vec_T, 11> const GaussQR<Tetrahedron, 11>::node =
{{
   GaussQR<Tetrahedron, 11>::Vec_T(0.25L, 0.25L, 0.25L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.0714285714285714L, 0.0714285714285714L, 0.0714285714285714L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.7857142857142857L, 0.0714285714285714L, 0.0714285714285714L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.0714285714285714L, 0.7857142857142857L, 0.0714285714285714L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.0714285714285714L, 0.0714285714285714L, 0.7857142857142857L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.399403576166799L, 0.399403576166799L, 0.100596423833201L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.399403576166799L, 0.100596423833201L, 0.399403576166799L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.100596423833201L, 0.399403576166799L, 0.399403576166799L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.399403576166799L, 0.100596423833201L, 0.100596423833201L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.100596423833201L, 0.399403576166799L, 0.100596423833201L),
   GaussQR<Tetrahedron, 11>::Vec_T(0.100596423833201L, 0.100596423833201L, 0.399403576166799L),
}};

// ----------------------------------------------------------------------------
template<> GaussQR<Hexahedron, 1>::Weights_T const GaussQR<Hexahedron, 1>::weight =
    GaussQR<Hexahedron, 1>::Weights_T::Constant(8.L);
template<> array<GaussQR<Hexahedron, 1>::Vec_T, 1> const GaussQR<Hexahedron, 1>::node =
{{
   GaussQR<Hexahedron, 1>::Vec_T(0.L, 0.L, 0.L),
}};

template<> GaussQR<Hexahedron, 8>::Weights_T const GaussQR<Hexahedron, 8>::weight =
    GaussQR<Hexahedron, 8>::Weights_T::Constant(1.L);
template<> array<GaussQR<Hexahedron, 8>::Vec_T, 8> const GaussQR<Hexahedron, 8>::node =
{{
   GaussQR<Hexahedron, 8>::Vec_T(-sqrt13rd, -sqrt13rd, -sqrt13rd),
   GaussQR<Hexahedron, 8>::Vec_T( sqrt13rd, -sqrt13rd, -sqrt13rd),
   GaussQR<Hexahedron, 8>::Vec_T(-sqrt13rd,  sqrt13rd, -sqrt13rd),
   GaussQR<Hexahedron, 8>::Vec_T( sqrt13rd,  sqrt13rd, -sqrt13rd),
   GaussQR<Hexahedron, 8>::Vec_T(-sqrt13rd, -sqrt13rd,  sqrt13rd),
   GaussQR<Hexahedron, 8>::Vec_T( sqrt13rd, -sqrt13rd,  sqrt13rd),
   GaussQR<Hexahedron, 8>::Vec_T(-sqrt13rd,  sqrt13rd,  sqrt13rd),
   GaussQR<Hexahedron, 8>::Vec_T( sqrt13rd,  sqrt13rd,  sqrt13rd),
}};

template<> GaussQR<Hexahedron, 27>::Weights_T const GaussQR<Hexahedron, 27>::weight =
    (GaussQR<Hexahedron, 27>::Weights_T() <<
     125.L, 200.L, 125.L, 200.L, 320.L, 200.L, 125.L, 200.L, 125.L,
     200.L, 320.L, 200.L, 320.L, 512.L, 320.L, 200.L, 320.L, 200.L,
     125.L, 200.L, 125.L, 200.L, 320.L, 200.L, 125.L, 200.L, 125.L
     ).finished() / 729.L;
template<> array<GaussQR<Hexahedron, 27>::Vec_T, 27> const GaussQR<Hexahedron, 27>::node =
{{
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th, -sqrt35th, -sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L, -sqrt35th, -sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th, -sqrt35th, -sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th,       0.L, -sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L,       0.L, -sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th,       0.L, -sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th,  sqrt35th, -sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L,  sqrt35th, -sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th,  sqrt35th, -sqrt35th),
   //
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th, -sqrt35th,       0.L),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L, -sqrt35th,       0.L),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th, -sqrt35th,       0.L),
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th,       0.L,       0.L),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L,       0.L,       0.L),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th,       0.L,       0.L),
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th,  sqrt35th,       0.L),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L,  sqrt35th,       0.L),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th,  sqrt35th,       0.L),
   //
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th, -sqrt35th,  sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L, -sqrt35th,  sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th, -sqrt35th,  sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th,       0.L,  sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L,       0.L,  sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th,       0.L,  sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(-sqrt35th,  sqrt35th,  sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T(      0.L,  sqrt35th,  sqrt35th),
   GaussQR<Hexahedron, 27>::Vec_T( sqrt35th,  sqrt35th,  sqrt35th),
}};

// ----------------------------------------------------------------------------
template <typename GeoElem>
struct TrapQR
{
  using GeoElem_T = GeoElem;
  using Vec_T = FVec<GeoElem::dim>;
  static uint const numPts = GeoElem::numPts;

  static FVec<numPts> const weight;
  static array<Vec_T, numPts> const node;
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

template<> FVec<8> const TrapQR<Hexahedron>::weight =
  FVec<8>::Constant(1.L);
template<> array<TrapQR<Hexahedron>::Vec_T,8> const TrapQR<Hexahedron>::node =
{{
   TrapQR<Hexahedron>::Vec_T(-1.,-1., -1.),
   TrapQR<Hexahedron>::Vec_T(-1., 1., -1.),
   TrapQR<Hexahedron>::Vec_T( 1.,-1., -1.),
   TrapQR<Hexahedron>::Vec_T( 1., 1., -1.),
   TrapQR<Hexahedron>::Vec_T(-1.,-1.,  1.),
   TrapQR<Hexahedron>::Vec_T(-1., 1.,  1.),
   TrapQR<Hexahedron>::Vec_T( 1.,-1.,  1.),
   TrapQR<Hexahedron>::Vec_T( 1., 1.,  1.)
}};

template <typename GeoElem>
struct SimpsonQR
{
  using Vec_T = FVec<GeoElem::dim>;
  static uint const numPts = pow(3, GeoElem::dim);

  static FVec<numPts> const weight;
  static array<Vec_T, numPts> const node;
};

template<> FVec<3> const SimpsonQR<Line>::weight =
    (FVec<3>() << 1., 4., 1.).finished() / 3.;
template<> array<SimpsonQR<Line>::Vec_T, 3> const SimpsonQR<Line>::node =
{{
   SimpsonQR<Line>::Vec_T::Constant(-1.L),
   SimpsonQR<Line>::Vec_T::Constant( 0.L),
   SimpsonQR<Line>::Vec_T::Constant( 1.L)
}};

template <typename QR>
struct SideQR {};

template <typename QR>
using SideQR_T = typename SideQR<QR>::type;

template <>
struct SideQR<GaussQR<Line,2>>
{
  using type = GaussQR<PointElem,1>;
};

template <>
struct SideQR<GaussQR<Line,3>>
{
  using type = GaussQR<PointElem,1>;
};

template <>
struct SideQR<GaussQR<Line,1>>
{
  using type = GaussQR<PointElem,1>;
};

template <>
struct SideQR<GaussQR<Triangle,1>>
{
  using type = GaussQR<Line,1>;
};

template <>
struct SideQR<GaussQR<Triangle,3>>
{
  using type = GaussQR<Line,2>;
};

template <>
struct SideQR<GaussQR<Triangle,4>>
{
  using type = GaussQR<Line,2>;
};

template <>
struct SideQR<GaussQR<Triangle,6>>
{
  using type = GaussQR<Line,3>;
};

template <>
struct SideQR<GaussQR<Triangle,7>>
{
  using type = GaussQR<Line,3>;
};

template <>
struct SideQR<GaussQR<Quad,1>>
{
  using type = GaussQR<Line,1>;
};

template <>
struct SideQR<GaussQR<Quad,4>>
{
  using type = GaussQR<Line,2>;
};

template <>
struct SideQR<GaussQR<Quad,9>>
{
  using type = GaussQR<Line,3>;
};

template <>
struct SideQR<TrapQR<Quad>>
{
  using type = TrapQR<Line>;
};

template <>
struct SideQR<GaussQR<Tetrahedron,4>>
{
  using type = GaussQR<Triangle,3>;
};

template <>
struct SideQR<GaussQR<Tetrahedron,11>>
{
  using type = GaussQR<Triangle,6>;
};

template <>
struct SideQR<GaussQR<Hexahedron,1>>
{
  using type = GaussQR<Quad,1>;
};

template <>
struct SideQR<GaussQR<Hexahedron,8>>
{
  using type = GaussQR<Quad,4>;
};

template <>
struct SideQR<GaussQR<Hexahedron,27>>
{
  using type = GaussQR<Quad,9>;
};

template <>
struct SideQR<TrapQR<Hexahedron>>
{
  using type = TrapQR<Quad>;
};
