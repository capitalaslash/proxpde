#pragma once

#include "def.hpp"

#include "geo.hpp"
#include "reffe.hpp"

// =====================================================================================

template <typename GeoElem, short_T N>
struct GaussQR
{
  using Real_T = double;
  using GeoElem_T = GeoElem;
  using Vec_T = FVec<GeoElem_T::dim>;
  using Weights_T = FVec<N>;
  // using Vec_T = Eigen::Matrix<long double, GeoElem::dim, 1>;
  // using Weights_T = Eigen::Matrix<long double, N, 1>;
  static short_T constexpr numPts = N;
  static short_T const bestPt;

  static Weights_T const weight;
  static std::array<Vec_T, N> const node;
};

// 0D ----------------------------------------------------------------------------------
template <>
GaussQR<NullElem, 0>::Weights_T const GaussQR<NullElem, 0>::weight;
template <>
std::array<GaussQR<NullElem, 0>::Vec_T, 0> const GaussQR<NullElem, 0>::node;
template <>
short_T const GaussQR<NullElem, 0>::bestPt;

// 1D ----------------------------------------------------------------------------------
template <>
GaussQR<Line, 1>::Weights_T const GaussQR<Line, 1>::weight;
template <>
std::array<GaussQR<Line, 1>::Vec_T, 1> const GaussQR<Line, 1>::node;
template <>
short_T const GaussQR<Line, 1>::bestPt;

template <>
GaussQR<Line, 2>::Weights_T const GaussQR<Line, 2>::weight;
template <>
std::array<GaussQR<Line, 2>::Vec_T, 2> const GaussQR<Line, 2>::node;
template <>
short_T const GaussQR<Line, 2>::bestPt;

template <>
GaussQR<Line, 3>::Weights_T const GaussQR<Line, 3>::weight;
template <>
std::array<GaussQR<Line, 3>::Vec_T, 3> const GaussQR<Line, 3>::node;
template <>
short_T const GaussQR<Line, 3>::bestPt;

template <>
GaussQR<Line, 4>::Weights_T const GaussQR<Line, 4>::weight;
template <>
std::array<GaussQR<Line, 4>::Vec_T, 4> const GaussQR<Line, 4>::node;
template <>
short_T const GaussQR<Line, 4>::bestPt;

// 2D Triangle -------------------------------------------------------------------------
template <>
GaussQR<Triangle, 1>::Weights_T const GaussQR<Triangle, 1>::weight;
template <>
std::array<GaussQR<Triangle, 1>::Vec_T, 1> const GaussQR<Triangle, 1>::node;
template <>
short_T const GaussQR<Triangle, 1>::bestPt;

template <>
GaussQR<Triangle, 3>::Weights_T const GaussQR<Triangle, 3>::weight;
template <>
std::array<GaussQR<Triangle, 3>::Vec_T, 3> const GaussQR<Triangle, 3>::node;
template <>
short_T const GaussQR<Triangle, 3>::bestPt;

template <>
GaussQR<Triangle, 6>::Weights_T const GaussQR<Triangle, 6>::weight;
template <>
std::array<GaussQR<Triangle, 6>::Vec_T, 6> const GaussQR<Triangle, 6>::node;
template <>
short_T const GaussQR<Triangle, 6>::bestPt;

template <>
GaussQR<Triangle, 7>::Weights_T const GaussQR<Triangle, 7>::weight;
template <>
std::array<GaussQR<Triangle, 7>::Vec_T, 7> const GaussQR<Triangle, 7>::node;
template <>
short_T const GaussQR<Triangle, 7>::bestPt;

// 2D Quad -----------------------------------------------------------------------------
template <>
GaussQR<Quad, 1>::Weights_T const GaussQR<Quad, 1>::weight;
template <>
std::array<GaussQR<Quad, 1>::Vec_T, 1> const GaussQR<Quad, 1>::node;
template <>
short_T const GaussQR<Quad, 1>::bestPt;

template <>
GaussQR<Quad, 4>::Weights_T const GaussQR<Quad, 4>::weight;
template <>
std::array<GaussQR<Quad, 4>::Vec_T, 4> const GaussQR<Quad, 4>::node;
template <>
short_T const GaussQR<Quad, 4>::bestPt;

template <>
GaussQR<Quad, 9>::Weights_T const GaussQR<Quad, 9>::weight;
template <>
std::array<GaussQR<Quad, 9>::Vec_T, 9> const GaussQR<Quad, 9>::node;
template <>
short_T const GaussQR<Quad, 9>::bestPt;

// 3D Tetra ----------------------------------------------------------------------------
template <>
GaussQR<Tetrahedron, 1>::Weights_T const GaussQR<Tetrahedron, 1>::weight;
template <>
std::array<GaussQR<Tetrahedron, 1>::Vec_T, 1> const GaussQR<Tetrahedron, 1>::node;
template <>
short_T const GaussQR<Tetrahedron, 1>::bestPt;

template <>
GaussQR<Tetrahedron, 4>::Weights_T const GaussQR<Tetrahedron, 4>::weight;
template <>
std::array<GaussQR<Tetrahedron, 4>::Vec_T, 4> const GaussQR<Tetrahedron, 4>::node;
template <>
short_T const GaussQR<Tetrahedron, 4>::bestPt;

template <>
GaussQR<Tetrahedron, 11>::Weights_T const GaussQR<Tetrahedron, 11>::weight;
template <>
std::array<GaussQR<Tetrahedron, 11>::Vec_T, 11> const GaussQR<Tetrahedron, 11>::node;
template <>
short_T const GaussQR<Tetrahedron, 11>::bestPt;

// 3D Hexahedron -----------------------------------------------------------------------
template <>
GaussQR<Hexahedron, 1>::Weights_T const GaussQR<Hexahedron, 1>::weight;
template <>
std::array<GaussQR<Hexahedron, 1>::Vec_T, 1> const GaussQR<Hexahedron, 1>::node;
template <>
short_T const GaussQR<Hexahedron, 1>::bestPt;

template <>
GaussQR<Hexahedron, 8>::Weights_T const GaussQR<Hexahedron, 8>::weight;
template <>
std::array<GaussQR<Hexahedron, 8>::Vec_T, 8> const GaussQR<Hexahedron, 8>::node;
template <>
short_T const GaussQR<Hexahedron, 8>::bestPt;

template <>
GaussQR<Hexahedron, 27>::Weights_T const GaussQR<Hexahedron, 27>::weight;
template <>
std::array<GaussQR<Hexahedron, 27>::Vec_T, 27> const GaussQR<Hexahedron, 27>::node;
template <>
short_T const GaussQR<Hexahedron, 27>::bestPt;

// =====================================================================================

template <typename GeoElem, short_T N>
struct SideGaussQR
{
  using GeoElem_T = GeoElem;
  using Vec_T = FVec<GeoElem_T::dim>;
  static short_T constexpr numPts = N * GeoElem_T::numFacets;
  using Weights_T = FVec<numPts>;

  static Weights_T const weight;
  static std::array<Vec_T, numPts> const node;
};

// 2D Triangle -------------------------------------------------------------------------
template <>
SideGaussQR<Triangle, 2>::Weights_T const SideGaussQR<Triangle, 2>::weight;
template <>
std::array<SideGaussQR<Triangle, 2>::Vec_T, 2 * 3> const SideGaussQR<Triangle, 2>::node;

// 2D Quad -------------------------------------------------------------------------
template <>
SideGaussQR<Quad, 2>::Weights_T const SideGaussQR<Quad, 2>::weight;
template <>
std::array<SideGaussQR<Quad, 2>::Vec_T, 2 * 4> const SideGaussQR<Quad, 2>::node;

template <>
SideGaussQR<Quad, 3>::Weights_T const SideGaussQR<Quad, 3>::weight;
template <>
std::array<SideGaussQR<Quad, 3>::Vec_T, 3 * 4> const SideGaussQR<Quad, 3>::node;

// =====================================================================================

template <typename GeoElem>
struct TrapQR
{
  using GeoElem_T = GeoElem;
  using Vec_T = FVec<GeoElem::dim>;
  static short_T const numPts = GeoElem::numPts;

  static FVec<numPts> const weight;
  static std::array<Vec_T, numPts> const node;
};

// 1D ----------------------------------------------------------------------------------
template <>
FVec<2> const TrapQR<Line>::weight;
template <>
std::array<TrapQR<Line>::Vec_T, 2> const TrapQR<Line>::node;

// 2D ----------------------------------------------------------------------------------
template <>
FVec<3> const TrapQR<Triangle>::weight;
template <>
std::array<TrapQR<Triangle>::Vec_T, 3> const TrapQR<Triangle>::node;

template <>
FVec<4> const TrapQR<Quad>::weight;
template <>
std::array<TrapQR<Quad>::Vec_T, 4> const TrapQR<Quad>::node;

// 3D ----------------------------------------------------------------------------------
template <>
FVec<8> const TrapQR<Hexahedron>::weight;
template <>
std::array<TrapQR<Hexahedron>::Vec_T, 8> const TrapQR<Hexahedron>::node;

// =====================================================================================

template <typename GeoElem>
struct SimpsonQR
{
  using GeoElem_T = GeoElem;
  using Vec_T = FVec<GeoElem::dim>;
  static short_T const numPts = cepow(3, GeoElem::dim);

  static FVec<numPts> const weight;
  static std::array<Vec_T, numPts> const node;
};

// 1D ----------------------------------------------------------------------------------
template <>
FVec<3> const SimpsonQR<Line>::weight;
template <>
std::array<SimpsonQR<Line>::Vec_T, 3> const SimpsonQR<Line>::node;

// 2D ----------------------------------------------------------------------------------
template <>
FVec<9> const SimpsonQR<Quad>::weight;
template <>
std::array<SimpsonQR<Quad>::Vec_T, 9> const SimpsonQR<Quad>::node;

// 3D ----------------------------------------------------------------------------------
template <>
FVec<27> const SimpsonQR<Hexahedron>::weight;
template <>
std::array<SimpsonQR<Hexahedron>::Vec_T, 27> const SimpsonQR<Hexahedron>::node;

// =====================================================================================

template <typename GeoElem, uint N>
struct MiniQR
{
  using GeoElem_T = GeoElem;
  using Vec_T = FVec<GeoElem::dim>;
  static uint constexpr numPts = cepow(N, GeoElem::dim);

  static FVec<numPts> const weight;
  static std::array<Vec_T, numPts> const node;
};

template <typename GeoElem, uint N>
FVec<MiniQR<GeoElem, N>::numPts> const
    MiniQR<GeoElem, N>::weight = FVec<MiniQR<GeoElem, N>::numPts>::Constant(
        GeoElem::refVolume / MiniQR<GeoElem, N>::numPts);

template <typename GeoElem, uint N>
static std::array<typename MiniQR<GeoElem, N>::Vec_T, MiniQR<GeoElem, N>::numPts>
miniNodes()
{
  std::array<typename MiniQR<GeoElem, N>::Vec_T, MiniQR<GeoElem, N>::numPts> pts;

  if constexpr (std::is_same_v<GeoElem, Line>)
  {
    auto const h = 2. / N;
    for (short_T k = 0; k < N; ++k)
    {
      pts[k] = FVec<1>::Constant(-1. + (k + .5) * h);
    }
  }
  else if constexpr (std::is_same_v<GeoElem, Triangle>)
  {
    // k is the number of stripes of triangles
    for (uint k = 0; k < N; ++k)
    {
      auto const f = 2 * (k + 1) - 1;
      for (uint i = 0; i < f; ++i)
      {
        pts[k * k + i] = FVec<2>{
            static_cast<double>(i + 1 + i / 2) / (3 * N),
            static_cast<double>(f - i + (f - i - 1) / 2) / (3 * N)};
      }
    }
  }
  else if constexpr (std::is_same_v<GeoElem, Quad>)
  {
    auto const h = 2. / N;
    for (uint i = 0; i < N; ++i)
    {
      for (uint j = 0; j < N; ++j)
      {
        pts[i * N + j] = FVec<2>(-1. + (i + .5) * h, -1. + (j + .5) * h);
      }
    }
  }
  else
  {
    // we should never reach this point
    std::abort();
  }
  return pts;
}

template <typename GeoElem, uint N>
std::array<typename MiniQR<GeoElem, N>::Vec_T, MiniQR<GeoElem, N>::numPts> const
    MiniQR<GeoElem, N>::node = miniNodes<GeoElem, N>();

// =====================================================================================

template <typename QR>
struct SideQR
{};

template <typename QR>
using SideQR_T = typename SideQR<QR>::type;

template <>
struct SideQR<GaussQR<Line, 2>>
{
  using type = GaussQR<PointElem, 1>;
};

template <>
struct SideQR<GaussQR<Line, 3>>
{
  using type = GaussQR<PointElem, 1>;
};

template <>
struct SideQR<GaussQR<Line, 1>>
{
  using type = GaussQR<PointElem, 1>;
};

template <>
struct SideQR<GaussQR<Triangle, 1>>
{
  using type = GaussQR<Line, 1>;
};

template <>
struct SideQR<GaussQR<Triangle, 3>>
{
  using type = GaussQR<Line, 2>;
};

template <>
struct SideQR<GaussQR<Triangle, 4>>
{
  using type = GaussQR<Line, 2>;
};

template <>
struct SideQR<GaussQR<Triangle, 6>>
{
  using type = GaussQR<Line, 3>;
};

template <>
struct SideQR<GaussQR<Triangle, 7>>
{
  using type = GaussQR<Line, 3>;
};

template <>
struct SideQR<GaussQR<Quad, 1>>
{
  using type = GaussQR<Line, 1>;
};

template <>
struct SideQR<GaussQR<Quad, 4>>
{
  using type = GaussQR<Line, 2>;
};

template <>
struct SideQR<GaussQR<Quad, 9>>
{
  using type = GaussQR<Line, 3>;
};

template <>
struct SideQR<TrapQR<Quad>>
{
  using type = TrapQR<Line>;
};

template <>
struct SideQR<GaussQR<Tetrahedron, 1>>
{
  using type = GaussQR<Triangle, 1>;
};

template <>
struct SideQR<GaussQR<Tetrahedron, 4>>
{
  using type = GaussQR<Triangle, 3>;
};

template <>
struct SideQR<GaussQR<Tetrahedron, 11>>
{
  using type = GaussQR<Triangle, 6>;
};

template <>
struct SideQR<GaussQR<Hexahedron, 1>>
{
  using type = GaussQR<Quad, 1>;
};

template <>
struct SideQR<GaussQR<Hexahedron, 8>>
{
  using type = GaussQR<Quad, 4>;
};

template <>
struct SideQR<GaussQR<Hexahedron, 27>>
{
  using type = GaussQR<Quad, 9>;
};

template <>
struct SideQR<TrapQR<Hexahedron>>
{
  using type = TrapQR<Quad>;
};
