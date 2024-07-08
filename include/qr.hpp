#pragma once

#include "def.hpp"

#include "geo.hpp"

namespace proxpde
{

// =====================================================================================
// TODO: define GaussQR by maximum order of the complete 1d polynomial, not number of
// points

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

// forward declarations

// 0D
template <>
GaussQR<NullElem, 0U>::Weights_T const GaussQR<NullElem, 0U>::weight;
template <>
std::array<GaussQR<NullElem, 0U>::Vec_T, 0U> const GaussQR<NullElem, 0U>::node;

// 1D
template <>
GaussQR<Line, 1U>::Weights_T const GaussQR<Line, 1U>::weight;
template <>
std::array<GaussQR<Line, 1>::Vec_T, 1U> const GaussQR<Line, 1>::node;
template <>
GaussQR<Line, 2U>::Weights_T const GaussQR<Line, 2U>::weight;
template <>
std::array<GaussQR<Line, 2U>::Vec_T, 2U> const GaussQR<Line, 2U>::node;
template <>
GaussQR<Line, 3U>::Weights_T const GaussQR<Line, 3U>::weight;
template <>
std::array<GaussQR<Line, 3U>::Vec_T, 3U> const GaussQR<Line, 3U>::node;
template <>
GaussQR<Line, 4U>::Weights_T const GaussQR<Line, 4U>::weight;
template <>
std::array<GaussQR<Line, 4>::Vec_T, 4U> const GaussQR<Line, 4U>::node;

// 2D
// - Triangle
template <>
GaussQR<Triangle, 1U>::Weights_T const GaussQR<Triangle, 1U>::weight;
template <>
std::array<GaussQR<Triangle, 1U>::Vec_T, 1U> const GaussQR<Triangle, 1U>::node;
template <>
GaussQR<Triangle, 3U>::Weights_T const GaussQR<Triangle, 3U>::weight;
template <>
std::array<GaussQR<Triangle, 3U>::Vec_T, 3U> const GaussQR<Triangle, 3U>::node;
template <>
GaussQR<Triangle, 4U>::Weights_T const GaussQR<Triangle, 4U>::weight;
template <>
std::array<GaussQR<Triangle, 4U>::Vec_T, 4U> const GaussQR<Triangle, 4U>::node;
template <>
GaussQR<Triangle, 6U>::Weights_T const GaussQR<Triangle, 6U>::weight;
template <>
std::array<GaussQR<Triangle, 6U>::Vec_T, 6U> const GaussQR<Triangle, 6U>::node;
template <>
GaussQR<Triangle, 7U>::Weights_T const GaussQR<Triangle, 7U>::weight;
template <>
std::array<GaussQR<Triangle, 7U>::Vec_T, 7U> const GaussQR<Triangle, 7U>::node;
// - Quad
template <>
GaussQR<Quad, 1U>::Weights_T const GaussQR<Quad, 1U>::weight;
template <>
std::array<GaussQR<Quad, 1U>::Vec_T, 1U> const GaussQR<Quad, 1U>::node;
template <>
GaussQR<Quad, 4U>::Weights_T const GaussQR<Quad, 4U>::weight;
template <>
std::array<GaussQR<Quad, 4U>::Vec_T, 4> const GaussQR<Quad, 4U>::node;
template <>
GaussQR<Quad, 9U>::Weights_T const GaussQR<Quad, 9U>::weight;
template <>
std::array<GaussQR<Quad, 9U>::Vec_T, 9> const GaussQR<Quad, 9U>::node;

// 3D
// - Tetrahedron
template <>
GaussQR<Tetrahedron, 1U>::Weights_T const GaussQR<Tetrahedron, 1U>::weight;
template <>
std::array<GaussQR<Tetrahedron, 1U>::Vec_T, 1U> const GaussQR<Tetrahedron, 1U>::node;
template <>
GaussQR<Tetrahedron, 4U>::Weights_T const GaussQR<Tetrahedron, 4U>::weight;
template <>
std::array<GaussQR<Tetrahedron, 4U>::Vec_T, 4U> const GaussQR<Tetrahedron, 4U>::node;
template <>
GaussQR<Tetrahedron, 11U>::Weights_T const GaussQR<Tetrahedron, 11U>::weight;
template <>
std::array<GaussQR<Tetrahedron, 11U>::Vec_T, 11U> const GaussQR<Tetrahedron, 11U>::node;
// - Hexahedron
template <>
GaussQR<Hexahedron, 1U>::Weights_T const GaussQR<Hexahedron, 1U>::weight;
template <>
std::array<GaussQR<Hexahedron, 1U>::Vec_T, 1U> const GaussQR<Hexahedron, 1U>::node;
template <>
GaussQR<Hexahedron, 8U>::Weights_T const GaussQR<Hexahedron, 8U>::weight;
template <>
std::array<GaussQR<Hexahedron, 8U>::Vec_T, 8U> const GaussQR<Hexahedron, 8U>::node;
template <>
GaussQR<Hexahedron, 27U>::Weights_T const GaussQR<Hexahedron, 27U>::weight;
template <>
std::array<GaussQR<Hexahedron, 27U>::Vec_T, 27U> const GaussQR<Hexahedron, 27U>::node;

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

// forward declarations

// 2D
// - Triangle
template <>
SideGaussQR<Triangle, 2U>::Weights_T const SideGaussQR<Triangle, 2U>::weight;
template <>
std::array<SideGaussQR<Triangle, 2U>::Vec_T, 2U * 3U> const
    SideGaussQR<Triangle, 2U>::node;
template <>
SideGaussQR<Triangle, 3U>::Weights_T const SideGaussQR<Triangle, 3U>::weight;
template <>
std::array<SideGaussQR<Triangle, 3U>::Vec_T, 3U * 3U> const
    SideGaussQR<Triangle, 3U>::node;
// - Quad
template <>
SideGaussQR<Quad, 2U>::Weights_T const SideGaussQR<Quad, 2U>::weight;
template <>
std::array<SideGaussQR<Quad, 2U>::Vec_T, 2U * 4U> const SideGaussQR<Quad, 2U>::node;
template <>
SideGaussQR<Quad, 3U>::Weights_T const SideGaussQR<Quad, 3U>::weight;
template <>
std::array<SideGaussQR<Quad, 3U>::Vec_T, 3U * 4U> const SideGaussQR<Quad, 3U>::node;

// 3D
// - Tetrahedron
template <>
SideGaussQR<Tetrahedron, 3U>::Weights_T const SideGaussQR<Tetrahedron, 3U>::weight;
template <>
std::array<SideGaussQR<Tetrahedron, 3U>::Vec_T, 3U * 4U> const
    SideGaussQR<Tetrahedron, 3U>::node;
template <>
SideGaussQR<Hexahedron, 4U>::Weights_T const SideGaussQR<Hexahedron, 4U>::weight;
template <>
std::array<SideGaussQR<Hexahedron, 4U>::Vec_T, 4U * 6U> const
    SideGaussQR<Hexahedron, 4U>::node;

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

// forward declarations

// 1D ----------------------------------------------------------------------------------
template <>
FVec<2U> const TrapQR<Line>::weight;
template <>
std::array<TrapQR<Line>::Vec_T, 2U> const TrapQR<Line>::node;

// 2D ----------------------------------------------------------------------------------
template <>
FVec<3U> const TrapQR<Triangle>::weight;
template <>
std::array<TrapQR<Triangle>::Vec_T, 3U> const TrapQR<Triangle>::node;

template <>
FVec<4U> const TrapQR<Quad>::weight;
template <>
std::array<TrapQR<Quad>::Vec_T, 4U> const TrapQR<Quad>::node;

// 3D ----------------------------------------------------------------------------------
template <>
FVec<4U> const TrapQR<Tetrahedron>::weight;
template <>
std::array<TrapQR<Tetrahedron>::Vec_T, 4U> const TrapQR<Tetrahedron>::node;

template <>
FVec<8U> const TrapQR<Hexahedron>::weight;
template <>
std::array<TrapQR<Hexahedron>::Vec_T, 8U> const TrapQR<Hexahedron>::node;

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

// forward declarations

// 1D ----------------------------------------------------------------------------------
template <>
FVec<3U> const SimpsonQR<Line>::weight;
template <>
std::array<SimpsonQR<Line>::Vec_T, 3U> const SimpsonQR<Line>::node;

// 2D ----------------------------------------------------------------------------------
template <>
FVec<9> const SimpsonQR<Quad>::weight;
template <>
std::array<SimpsonQR<Quad>::Vec_T, 9U> const SimpsonQR<Quad>::node;

// 3D ----------------------------------------------------------------------------------
template <>
FVec<27> const SimpsonQR<Hexahedron>::weight;
template <>
std::array<SimpsonQR<Hexahedron>::Vec_T, 27U> const SimpsonQR<Hexahedron>::node;

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
FVec<MiniQR<GeoElem, N>::numPts> const MiniQR<GeoElem, N>::weight =
    FVec<MiniQR<GeoElem, N>::numPts>::Constant(
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

template <typename GeoElem, short_T N>
struct DynamicQR
{
  using Real_T = double;
  using GeoElem_T = GeoElem;
  using Vec_T = FVec<GeoElem_T::dim>;
  using Weights_T = FVec<N>;
  // using Vec_T = Eigen::Matrix<long double, GeoElem::dim, 1>;
  // using Weights_T = Eigen::Matrix<long double, N, 1>;
  static short_T constexpr numPts = N;
  static short_T const bestPt;

  static Weights_T weight;
  static std::array<Vec_T, N> node;
};

template <typename GeoElem, short_T N>
typename DynamicQR<GeoElem, N>::Weights_T DynamicQR<GeoElem, N>::weight =
    FVec<N>::Constant(1.0 / N);

namespace details
{
template <typename GeoElem, short_T N>
static constexpr auto generateFixedNodes()
{
  auto data = std::array<FVec<GeoElem::dim>, N>{};
  for (uint n = 0; n < N; ++n)
  {
    data[n] = FVec<GeoElem::dim>::Constant(0.0);
  }
  return data;
}
} // namespace details

template <typename GeoElem, short_T N>
std::array<typename DynamicQR<GeoElem, N>::Vec_T, N> DynamicQR<GeoElem, N>::node =
    details::generateFixedNodes<GeoElem, N>();

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

} // namespace proxpde
