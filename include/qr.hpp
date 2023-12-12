#pragma once

#include "def.hpp"

#include "geo.hpp"

namespace proxpde
{

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
DynamicQR<GeoElem, N>::Weights_T DynamicQR<GeoElem, N>::weight =
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
