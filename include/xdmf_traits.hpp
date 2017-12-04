#pragma once

#include "def.hpp"
#include "reffe.hpp"

template <typename RefFE>
struct XDMFTraits
{
  static constexpr char const * shape_name = "None";
};

template <>
struct XDMFTraits<RefLineP1>
{
  static constexpr char const * shape_name = "Polyline";
};

template <>
struct XDMFTraits<RefTriangleP1>
{
  static constexpr char const * shape_name = "Triangle";
};

template <>
struct XDMFTraits<RefTriangleP2>
{
  static constexpr char const * shape_name = "Triangle_6";
};

template <>
struct XDMFTraits<RefQuadQ1>
{
  static constexpr char const * shape_name = "Quadrilateral";
};

template <>
struct XDMFTraits<RefQuadP2>
{
  static constexpr char const * shape_name = "Quadrilateral_8";
};

template <>
struct XDMFTraits<RefQuadQ2>
{
  static constexpr char const * shape_name = "Quadrilateral_9";
};
