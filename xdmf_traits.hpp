#pragma once

#include "def.hpp"
#include "geo.hpp"

template <typename Elem>
struct XDMFTraits
{
  static constexpr char const * shape_name = "None";
};

template <>
struct XDMFTraits<Line>
{
  static constexpr char const * shape_name = "Polyline";
};

template <>
struct XDMFTraits<Triangle>
{
  static constexpr char const * shape_name = "Triangle";
};

template <>
struct XDMFTraits<Quad>
{
  static constexpr char const * shape_name = "Quadrilateral";
};
