#pragma once

#include "def.hpp"
#include "reffe.hpp"

template <typename RefFE>
struct XDMFTraits {};

template <>
struct XDMFTraits<RefLineP0>
{
  static constexpr char const * shapeName = "Polyline";
  static constexpr char const * attributeType = "Cell";
};

template <>
struct XDMFTraits<RefLineP1>
{
  static constexpr char const * shapeName = "Polyline";
  static constexpr char const * attributeType = "Node";
};

template <>
struct XDMFTraits<RefLineP2>
{
  static constexpr char const * shapeName = "Polyline";
  static constexpr char const * attributeType = "Node";
};

template <>
struct XDMFTraits<RefTriangleP0>
{
  static constexpr char const * shapeName = "Triangle";
  static constexpr char const * attributeType = "Cell";
};

template <>
struct XDMFTraits<RefTriangleP1>
{
  static constexpr char const * shapeName = "Triangle";
  static constexpr char const * attributeType = "Node";
};

template <>
struct XDMFTraits<RefTriangleP2>
{
  static constexpr char const * shapeName = "Triangle_6";
  static constexpr char const * attributeType = "Node";
};

template <>
struct XDMFTraits<RefQuadQ1>
{
  static constexpr char const * shapeName = "Quadrilateral";
  static constexpr char const * attributeType = "Node";
};

template <>
struct XDMFTraits<RefQuadP2>
{
  static constexpr char const * shapeName = "Quadrilateral_8";
  static constexpr char const * attributeType = "Node";
};

template <>
struct XDMFTraits<RefQuadQ2>
{
  static constexpr char const * shapeName = "Quadrilateral_9";
  static constexpr char const * attributeType = "Node";
};

template <>
struct XDMFTraits<RefTetrahedronP1>
{
  static constexpr char const * shapeName = "Tetrahedron";
  static constexpr char const * attributeType = "Node";
};

template <>
struct XDMFTraits<RefHexahedronQ1>
{
  static constexpr char const * shapeName = "Hexahedron";
  static constexpr char const * attributeType = "Node";
};
