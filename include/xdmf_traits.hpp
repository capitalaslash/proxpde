#pragma once

#include "def.hpp"

#include "reffe.hpp"

enum class XDMFGridType : int8_t
{
  SINGLE,
  COLLECTION
};

constexpr std::string_view to_string(XDMFGridType const type)
{
  switch (type)
  {
    // using enum XDMFGridType; c++20
  case XDMFGridType::SINGLE:
    return "Uniform";
  case XDMFGridType::COLLECTION:
    return "Collection";
  }
  return "ERROR";
}

enum class XDMFNumberType : int8_t
{
  INT,
  FLOAT
};

constexpr std::string_view to_string(XDMFNumberType const type)
{
  switch (type)
  {
    // using enum XDMFNumberType; c++20
  case XDMFNumberType::INT:
    return "Int";
  case XDMFNumberType::FLOAT:
    return "Float";
  }
  return "ERROR";
}

enum class XDMFFormat : int8_t
{
  HDF,
  INLINE
};

constexpr std::string_view to_string(XDMFFormat const format)
{
  switch (format)
  {
    // using enum XDMFFormat; c++20
  case XDMFFormat::HDF:
    return "HDF";
  case XDMFFormat::INLINE:
    return "XML";
  }
  return "ERROR";
}

enum class XDMFCenter : int8_t
{
  CELL,
  NODE
};

constexpr std::string_view to_string(XDMFCenter const center)
{
  switch (center)
  {
    // using enum XDMFCenter; c++20
  case XDMFCenter::CELL:
    return "Cell";
  case XDMFCenter::NODE:
    return "Node";
  }
  return "ERROR";
}

struct XDMFVar
{
  std::string const name;
  XDMFCenter const center;
  ulong const size;
};

template <typename RefFE>
struct XDMFTraits
{};

template <>
struct XDMFTraits<RefLineP0>
{
  static constexpr char const * shapeName = "Polyline";
  static constexpr XDMFCenter attributeType = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefLineP1>
{
  static constexpr char const * shapeName = "Polyline";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefLineP2>
{
  static constexpr char const * shapeName = "Edge_3";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTriangleP0>
{
  static constexpr char const * shapeName = "Triangle";
  static constexpr XDMFCenter attributeType = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTriangleP1>
{
  static constexpr char const * shapeName = "Triangle";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTriangleP2>
{
  static constexpr char const * shapeName = "Triangle_6";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefQuadP0>
{
  static constexpr char const * shapeName = "Quadrilateral";
  static constexpr XDMFCenter attributeType = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefQuadQ1>
{
  static constexpr char const * shapeName = "Quadrilateral";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefQuadP2>
{
  static constexpr char const * shapeName = "Quadrilateral_8";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefQuadQ2>
{
  static constexpr char const * shapeName = "Quadrilateral_9";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTetrahedronP1>
{
  static constexpr char const * shapeName = "Tetrahedron";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTetrahedronP0>
{
  static constexpr char const * shapeName = "Tetrahedron";
  static constexpr XDMFCenter attributeType = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTetrahedronP2>
{
  static constexpr char const * shapeName = "Tetrahedron_10";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefHexahedronP0>
{
  static constexpr char const * shapeName = "Hexahedron";
  static constexpr XDMFCenter attributeType = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefHexahedronQ1>
{
  static constexpr char const * shapeName = "Hexahedron";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefHexahedronQ2>
{
  static constexpr char const * shapeName = "Hexahedron_27";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
  static constexpr bool needsMapping = true;
  // from VTK_TRIQUADRATIC_HEXAHEDRON documentation
  // top
  //  7--14--6
  //  |      |
  // 15  25  13
  //  |      |
  //  4--12--5
  //
  //  middle
  // 19--23--18
  //  |      |
  // 20  26  21
  //  |      |
  // 16--22--17
  //
  // bottom
  //  3--10--2
  //  |      |
  // 11  24  9
  //  |      |
  //  0-- 8--1
  // clang-format off
  static constexpr std::array<short_T, 27> mapping = {
    0, 1, 2, 3, 4, 5, 6, 7,
    8, 9, 10, 11,
    16, 17, 18, 19,
    12, 13, 14, 15,
    24, 22, 21, 23, 20, 25,
    26};
  // clang-format on
};
