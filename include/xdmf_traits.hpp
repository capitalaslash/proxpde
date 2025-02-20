#pragma once

#include "def.hpp"

#include "reffe.hpp"

namespace proxpde
{

enum class XDMFGridType : uint8_t
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
  default:
    abort();
  }
}

enum class XDMFNumberType : uint8_t
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
  default:
    abort();
  }
}

enum class XDMFFormat : uint8_t
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
  default:
    abort();
  }
}

enum class XDMFCenter : uint8_t
{
  CELL,
  NODE,
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
  default:
    abort();
  }
}

enum class XDMFType : uint8_t
{
  SCALAR,
  VECTOR,
};

constexpr std::string_view to_string(XDMFType const type)
{
  switch (type)
  {
  case XDMFType::SCALAR:
    return "Scalar";
  case XDMFType::VECTOR:
    return "Vector";
  default:
    abort();
  }
}

struct XDMFVar
{
  std::string const name;
  XDMFType const type;
  XDMFCenter const center;
  XDMFNumberType const numberType;
  ulong const size;
  uint const dim;
};

template <typename RefFE>
struct XDMFTraits
{};

template <>
struct XDMFTraits<RefLineP0>
{
  static constexpr char const * shapeName = "Polyline";
  static constexpr XDMFCenter center = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefLineP1>
{
  static constexpr char const * shapeName = "Polyline";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefLineP2>
{
  static constexpr char const * shapeName = "Edge_3";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTriangleP0>
{
  static constexpr char const * shapeName = "Triangle";
  static constexpr XDMFCenter center = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTriangleP1>
{
  static constexpr char const * shapeName = "Triangle";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTriangleP2>
{
  static constexpr char const * shapeName = "Triangle_6";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefQuadP0>
{
  static constexpr char const * shapeName = "Quadrilateral";
  static constexpr XDMFCenter center = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefQuadQ1>
{
  static constexpr char const * shapeName = "Quadrilateral";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefQuadP2>
{
  static constexpr char const * shapeName = "Quadrilateral_8";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefQuadQ2>
{
  static constexpr char const * shapeName = "Quadrilateral_9";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTetrahedronP1>
{
  static constexpr char const * shapeName = "Tetrahedron";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTetrahedronP0>
{
  static constexpr char const * shapeName = "Tetrahedron";
  static constexpr XDMFCenter center = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefTetrahedronP2>
{
  static constexpr char const * shapeName = "Tetrahedron_10";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefHexahedronP0>
{
  static constexpr char const * shapeName = "Hexahedron";
  static constexpr XDMFCenter center = XDMFCenter::CELL;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefHexahedronQ1>
{
  static constexpr char const * shapeName = "Hexahedron";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
  static constexpr bool needsMapping = false;
};

template <>
struct XDMFTraits<RefHexahedronQ2>
{
  static constexpr char const * shapeName = "Hexahedron_27";
  static constexpr XDMFCenter center = XDMFCenter::NODE;
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
  // order: bottom corners, top corners, bottom middles, top middles, middle corners,
  // face centers (x-, x+, y-, y+, z-, z+), middle
  // clang-format off
  static constexpr std::array<short_T, 27> mapping = {
    // bottom corners
    0, 1, 2, 3,
    // top corners
    4, 5, 6, 7,
    // bottom middles
    8, 9, 10, 11,
    // middle corners
    16, 17, 18, 19,
    // top middles
    12, 13, 14, 15,
    // face centers (z-, x-, y-, x+, y+, z+)
    21, 23, 22, 24, 20, 25,
    // middle
    26,
  };
  // clang-format on
};

} // namespace proxpde
