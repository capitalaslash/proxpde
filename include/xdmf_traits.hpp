#pragma once

#include "def.hpp"
#include "reffe.hpp"

enum class XDMFGridType : int8_t
{
  SINGLE,
  COLLECTION
};

static const std::unordered_map<XDMFGridType,std::string> XDMFGridTypeToString =
{
  {XDMFGridType::SINGLE, "Uniform"},
  {XDMFGridType::COLLECTION, "Collection"},
};

enum class XDMFNumberType : int8_t
{
  INT,
  FLOAT
};

static const std::unordered_map<XDMFNumberType, char const *> XDMFNumberTypeToString =
{
  {XDMFNumberType::INT, "Int"},
  {XDMFNumberType::FLOAT, "Float"},
};

enum class XDMFFormat : int8_t
{
  HDF,
  // INLINE
};

static const std::unordered_map<XDMFFormat, char const *> XDMFFormatToString =
{
  {XDMFFormat::HDF, "HDF"},
  // {XDMFFormat::INLINE, "XML"},
};

enum class XDMFCenter : int8_t
{
  CELL,
  NODE
};

static const std::unordered_map<XDMFCenter, std::string> XDMFCenterToString =
{
  {XDMFCenter::CELL, "Cell"},
  {XDMFCenter::NODE, "Node"},
};

struct XDMFVar
{
  std::string const name;
  XDMFCenter const center;
  ulong const size;
};

template <typename RefFE>
struct XDMFTraits {};

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
  static constexpr char const * shapeName = "Polyline";
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
struct XDMFTraits<RefTetrahedronP2>
{
  static constexpr char const * shapeName = "Tetrahedron_10";
  static constexpr XDMFCenter attributeType = XDMFCenter::NODE;
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
  static constexpr array<char,27> mapping = {
    0, 1, 2, 3, 4, 5, 6, 7,
    8, 9, 10, 11,
    16, 17, 18, 19,
    12, 13, 14, 15,
    24, 22, 21, 23, 20, 25,
    26};
};
