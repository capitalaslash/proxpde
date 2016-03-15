#pragma once

#include "def.hpp"
#include "reffe.hpp"

template <typename Elem, uint order>
struct FEType {};

template <>
struct FEType<Line,1>
{
  typedef RefLineP1 RefFE_T;
};

template <>
struct FEType<Line,2>
{
  typedef RefLineP2 RefFE_T;
};

template <>
struct FEType<Triangle,1>
{
  typedef RefTriangleP1 RefFE_T;
};

template <>
struct FEType<Quad,1>
{
  typedef RefQuadQ1 RefFE_T;
};
