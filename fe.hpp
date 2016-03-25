#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename Elem, uint order>
struct FEType {};

template <>
struct FEType<Line,1>
{
  typedef RefLineP1 RefFE_T;
  typedef GaussQR<Line,3> RecommendedQR;
};

template <>
struct FEType<Line,2>
{
  typedef RefLineP2 RefFE_T;
  typedef GaussQR<Line,3> RecommendedQR;
};

template <>
struct FEType<Triangle,1>
{
  typedef RefTriangleP1 RefFE_T;
  typedef GaussQR<Triangle,3> RecommendedQR;
};

template <>
struct FEType<Triangle,2>
{
  typedef RefTriangleP2 RefFE_T;
  typedef GaussQR<Triangle,4> RecommendedQR;
};

template <>
struct FEType<Quad,1>
{
  typedef RefQuadQ1 RefFE_T;
  typedef GaussQR<Quad,9> RecommendedQR;
};

template <>
struct FEType<Quad,2>
{
  typedef RefQuadQ2 RefFE_T;
  typedef GaussQR<Quad,9> RecommendedQR;
};
