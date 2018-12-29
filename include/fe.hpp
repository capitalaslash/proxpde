#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename Elem, uint order>
struct FEType {};

template <>
struct FEType<Line,0>
{
  using RefFE_T = RefLineP0;
  using RecommendedQR = GaussQR<Line,1>;
};

template <>
struct FEType<Line,1>
{
  using RefFE_T = RefLineP1;
  using RecommendedQR = GaussQR<Line,3>;
};

template <>
struct FEType<Line,2>
{
  using RefFE_T = RefLineP2;
  using RecommendedQR = GaussQR<Line,3>;
};

template <>
struct FEType<Triangle,1>
{
  using RefFE_T = RefTriangleP1;
  using RecommendedQR = GaussQR<Triangle,3>;
};

template <>
struct FEType<Triangle,2>
{
  using RefFE_T = RefTriangleP2;
  using RecommendedQR = GaussQR<Triangle,7>;
};

template <>
struct FEType<Quad,1>
{
  using RefFE_T = RefQuadQ1;
  using RecommendedQR = GaussQR<Quad,9>;
};

template <>
struct FEType<Quad,2>
{
  using RefFE_T = RefQuadQ2;
  using RecommendedQR = GaussQR<Quad,9>;
};

template <>
struct FEType<Tetrahedron,1>
{
  using RefFE_T = RefTetrahedronP1;
  using RecommendedQR = GaussQR<Tetrahedron,4>;
};
