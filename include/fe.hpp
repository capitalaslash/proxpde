#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename Elem, uint order>
struct LagrangeFE {};

template <>
struct LagrangeFE<PointElem,0>
{
  using RefFE_T = RefPointP1;
  using RecommendedQR = GaussQR<PointElem,1>;
};

template <>
struct LagrangeFE<Line,0>
{
  using RefFE_T = RefLineP0;
  using RecommendedQR = GaussQR<Line,1>;
};

template <>
struct LagrangeFE<Line,1>
{
  using RefFE_T = RefLineP1;
  using RecommendedQR = GaussQR<Line,2>;
  using ReconstructionQR = TrapQR<Line>;
};

template <>
struct LagrangeFE<Line,2>
{
  using RefFE_T = RefLineP2;
  using RecommendedQR = GaussQR<Line,3>;
  using ReconstructionQR = SimpsonQR<Line>;
};

template <>
struct LagrangeFE<Triangle,0>
{
  using RefFE_T = RefTriangleP0;
  using RecommendedQR = GaussQR<Triangle,1>;
};

template <>
struct LagrangeFE<Triangle,1>
{
  using RefFE_T = RefTriangleP1;
  using RecommendedQR = GaussQR<Triangle,3>;
};

template <>
struct LagrangeFE<Triangle,2>
{
  using RefFE_T = RefTriangleP2;
  using RecommendedQR = GaussQR<Triangle,6>;
};

template <>
struct LagrangeFE<Quad,0>
{
  using RefFE_T = RefQuadP0;
  using RecommendedQR = GaussQR<Quad,1>;
};

template <>
struct LagrangeFE<Quad,1>
{
  using RefFE_T = RefQuadQ1;
  using RecommendedQR = GaussQR<Quad,4>;
  using ReconstructionQR = TrapQR<Quad>;
};

template <>
struct LagrangeFE<Quad,2>
{
  using RefFE_T = RefQuadQ2;
  using RecommendedQR = GaussQR<Quad,9>;
  using ReconstructionQR = SimpsonQR<Quad>;
};

template <>
struct LagrangeFE<Tetrahedron,0>
{
  using RefFE_T = RefTetrahedronP0;
  using RecommendedQR = GaussQR<Tetrahedron,1>;
};

template <>
struct LagrangeFE<Tetrahedron,1>
{
  using RefFE_T = RefTetrahedronP1;
  using RecommendedQR = GaussQR<Tetrahedron,4>;
};

template <>
struct LagrangeFE<Tetrahedron,2>
{
  using RefFE_T = RefTetrahedronP2;
  using RecommendedQR = GaussQR<Tetrahedron,11>;
};

template <>
struct LagrangeFE<Hexahedron,0>
{
  using RefFE_T = RefHexahedronP0;
  using RecommendedQR = GaussQR<Hexahedron,1>;
};

template <>
struct LagrangeFE<Hexahedron,1>
{
  using RefFE_T = RefHexahedronQ1;
  using RecommendedQR = GaussQR<Hexahedron,8>;
  using ReconstructionQR = TrapQR<Hexahedron>;
};

template <>
struct LagrangeFE<Hexahedron,2>
{
  using RefFE_T = RefHexahedronQ2;
  using RecommendedQR = GaussQR<Hexahedron,27>;
  using ReconstructionQR = SimpsonQR<Hexahedron>;
};
