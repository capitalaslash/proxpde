#include "geo.hpp"

array<array<id_T,1>,2> constexpr Line::elemToFacet;
array<array<id_T,0>,0> constexpr Line::elemToFace;
array<array<id_T,2>,1> constexpr Line::elemToEdge;

array<array<id_T,2>,3> constexpr Triangle::elemToFacet;
array<array<id_T,3>,1> constexpr Triangle::elemToFace;
array<array<id_T,2>,3> constexpr Triangle::elemToEdge;

array<array<id_T,2>,4> constexpr Quad::elemToFacet;
array<array<id_T,4>,1> constexpr Quad::elemToFace;
array<array<id_T,2>,4> constexpr Quad::elemToEdge;
