#include "geo.hpp"

std::array<std::array<id_T,1>,2> constexpr Line::elemToFacet;
std::array<std::array<id_T,0>,0> constexpr Line::elemToFace;
std::array<std::array<id_T,2>,1> constexpr Line::elemToEdge;

std::array<std::array<id_T,2>,3> constexpr Triangle::elemToFacet;
std::array<std::array<id_T,3>,1> constexpr Triangle::elemToFace;
std::array<std::array<id_T,2>,3> constexpr Triangle::elemToEdge;

std::array<std::array<id_T,2>,4> constexpr Quad::elemToFacet;
std::array<std::array<id_T,4>,1> constexpr Quad::elemToFace;
std::array<std::array<id_T,2>,4> constexpr Quad::elemToEdge;
