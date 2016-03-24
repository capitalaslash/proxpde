#include "geo.hpp"

std::array<std::array<id_T,1>,2> constexpr Line::elemToFacet;

std::array<std::array<id_T,2>,3> constexpr Triangle::elemToFacet;

std::array<std::array<id_T,2>,4> constexpr Quad::elemToFacet;
