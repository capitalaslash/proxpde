#pragma once

#include "def.hpp"

#include "geo.hpp"

template <typename RefElem>
uint constexpr numDOFs()
{
  return 1 * RefElem::dofPlace[0] +
         RefElem::GeoElem_T::numFaces * RefElem::dofPlace[1] +
         RefElem::GeoElem_T::numEdges * RefElem::dofPlace[2] +
         RefElem::GeoElem_T::numPts * RefElem::dofPlace[3];
}

struct RefPointP1
{
  using GeoElem_T = PointElem;
  using FacetFE_T = NullElem;
  static GeoElem_T const geoElem;
  static int constexpr dim = 1;
  static uint constexpr numDOFs = 1U;
  static uint constexpr numGeoDOFs = 1U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 1, 1> constexpr dofOnFacet = {{{{0}}}};
  using Vec_T = FVec<dim>;

  static std::array<Vec_T, numDOFs> const points;
  static std::array<scalarOnedFun_T, numDOFs> const phiFun;
  static std::array<onedFun_T, numDOFs> const dphiFun;
  static std::array<onedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 1.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts = {{e.pointList[0]->coord}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts = {{e.pointList[0]->coord}};
    return mappingPts;
  }

  static bool inside(Vec_T const & p) { return p[0] == 0.; }
};

struct RefLineP0
{
  using GeoElem_T = Line;
  using FacetFE_T = RefPointP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 1;
  static uint constexpr numDOFs = 1U;
  static uint constexpr numGeoDOFs = 2U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 1, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 2, 1> constexpr dofOnFacet = {{{{0}}, {{0}}}};
  static Vec1 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<Vec_T, numDOFs> const points;
  static std::array<scalarOnedFun_T, numDOFs> const phiFun;
  static std::array<onedFun_T, numDOFs> const dphiFun;
  static std::array<onedFun_T, numDOFs * dim> const d2phiFun;
  static std::array<onedFun_T, numGeoDOFs> const mapping;

  static double constexpr volume = 2.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{
        {.5 * (e.pointList[0]->coord + e.pointList[1]->coord)}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts = {
        {e.pointList[0]->coord, e.pointList[1]->coord}};
    return mappingPts;
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16;
  }
};

struct RefLineP1
{
  using GeoElem_T = Line;
  using FacetFE_T = RefPointP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 1;
  static uint constexpr numDOFs = 2U;
  static uint constexpr numGeoDOFs = 2U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 2, 1> constexpr dofOnFacet = {{{{0}}, {{1}}}};
  static Vec1 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<Vec_T, numDOFs> const points;
  static std::array<scalarOnedFun_T, numDOFs> const phiFun;
  static std::array<onedFun_T, numDOFs> const dphiFun;
  static std::array<onedFun_T, numDOFs * dim> const d2phiFun;
  static std::array<onedFun_T, numGeoDOFs> const mapping;

  static double constexpr volume = 2.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts = {{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
    }};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16;
  }
};

struct RefLineP2
{
  using GeoElem_T = Line;
  using FacetFE_T = RefPointP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 1;
  static uint constexpr numDOFs = 3U;
  static uint constexpr numGeoDOFs = 3U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 1, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 2, 1> constexpr dofOnFacet = {{{{0}}, {{1}}}};
  static Vec1 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<Vec_T, numDOFs> const points;
  static std::array<scalarOnedFun_T, numDOFs> const phiFun;
  static std::array<onedFun_T, numDOFs> const dphiFun;
  static std::array<onedFun_T, numDOFs * dim> const d2phiFun;
  static std::array<onedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 2.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{
        {e.pointList[0]->coord, e.pointList[1]->coord, e.midpoint()}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16;
  }
};

struct RefTriangleP0
{
  using GeoElem_T = Triangle;
  using FacetFE_T = RefLineP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 1U;
  static uint constexpr numGeoDOFs = 3U;
  static std::array<uint, 4> constexpr dofPlace{{0, 1, 0, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 3, 1> constexpr dofOnFacet = {{{{0}}, {{0}}, {{0}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarTwodFun_T, numDOFs> const phiFun;
  static std::array<twodFun_T, numDOFs> const dphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 0.5L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{{e.midpoint()}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts{
        {e.pointList[0]->coord, e.pointList[1]->coord, e.pointList[2]->coord}};
    return mappingPts;
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > 0. - 1.e-16 && p[1] > 0. - 1.e-16 && p[0] + p[1] < 1. + 1.e-16;
  }
};

struct RefTriangleP1
{
  using GeoElem_T = Triangle;
  using FacetFE_T = RefLineP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 3U;
  static uint constexpr numGeoDOFs = 3U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 2U;
  static array2d<uint, 3, 2> constexpr dofOnFacet = {{{{0, 1}}, {{1, 2}}, {{2, 0}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarTwodFun_T, numDOFs> const phiFun;
  static std::array<twodFun_T, numDOFs> const dphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static short_T constexpr numChildren = 4U;
  static array2d<short_T, numChildren, numDOFs> constexpr childrenIndices = {
      {{{0, 3, 5}}, {{1, 4, 3}}, {{2, 5, 4}}, {{3, 4, 5}}}};
  static std::array<FMat<numDOFs, numDOFs>, numChildren> const embeddingMatrix;

  static double constexpr volume = 0.5L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{
        {e.pointList[0]->coord, e.pointList[1]->coord, e.pointList[2]->coord}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > 0. - 1.e-16 && p[1] > 0. - 1.e-16 && p[0] + p[1] < 1. + 1.e-16;
  }
};

struct RefTriangleP2
{
  using GeoElem_T = Triangle;
  using FacetFE_T = RefLineP2;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 6U;
  static uint constexpr numGeoDOFs = 6U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 1, 1}};
  static uint constexpr dofPerFacet = 3U;
  static uint constexpr numEdges = 3U;
  static uint constexpr dofPerEdge = 3U;
  static array2d<uint, 3, 3> constexpr dofOnFacet = {
      {{{0, 1, 3}}, {{1, 2, 4}}, {{2, 0, 5}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarTwodFun_T, numDOFs> const phiFun;
  static std::array<twodFun_T, numDOFs> const dphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 0.5L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{
        {e.pointList[0]->coord,
         e.pointList[1]->coord,
         e.pointList[2]->coord,
         0.5 * (e.pointList[0]->coord + e.pointList[1]->coord),
         0.5 * (e.pointList[1]->coord + e.pointList[2]->coord),
         0.5 * (e.pointList[2]->coord + e.pointList[0]->coord)}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > 0. - 1.e-16 && p[1] > 0. - 1.e-16 && p[0] + p[1] < 1. + 1.e-16;
  }
};

struct RefTriangleCR1
{
  using GeoElem_T = Triangle;
  using FacetFE_T = RefLineP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 3U;
  static uint constexpr numGeoDOFs = 3U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 1, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 3, 1> constexpr dofOnFacet = {{{{0}}, {{1}}, {{2}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarTwodFun_T, numDOFs> const phiFun;
  static std::array<twodFun_T, numDOFs> const dphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 0.5L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    return std::array<Vec3, numDOFs>{
        0.5 * (e.pointList[0]->coord + e.pointList[1]->coord),
        0.5 * (e.pointList[1]->coord + e.pointList[2]->coord),
        0.5 * (e.pointList[2]->coord + e.pointList[0]->coord)};
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts{
        {e.pointList[0]->coord, e.pointList[1]->coord, e.pointList[2]->coord}};
    return mappingPts;
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > 0. - 1.e-16 && p[1] > 0. - 1.e-16 && p[0] + p[1] < 1. + 1.e-16;
  }
};

struct RefTriangleRT0
{
  using GeoElem_T = Triangle;
  using FacetFE_T = RefLineP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 3U;
  static uint constexpr numGeoDOFs = 3U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 1, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 3, 1> constexpr dofOnFacet = {{{{0}}, {{1}}, {{2}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<twodFun_T, numDOFs> const phiVectFun;
  static std::array<scalarTwodFun_T, numDOFs> const divphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 0.5L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    return std::array<Vec3, numDOFs>{
        0.5 * (e.pointList[0]->coord + e.pointList[1]->coord),
        0.5 * (e.pointList[1]->coord + e.pointList[2]->coord),
        0.5 * (e.pointList[2]->coord + e.pointList[0]->coord)};
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts{
        {e.pointList[0]->coord, e.pointList[1]->coord, e.pointList[2]->coord}};
    return mappingPts;
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > 0. - 1.e-16 && p[1] > 0. - 1.e-16 && p[0] + p[1] < 1. + 1.e-16;
  }
};

struct RefQuadP0
{
  using GeoElem_T = Quad;
  using FacetFE_T = RefLineP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 1U;
  static uint constexpr numGeoDOFs = 4U;
  static std::array<uint, 4> constexpr dofPlace{{0, 1, 0, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 4, 1> constexpr dofOnFacet = {{{{0}}, {{0}}, {{0}}, {{0}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarTwodFun_T, numDOFs> const phiFun;
  static std::array<twodFun_T, numDOFs> const dphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 4.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{{e.midpoint()}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts{{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
        e.pointList[2]->coord,
        e.pointList[3]->coord,
    }};
    return mappingPts;
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16;
  }
};

struct RefQuadQ1
{
  using GeoElem_T = Quad;
  using FacetFE_T = RefLineP1;

  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 4U;
  static uint constexpr numGeoDOFs = 4U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 2U;
  static array2d<uint, 4, 2> constexpr dofOnFacet = {
      {{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarTwodFun_T, numDOFs> const phiFun;
  static std::array<twodFun_T, numDOFs> const dphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static short_T constexpr numChildren = 4U;
  static std::array<FMat<numDOFs, numDOFs>, numChildren> const embeddingMatrix;

  static double constexpr volume = 4.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{
        {e.pointList[0]->coord,
         e.pointList[1]->coord,
         e.pointList[2]->coord,
         e.pointList[3]->coord}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16;
  }
};

struct RefQuadP2
{
  using GeoElem_T = Quad;
  using FacetFE_T = RefLineP2;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 8U;
  static uint constexpr numGeoDOFs = 8U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 1, 1}};
  static uint constexpr dofPerFacet = 3U;
  static array2d<uint, 4, 3> constexpr dofOnFacet = {
      {{{0, 1, 4}}, {{1, 2, 5}}, {{2, 3, 6}}, {{3, 0, 7}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarTwodFun_T, numDOFs> const phiFun;
  static std::array<twodFun_T, numDOFs> const dphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 4.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{
        {e.pointList[0]->coord,
         e.pointList[1]->coord,
         e.pointList[2]->coord,
         e.pointList[3]->coord,
         0.5 * (e.pointList[0]->coord + e.pointList[1]->coord),
         0.5 * (e.pointList[1]->coord + e.pointList[2]->coord),
         0.5 * (e.pointList[2]->coord + e.pointList[3]->coord),
         0.5 * (e.pointList[3]->coord + e.pointList[0]->coord)}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16;
  }
};

struct RefQuadQ2
{
  using GeoElem_T = Quad;
  using FacetFE_T = RefLineP2;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 9U;
  static uint constexpr numGeoDOFs = 9U;
  static std::array<uint, 4> constexpr dofPlace{{0, 1, 1, 1}};
  static uint constexpr dofPerFacet = 3U;
  static array2d<uint, 4, 3> constexpr dofOnFacet = {
      {{{0, 1, 4}}, {{1, 2, 5}}, {{2, 3, 6}}, {{3, 0, 7}}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarTwodFun_T, numDOFs> const phiFun;
  static std::array<twodFun_T, numDOFs> const dphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 4.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{
        {e.pointList[0]->coord,
         e.pointList[1]->coord,
         e.pointList[2]->coord,
         e.pointList[3]->coord,
         0.5 * (e.pointList[0]->coord + e.pointList[1]->coord),
         0.5 * (e.pointList[1]->coord + e.pointList[2]->coord),
         0.5 * (e.pointList[2]->coord + e.pointList[3]->coord),
         0.5 * (e.pointList[3]->coord + e.pointList[0]->coord),
         e.midpoint()}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16;
  }
};

struct RefQuadRT0
{
  using GeoElem_T = Quad;
  using FacetFE_T = RefLineP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numDOFs = 4U;
  static uint constexpr numGeoDOFs = 4U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 1, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 4, 1> constexpr dofOnFacet = {{{0}, {1}, {2}, {3}}};
  static Vec2 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<twodFun_T, numDOFs> const phiVectFun;
  static std::array<scalarTwodFun_T, numDOFs> const divphiFun;
  static std::array<twodFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 4.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    return std::array<Vec3, numDOFs>{
        0.5 * (e.pointList[0]->coord + e.pointList[1]->coord),
        0.5 * (e.pointList[1]->coord + e.pointList[2]->coord),
        0.5 * (e.pointList[2]->coord + e.pointList[3]->coord),
        0.5 * (e.pointList[3]->coord + e.pointList[0]->coord),
    };
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return std::array<Vec3, numGeoDOFs>{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
        e.pointList[2]->coord,
        e.pointList[3]->coord,
    };
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16;
  }
};

struct RefTetrahedronP0
{
  using GeoElem_T = Tetrahedron;
  using FacetFE_T = RefTriangleP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 3;
  static uint constexpr numDOFs = 1U;
  static uint constexpr numGeoDOFs = 4U;
  static std::array<uint, 4> constexpr dofPlace{{1, 0, 0, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 4, 1> constexpr dofOnFacet = {{{{0}}, {{0}}, {{0}}, {{0}}}};
  static Vec3 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarThreedFun_T, numDOFs> const phiFun;
  static std::array<threedFun_T, numDOFs> const dphiFun;
  static std::array<threedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 1. / 6;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{{e.midpoint()}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts{{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
        e.pointList[2]->coord,
        e.pointList[3]->coord,
    }};
    return mappingPts;
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > 0. - 1.e-16 && p[1] > 0. - 1.e-16 && p[2] > 0. - 1.e-16 &&
           p[0] + p[1] + p[2] < 1. + 1.e-16;
  }
};

struct RefTetrahedronP1
{
  using GeoElem_T = Tetrahedron;
  using FacetFE_T = RefTriangleP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 3;
  static uint constexpr numDOFs = 4U;
  static uint constexpr numGeoDOFs = 4U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 3U;
  static array2d<uint, 4, 3> constexpr dofOnFacet = {
      {{{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}}, {{1, 2, 3}}}};
  static Vec3 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarThreedFun_T, numDOFs> const phiFun;
  static std::array<threedFun_T, numDOFs> const dphiFun;
  static std::array<threedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 1. / 6;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{
        {e.pointList[0]->coord,
         e.pointList[1]->coord,
         e.pointList[2]->coord,
         e.pointList[3]->coord}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > 0. - 1.e-16 && p[1] > 0. - 1.e-16 && p[2] > 0. - 1.e-16 &&
           p[0] + p[1] + p[2] < 1. + 1.e-16;
  }
};

struct RefTetrahedronP2
{
  using GeoElem_T = Tetrahedron;
  using FacetFE_T = RefTriangleP2;
  static GeoElem_T const geoElem;
  static int constexpr dim = 3;
  static uint constexpr numDOFs = 10U;
  static uint constexpr numGeoDOFs = 10U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 1, 1}};
  static uint constexpr dofPerFacet = 6U;
  static array2d<uint, 4, 6> constexpr dofOnFacet = {
      {{{0, 2, 1, 6, 5, 4}},
       {{0, 1, 3, 4, 8, 7}},
       {{0, 3, 2, 7, 9, 6}},
       {{1, 2, 3, 5, 9, 8}}}};
  static Vec3 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarThreedFun_T, numDOFs> const phiFun;
  static std::array<threedFun_T, numDOFs> const dphiFun;
  static std::array<threedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 1. / 6;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts = {{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
        e.pointList[2]->coord,
        e.pointList[3]->coord,
        .5 * (e.pointList[0]->coord + e.pointList[1]->coord),
        .5 * (e.pointList[1]->coord + e.pointList[2]->coord),
        .5 * (e.pointList[2]->coord + e.pointList[0]->coord),
        .5 * (e.pointList[0]->coord + e.pointList[3]->coord),
        .5 * (e.pointList[1]->coord + e.pointList[3]->coord),
        .5 * (e.pointList[2]->coord + e.pointList[3]->coord),
    }};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > 0. - 1.e-16 && p[1] > 0. - 1.e-16 && p[2] > 0. - 1.e-16 &&
           p[0] + p[1] + p[2] < 1. + 1.e-16;
  }
};

struct RefTetrahedronRT0
{
  using GeoElem_T = Tetrahedron;
  using FacetFE_T = RefTriangleP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 3;
  static uint constexpr numDOFs = 4U;
  static uint constexpr numGeoDOFs = 4U;
  static std::array<uint, 4> constexpr dofPlace{{0, 1, 0, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 4, 1> constexpr dofOnFacet = {{{{0}}, {{1}}, {{2}}, {{3}}}};
  static Vec3 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<threedFun_T, numDOFs> const phiVectFun;
  static std::array<scalarThreedFun_T, numDOFs> const divphiFun;
  static std::array<threedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 0.1666666666666667L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    return std::array<Vec3, numDOFs>{
        (e.pointList[0]->coord + e.pointList[2]->coord + e.pointList[1]->coord) / 3.,
        (e.pointList[0]->coord + e.pointList[1]->coord + e.pointList[3]->coord) / 3.,
        (e.pointList[0]->coord + e.pointList[3]->coord + e.pointList[2]->coord) / 3.,
        (e.pointList[1]->coord + e.pointList[2]->coord + e.pointList[3]->coord) / 3.,
    };
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts{
        {e.pointList[0]->coord,
         e.pointList[1]->coord,
         e.pointList[2]->coord,
         e.pointList[3]->coord}};
    return mappingPts;
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > (0. - 1.e-16) && p[1] > (0. - 1.e-16) && p[2] > (0. - 1.e-16) &&
           p[0] + p[1] + p[2] < (1. + 1.e-16);
  }
};

struct RefHexahedronP0
{
  using GeoElem_T = Hexahedron;
  using FacetFE_T = RefQuadP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 3;
  static uint constexpr numDOFs = 1U;
  static uint constexpr numGeoDOFs = 8U;
  static std::array<uint, 4> constexpr dofPlace{{1, 0, 0, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 6, 1> constexpr dofOnFacet = {
      {{{0}}, {{0}}, {{0}}, {{0}}, {{0}}, {{0}}}};
  static Vec3 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarThreedFun_T, numDOFs> const phiFun;
  static std::array<threedFun_T, numDOFs> const dphiFun;
  static std::array<threedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 8.;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts{{e.midpoint()}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    std::array<Vec3, numGeoDOFs> mappingPts{{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
        e.pointList[2]->coord,
        e.pointList[3]->coord,
        e.pointList[4]->coord,
        e.pointList[5]->coord,
        e.pointList[6]->coord,
        e.pointList[7]->coord,
    }};
    return mappingPts;
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16 && p[2] > -1. - 1.e-16 && p[2] < 1. + 1.e-16;
  }
};

struct RefHexahedronQ1
{
  using GeoElem_T = Hexahedron;
  using FacetFE_T = RefQuadQ1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 3;
  static uint constexpr numDOFs = 8U;
  static uint constexpr numGeoDOFs = 8U;
  static std::array<uint, 4> constexpr dofPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 4U;
  static array2d<uint, 6, 4> constexpr dofOnFacet = {
      {{{0, 3, 2, 1}},
       {{0, 1, 5, 4}},
       {{1, 2, 6, 5}},
       {{2, 3, 7, 6}},
       {{3, 0, 4, 7}},
       {{4, 5, 6, 7}}}};
  static Vec3 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarThreedFun_T, numDOFs> const phiFun;
  static std::array<threedFun_T, numDOFs> const dphiFun;
  static std::array<threedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 8.;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts = {{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
        e.pointList[2]->coord,
        e.pointList[3]->coord,
        e.pointList[4]->coord,
        e.pointList[5]->coord,
        e.pointList[6]->coord,
        e.pointList[7]->coord,
    }};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16 && p[2] > -1. - 1.e-16 && p[2] < 1. + 1.e-16;
  }
};

struct RefHexahedronQ2
{
  using GeoElem_T = Hexahedron;
  using FacetFE_T = RefQuadQ2;
  static GeoElem_T const geoElem;
  static int constexpr dim = 3;
  static uint constexpr numDOFs = 27U;
  static uint constexpr numGeoDOFs = 27U;
  static std::array<uint, 4> constexpr dofPlace{{1, 1, 1, 1}};
  static uint constexpr dofPerFacet = 9U;
  static array2d<uint, 6, 9> constexpr dofOnFacet = {
      {{{0, 3, 2, 1, 11, 10, 9, 8, 20}},
       {{0, 1, 5, 4, 8, 13, 16, 12, 21}},
       {{1, 2, 6, 5, 9, 14, 17, 13, 22}},
       {{2, 3, 7, 6, 10, 15, 18, 14, 23}},
       {{3, 0, 4, 7, 11, 12, 19, 15, 24}},
       {{4, 5, 6, 7, 16, 17, 18, 19, 25}}}};
  static Vec3 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<scalarThreedFun_T, numDOFs> const phiFun;
  static std::array<threedFun_T, numDOFs> const dphiFun;
  static std::array<threedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 8.;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    std::array<Vec3, numDOFs> dofPts = {
        {e.pointList[0]->coord,
         e.pointList[1]->coord,
         e.pointList[2]->coord,
         e.pointList[3]->coord,
         e.pointList[4]->coord,
         e.pointList[5]->coord,
         e.pointList[6]->coord,
         e.pointList[7]->coord,
         0.5 * (e.pointList[0]->coord + e.pointList[1]->coord),
         0.5 * (e.pointList[1]->coord + e.pointList[2]->coord),
         0.5 * (e.pointList[2]->coord + e.pointList[3]->coord),
         0.5 * (e.pointList[3]->coord + e.pointList[0]->coord),
         0.5 * (e.pointList[0]->coord + e.pointList[4]->coord),
         0.5 * (e.pointList[1]->coord + e.pointList[5]->coord),
         0.5 * (e.pointList[2]->coord + e.pointList[6]->coord),
         0.5 * (e.pointList[3]->coord + e.pointList[7]->coord),
         0.5 * (e.pointList[4]->coord + e.pointList[5]->coord),
         0.5 * (e.pointList[5]->coord + e.pointList[6]->coord),
         0.5 * (e.pointList[6]->coord + e.pointList[7]->coord),
         0.5 * (e.pointList[7]->coord + e.pointList[4]->coord),
         0.25 * (e.pointList[0]->coord + e.pointList[1]->coord + e.pointList[2]->coord +
                 e.pointList[3]->coord),
         0.25 * (e.pointList[0]->coord + e.pointList[1]->coord + e.pointList[5]->coord +
                 e.pointList[4]->coord),
         0.25 * (e.pointList[1]->coord + e.pointList[2]->coord + e.pointList[6]->coord +
                 e.pointList[5]->coord),
         0.25 * (e.pointList[2]->coord + e.pointList[3]->coord + e.pointList[7]->coord +
                 e.pointList[6]->coord),
         0.25 * (e.pointList[3]->coord + e.pointList[0]->coord + e.pointList[4]->coord +
                 e.pointList[7]->coord),
         0.25 * (e.pointList[4]->coord + e.pointList[5]->coord + e.pointList[6]->coord +
                 e.pointList[7]->coord),
         e.midpoint()}};
    return dofPts;
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return dofPts(e);
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16 && p[2] > -1. - 1.e-16 && p[2] < 1. + 1.e-16;
  }
};

struct RefHexahedronRT0
{
  using GeoElem_T = Hexahedron;
  using FacetFE_T = RefQuadP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 3;
  static uint constexpr numDOFs = 6U;
  static uint constexpr numGeoDOFs = 8U;
  static std::array<uint, 4> constexpr dofPlace{{0, 1, 0, 0}};
  static std::array<uint, 4> inline constexpr geoPlace{{0, 0, 0, 1}};
  static uint constexpr dofPerFacet = 1U;
  static array2d<uint, 6, 1> constexpr dofOnFacet = {{{0}, {1}, {2}, {3}, {4}, {5}}};
  static Vec3 const refMidpoint;
  using Vec_T = FVec<dim>;

  static std::array<threedFun_T, numDOFs> const phiVectFun;
  static std::array<scalarThreedFun_T, numDOFs> const divphiFun;
  static std::array<threedFun_T, numGeoDOFs> const mapping;
  static double constexpr volume = 8.L;

  static std::array<Vec3, numDOFs> dofPts(GeoElem const & e)
  {
    return std::array<Vec3, numDOFs>{
        0.25 * (e.pointList[0]->coord + e.pointList[3]->coord + e.pointList[2]->coord +
                e.pointList[1]->coord),
        0.25 * (e.pointList[0]->coord + e.pointList[1]->coord + e.pointList[5]->coord +
                e.pointList[4]->coord),
        0.25 * (e.pointList[1]->coord + e.pointList[2]->coord + e.pointList[6]->coord +
                e.pointList[5]->coord),
        0.25 * (e.pointList[2]->coord + e.pointList[3]->coord + e.pointList[7]->coord +
                e.pointList[6]->coord),
        0.25 * (e.pointList[3]->coord + e.pointList[0]->coord + e.pointList[4]->coord +
                e.pointList[7]->coord),
        0.25 * (e.pointList[4]->coord + e.pointList[5]->coord + e.pointList[6]->coord +
                e.pointList[7]->coord),
    };
  }

  static std::array<Vec3, numGeoDOFs> mappingPts(GeoElem const & e)
  {
    return std::array<Vec3, numGeoDOFs>{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
        e.pointList[2]->coord,
        e.pointList[3]->coord,
        e.pointList[4]->coord,
        e.pointList[5]->coord,
        e.pointList[6]->coord,
        e.pointList[7]->coord,
    };
  }

  static bool inside(Vec_T const & p)
  {
    return p[0] > -1. - 1.e-16 && p[0] < 1. + 1.e-16 && p[1] > -1. - 1.e-16 &&
           p[1] < 1. + 1.e-16 && p[2] > -1. - 1.e-16 && p[2] < 1. + 1.e-16;
  }
};

// ----------------------------------------------------------------------------
template <typename RefFE>
struct Order
{
  static constexpr uint value = 0;
};

template <>
struct Order<RefLineP0>
{
  static constexpr uint value = 0;
};
template <>
struct Order<RefPointP1>
{
  static constexpr uint value = 1;
};
template <>
struct Order<RefLineP1>
{
  static constexpr uint value = 1;
};
template <>
struct Order<RefTriangleP1>
{
  static constexpr uint value = 1;
};
template <>
struct Order<RefQuadQ1>
{
  static constexpr uint value = 1;
};
template <>
struct Order<RefLineP2>
{
  static constexpr uint value = 2;
};
template <>
struct Order<RefTriangleP2>
{
  static constexpr uint value = 2;
};
template <>
struct Order<RefQuadP2>
{
  static constexpr uint value = 2;
};
template <>
struct Order<RefQuadQ2>
{
  static constexpr uint value = 2;
};
template <>
struct Order<RefTetrahedronP1>
{
  static constexpr uint value = 1;
};
template <>
struct Order<RefTetrahedronP2>
{
  static constexpr uint value = 2;
};
template <>
struct Order<RefHexahedronQ1>
{
  static constexpr uint value = 1;
};
template <>
struct Order<RefHexahedronQ2>
{
  static constexpr uint value = 2;
};

template <typename RefFE>
inline constexpr uint order_v = Order<RefFE>::value;

// ----------------------------------------------------------------------------
enum class FamilyType : int8_t
{
  NONE = 0,
  LAGRANGE = 1,
  RAVIART_THOMAS = 2,
  CROUZEIX_RAVIART = 3,
};

template <typename RefFE>
struct Family
{};

template <>
struct Family<RefPointP1>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefLineP0>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefLineP1>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefLineP2>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefTriangleP0>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefTriangleP1>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefTriangleP2>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefTriangleCR1>
{
  static constexpr FamilyType value = FamilyType::CROUZEIX_RAVIART;
};
template <>
struct Family<RefTriangleRT0>
{
  static constexpr FamilyType value = FamilyType::RAVIART_THOMAS;
};
template <>
struct Family<RefQuadP0>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefQuadQ1>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefQuadP2>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefQuadQ2>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefQuadRT0>
{
  static constexpr FamilyType value = FamilyType::RAVIART_THOMAS;
};
template <>
struct Family<RefTetrahedronP0>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefTetrahedronP1>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefTetrahedronP2>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefTetrahedronRT0>
{
  static constexpr FamilyType value = FamilyType::RAVIART_THOMAS;
};
template <>
struct Family<RefHexahedronP0>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefHexahedronQ1>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefHexahedronQ2>
{
  static constexpr FamilyType value = FamilyType::LAGRANGE;
};
template <>
struct Family<RefHexahedronRT0>
{
  static constexpr FamilyType value = FamilyType::RAVIART_THOMAS;
};

template <typename RefFE>
inline constexpr FamilyType family_v = Family<RefFE>::value;

// ----------------------------------------------------------------------------
enum class FEDimType : int8_t
{
  NONE = 0,
  SCALAR = 1,
  VECTOR = 2
};

template <typename RefFE>
struct FEDim
{};

template <>
struct FEDim<RefPointP1>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefLineP0>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefLineP1>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefLineP2>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefTriangleP0>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefTriangleP1>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefTriangleP2>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefTriangleCR1>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefTriangleRT0>
{
  static constexpr FEDimType value = FEDimType::VECTOR;
};
template <>
struct FEDim<RefQuadP0>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefQuadQ1>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefQuadP2>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefQuadQ2>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefQuadRT0>
{
  static constexpr FEDimType value = FEDimType::VECTOR;
};
template <>
struct FEDim<RefTetrahedronP0>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefTetrahedronP1>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefTetrahedronP2>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefTetrahedronRT0>
{
  static constexpr FEDimType value = FEDimType::VECTOR;
};
template <>
struct FEDim<RefHexahedronP0>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefHexahedronQ1>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefHexahedronQ2>
{
  static constexpr FEDimType value = FEDimType::SCALAR;
};
template <>
struct FEDim<RefHexahedronRT0>
{
  static constexpr FEDimType value = FEDimType::VECTOR;
};

template <typename RefFE>
inline constexpr FEDimType fedim_v = FEDim<RefFE>::value;

template <typename RefFE>
static constexpr int feDimValue()
{
  if constexpr (fedim_v<RefFE> == FEDimType::SCALAR)
    return 1;
  else if constexpr (fedim_v<RefFE> == FEDimType::VECTOR)
    return RefFE::dim;
  else
    abort();
  return 0;
}

// ----------------------------------------------------------------------------
template <typename RefFE>
struct MappingIsSeparate
{
  static constexpr bool value = false;
};

template <>
struct MappingIsSeparate<RefLineP0>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefTriangleP0>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefTriangleCR1>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefTriangleRT0>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefQuadP0>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefQuadRT0>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefTetrahedronP0>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefTetrahedronRT0>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefHexahedronP0>
{
  static constexpr bool value = true;
};
template <>
struct MappingIsSeparate<RefHexahedronRT0>
{
  static constexpr bool value = true;
};

template <typename RefFE>
inline constexpr bool mappingIsSeparate_v = MappingIsSeparate<RefFE>::value;

// ----------------------------------------------------------------------------
template <typename RefFE>
struct EnableSecondDeriv
{
  bool static constexpr value = false;
};

#ifdef MINIFEM_ENABLE_SECONDDERIV
template <>
struct EnableSecondDeriv<RefLineP1>
{
  bool static constexpr value = true;
};

template <>
struct EnableSecondDeriv<RefLineP2>
{
  bool static constexpr value = true;
};
#endif

template <typename RefFE>
inline bool static constexpr enableSecondDeriv_v = EnableSecondDeriv<RefFE>::value;
