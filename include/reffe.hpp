#pragma once

#include "def.hpp"
#include "geo.hpp"

template <typename RefElem>
uint constexpr numDOFs()
{
  return 1 * RefElem::dof_place[0] +
    RefElem::GeoElem_T::numFaces * RefElem::dof_place[1] +
    RefElem::GeoElem_T::numEdges * RefElem::dof_place[2] +
    RefElem::GeoElem_T::numPts * RefElem::dof_place[3];
}

struct RefPointP1
{
  using GeoElem_T = PointElem;
  using RefFacet_T = NullElem;
  static GeoElem_T const geoElem;
  static int constexpr dim = 1;
  static uint constexpr numFuns = 1U;
  static uint constexpr numGeoFuns = 1U;
  static array<uint,4> constexpr dof_place{{0,0,0,1}};
  static uint constexpr dofPerFacet = 0U;
  static array<array<uint,0>,0> constexpr dofOnFacet = {};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<1>;

  static array<Vec_T,numFuns> const points;
  static array<scalarOnedFun_T,numFuns> const phiFun;
  static array<onedFun_T,numFuns> const phiVectFun;
  static array<onedFun_T,numFuns> const dphiFun;
  static array<onedFun_T,numGeoFuns> const mapping;
  static double constexpr volume = 1.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts =
     {{
       e.pointList[0]->coord
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
    array<Vec3,numGeoFuns> mappingPts = {{
        e.pointList[0]->coord
    }};
    return mappingPts;
  }
};

struct RefLineP0
{
  using GeoElem_T = Line;
  using RefFacet_T = NullElem;
  static GeoElem_T const geoElem;
  static int constexpr dim = 1;
  static uint constexpr numFuns = 1U;
  static uint constexpr numGeoFuns = 2U;
  static array<uint,4> constexpr dof_place{{0,0,1,0}};
  static uint constexpr dofPerFacet = 0U;
  static array<array<uint,0>,0> constexpr dofOnFacet = {};
  using Vec_T = FVec<dim>;
//  using LocalVec_T = FVec<2>;
//  using LocalMat_T = FMat<2,2>;

  static array<Vec_T,numFuns> const points;
  static array<scalarOnedFun_T,numFuns> const phiFun;
  static array<onedFun_T,numFuns> const phiVectFun;
  static array<onedFun_T,numFuns> const dphiFun;
  static array<onedFun_T,numGeoFuns> const mapping;
  // static LocalMat_T const massMat;
  // static LocalMat_T const gradMat;
  static double constexpr volume = 2.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts {{
       .5*(e.pointList[0]->coord + e.pointList[1]->coord)
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
    array<Vec3,numGeoFuns> mappingPts = {{
        e.pointList[0]->coord,
        e.pointList[1]->coord
    }};
    return mappingPts;
  }
};

struct RefLineP1
{
  using GeoElem_T = Line;
  using RefFacet_T = RefPointP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 1;
  static uint constexpr numFuns = 2U;
  static uint constexpr numGeoFuns = 2U;
  static array<uint,4> constexpr dof_place{{0,0,0,1}};
  static uint constexpr dofPerFacet = 1U;
  static array<array<uint,1>,2> constexpr dofOnFacet = {{
    {{0}}, {{1}}
  }};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<2>;
  using LocalMat_T = FMat<2,2>;

  static array<Vec_T,numFuns> const points;
  static array<scalarOnedFun_T,numFuns> const phiFun;
  static array<onedFun_T,numFuns> const phiVectFun;
  static array<onedFun_T,numFuns> const dphiFun;
  static array<onedFun_T,numGeoFuns> const mapping;
  // static LocalMat_T const massMat;
  // static LocalMat_T const gradMat;
  static double constexpr volume = 2.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts = {{
       e.pointList[0]->coord,
       e.pointList[1]->coord,
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
     return dofPts(e);
  }
};

struct RefLineP2
{
  using GeoElem_T = Line;
  using RefFacet_T = RefPointP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 1;
  static uint constexpr numFuns = 3U;
  static uint constexpr numGeoFuns = 3U;
  static array<uint,4> constexpr dof_place{{0,0,1,1}};
  static uint constexpr dofPerFacet = 1U;
  static array<array<uint,1>,2> constexpr dofOnFacet = {{
    {{0}}, {{1}}
  }};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<numFuns>;
  using LocalMat_T = FMat<numFuns,numFuns>;

  static array<Vec_T,numFuns> const points;
  static array<scalarOnedFun_T,numFuns> const phiFun;
  static array<onedFun_T,numFuns> const phiVectFun;
  static array<onedFun_T,numFuns> const dphiFun;
  static array<onedFun_T,numGeoFuns> const mapping;
  static LocalMat_T const massMat;
  static LocalMat_T const gradMat;
  static double constexpr volume = 2.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts {{
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.midpoint()
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
     return dofPts(e);
  }
};

struct RefTriangleP0
{
  using GeoElem_T = Triangle;
  using RefFacet_T = RefLineP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numFuns = 1U;
  static uint constexpr numGeoFuns = 3U;
  static array<uint,4> constexpr dof_place{{0,1,0,0}};
  static uint constexpr dofPerFacet = 0U;
  static array<array<uint,0>,0> constexpr dofOnFacet = {};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<numFuns>;
  using LocalMat_T = FMat<numFuns,numFuns>;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const phiVectFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static array<twodFun_T,numGeoFuns> const mapping;
  static double constexpr volume = 0.5L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
    array<Vec3,numFuns> dofPts {{
        e.midpoint()
    }};
    return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
    array<Vec3,numGeoFuns> mappingPts {{
        e.pointList[0]->coord,
        e.pointList[1]->coord,
        e.pointList[2]->coord
    }};
    return mappingPts;
  }
};

struct RefTriangleP1
{
  using GeoElem_T = Triangle;
  using RefFacet_T = RefLineP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numFuns = 3U;
  static uint constexpr numGeoFuns = 3U;
  static array<uint,4> constexpr dof_place{{0,0,0,1}};
  static uint constexpr dofPerFacet = 2U;
  static array<array<uint,2>,3> constexpr dofOnFacet = {{
    {{0,1}}, {{1,2}}, {{2,0}}
  }};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<numFuns>;
  using LocalMat_T = FMat<numFuns,numFuns>;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const phiVectFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static array<twodFun_T,numGeoFuns> const mapping;
  static double constexpr volume = 0.5L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts {{
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
     return dofPts(e);
  }
};

struct RefTriangleP2
{
  using GeoElem_T = Triangle;
  using RefFacet_T = RefLineP2;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numFuns = 6U;
  static uint constexpr numGeoFuns = 6U;
  static array<uint,4> constexpr dof_place{{0,0,1,1}};
  static uint constexpr dofPerFacet = 3U;
  static array<array<uint,3>,3> constexpr dofOnFacet = {{
    {{0,1,3}}, {{1,2,4}}, {{2,0,5}}
  }};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<numFuns>;
  using LocalMat_T = FMat<numFuns,numFuns>;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const phiVectFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static array<twodFun_T,numGeoFuns> const mapping;
  static double constexpr volume = 0.5L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts {{
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord,
       0.5*(e.pointList[0]->coord+e.pointList[1]->coord),
       0.5*(e.pointList[1]->coord+e.pointList[2]->coord),
       0.5*(e.pointList[2]->coord+e.pointList[0]->coord)
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
     return dofPts(e);
  }
};

struct RefTriangleRT0
{
  using GeoElem_T = Triangle;
  using RefFacet_T = RefLineP0;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numFuns = 3U;
  static uint constexpr numGeoFuns = 3U;
  static array<uint,4> constexpr dof_place{{0,0,1,0}};
  static uint constexpr dofPerFacet = 1U;
  static array<array<uint,1>,3> constexpr dofOnFacet = {{
    {{0}}, {{1}}, {{2}}
  }};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<numFuns>;
  using LocalMat_T = FMat<numFuns,numFuns>;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const phiVectFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static array<twodFun_T,numGeoFuns> const mapping;
  static double constexpr volume = 0.5L;

  static array<Vec3,numFuns> dofPts(GeoElem const &)
  {
    return array<Vec3,numFuns>{};
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
     array<Vec3,numGeoFuns> mappingPts {{
         e.pointList[0]->coord,
         e.pointList[1]->coord,
         e.pointList[2]->coord
     }};
     return mappingPts;
  }
};

struct RefQuadQ1
{
  using GeoElem_T = Quad;
  using RefFacet_T = RefLineP1;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numFuns = 4U;
  static uint constexpr numGeoFuns = 4U;
  static array<uint,4> constexpr dof_place{{0,0,0,1}};
  static uint constexpr dofPerFacet = 2U;
  static array<array<uint,2>,4> constexpr dofOnFacet = {{
    {{0,1}}, {{1,2}}, {{2,3}}, {{3,0}}
  }};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<numFuns>;
  using LocalMat_T = FMat<numFuns,numFuns>;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const phiVectFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static array<twodFun_T,numGeoFuns> const mapping;
  static double constexpr volume = 4.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts {{
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord,
       e.pointList[3]->coord
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
     return dofPts(e);
  }
};

struct RefQuadP2
{
  using GeoElem_T = Quad;
  using RefFacet_T = RefLineP2;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numFuns = 8U;
  static uint constexpr numGeoFuns = 8U;
  static array<uint,4> constexpr dof_place{{0,0,1,1}};
  static uint constexpr dofPerFacet = 3U;
  static array<array<uint,3>,4> constexpr dofOnFacet = {{
    {{0,1,4}}, {{1,2,5}}, {{2,3,6}}, {{3,0,7}}
  }};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<numFuns>;
  using LocalMat_T = FMat<numFuns,numFuns>;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const phiVectFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static array<twodFun_T,numGeoFuns> const mapping;
  static double constexpr volume = 4.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts {{
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord,
       e.pointList[3]->coord,
       0.5*(e.pointList[0]->coord+e.pointList[1]->coord),
       0.5*(e.pointList[1]->coord+e.pointList[2]->coord),
       0.5*(e.pointList[2]->coord+e.pointList[3]->coord),
       0.5*(e.pointList[3]->coord+e.pointList[0]->coord)
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
     return dofPts(e);
  }
};

struct RefQuadQ2
{
  using GeoElem_T = Quad;
  using RefFacet_T = RefLineP2;
  static GeoElem_T const geoElem;
  static int constexpr dim = 2;
  static uint constexpr numFuns = 9U;
  static uint constexpr numGeoFuns = 9U;
  static array<uint,4> constexpr dof_place{{0,1,1,1}};
  static uint constexpr dofPerFacet = 3U;
  static array<array<uint,3>,4> constexpr dofOnFacet = {{
    {{0,1,4}}, {{1,2,5}}, {{2,3,6}}, {{3,0,7}}
  }};
  using Vec_T = FVec<dim>;
  using LocalVec_T = FVec<numFuns>;
  using LocalMat_T = FMat<numFuns,numFuns>;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const phiVectFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static array<twodFun_T,numGeoFuns> const mapping;
  static double constexpr volume = 4.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
    array<Vec3,numFuns> dofPts {{
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord,
       e.pointList[3]->coord,
       0.5*(e.pointList[0]->coord+e.pointList[1]->coord),
       0.5*(e.pointList[1]->coord+e.pointList[2]->coord),
       0.5*(e.pointList[2]->coord+e.pointList[3]->coord),
       0.5*(e.pointList[3]->coord+e.pointList[0]->coord),
       e.midpoint()
     }};
     return dofPts;
  }

  static array<Vec3,numGeoFuns> mappingPts(GeoElem const & e)
  {
     return dofPts(e);
  }
};

template <typename RefFE>
struct Order{ static constexpr uint value = 0; };

template <>
struct Order<RefLineP0>{ static constexpr uint value = 0; };
template <>
struct Order<RefPointP1>{ static constexpr uint value = 1; };
template <>
struct Order<RefLineP1>{ static constexpr uint value = 1; };
template <>
struct Order<RefTriangleP1>{ static constexpr uint value = 1; };
template <>
struct Order<RefQuadQ1>{ static constexpr uint value = 1; };
template <>
struct Order<RefLineP2>{ static constexpr uint value = 2; };
template <>
struct Order<RefTriangleP2>{ static constexpr uint value = 2; };
template <>
struct Order<RefQuadP2>{ static constexpr uint value = 2; };
template <>
struct Order<RefQuadQ2>{ static constexpr uint value = 2; };
