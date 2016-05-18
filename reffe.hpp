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
  static uint constexpr dim = 0U;
  static uint constexpr numFuns = 1U;
  typedef FVec<dim> Vec_T;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts =
     {
       e.pointList[0]->coord
     };
     return std::move(dofPts);
  }
};

struct RefLineP1
{
  typedef Line GeoElem_T;
  typedef RefPointP1 RefFacet_T;
  static GeoElem_T const geoElem;
  static uint constexpr dim = 1U;
  static uint constexpr numFuns = 2U;
  static array<uint,4> constexpr dof_place{0,0,0,1};
  static uint constexpr dofPerFacet = 1U;
  static array<array<uint,1>,2> constexpr dofOnFacet = {
    {{0}, {1}}
  };
  typedef FVec<dim> Vec_T;
  typedef FVec<2> LocalVec_T;
  typedef FMat<2,2> LocalMat_T;

  static array<Vec_T,numFuns> const points;
  static array<ScalarFun<dim>,numFuns> const phiFun;
  static array<Fun<dim,dim>,numFuns> const dphiFun;
  static LocalMat_T const massMat;
  static LocalMat_T const gradMat;
  static double constexpr volume = 2.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts =
     {
       e.pointList[0]->coord,
       e.pointList[1]->coord,
     };
     return std::move(dofPts);
  }
};

struct RefLineP2
{
  typedef Line GeoElem_T;
  typedef RefPointP1 RefFacet_T;
  static GeoElem_T const geoElem;
  static uint constexpr dim = 1U;
  static uint constexpr numFuns = 3U;
  static array<uint,4> constexpr dof_place{0,0,1,1};
  static uint constexpr dofPerFacet = 1U;
  static array<array<uint,1>,2> constexpr dofOnFacet = {
    {{0}, {1}}
  };
  typedef FVec<dim> Vec_T;
  typedef FVec<numFuns> LocalVec_T;
  typedef FMat<numFuns,numFuns> LocalMat_T;

  static array<Vec_T,numFuns> const points;
  static array<ScalarFun<dim>,numFuns> const phiFun;
  static array<Fun<dim,dim>,numFuns> const dphiFun;
  static LocalMat_T const massMat;
  static LocalMat_T const gradMat;
  static double constexpr volume = 2.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts =
     {
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.midpoint()
     };
     return std::move(dofPts);
  }
};

struct RefTriangleP1
{
  typedef Triangle GeoElem_T;
  typedef RefLineP1 RefFacet_T;
  static GeoElem_T const geoElem;
  static uint constexpr dim = 2U;
  static uint constexpr numFuns = 3U;
  static array<uint,4> constexpr dof_place{0,0,0,1};
  static uint constexpr dofPerFacet = 2U;
  static array<array<uint,2>,3> constexpr dofOnFacet = {
    {{0,1}, {1,2}, {2,0}}
  };
  typedef FVec<dim> Vec_T;
  typedef FVec<numFuns> LocalVec_T;
  typedef FMat<numFuns,numFuns> LocalMat_T;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static double constexpr volume = 0.5L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts =
     {
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord
     };
     return std::move(dofPts);
  }
};

struct RefTriangleP2
{
  typedef Triangle GeoElem_T;
  typedef RefLineP2 RefFacet_T;
  static GeoElem_T const geoElem;
  static uint constexpr dim = 2U;
  static uint constexpr numFuns = 6U;
  static array<uint,4> constexpr dof_place{0,0,1,1};
  static uint constexpr dofPerFacet = 3U;
  static array<array<uint,3>,3> constexpr dofOnFacet = {
    {{0,1,3}, {1,2,4}, {2,0,5}}
  };
  typedef FVec<dim> Vec_T;
  typedef FVec<numFuns> LocalVec_T;
  typedef FMat<numFuns,numFuns> LocalMat_T;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static double constexpr volume = 0.5L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts =
     {
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord,
       0.5*(e.pointList[0]->coord+e.pointList[1]->coord),
       0.5*(e.pointList[1]->coord+e.pointList[2]->coord),
       0.5*(e.pointList[2]->coord+e.pointList[0]->coord)
     };
     return std::move(dofPts);
  }
};

struct RefQuadQ1
{
  typedef Quad GeoElem_T;
  typedef RefLineP1 RefFacet_T;
  static GeoElem_T const geoElem;
  static uint constexpr dim = 2U;
  static uint constexpr numFuns = 4U;
  static array<uint,4> constexpr dof_place{0,0,0,1};
  static uint constexpr dofPerFacet = 2U;
  static array<array<uint,2>,4> constexpr dofOnFacet = {
    {{0,1}, {1,2}, {2,3}, {3,0}}
  };
  typedef FVec<dim> Vec_T;
  typedef FVec<numFuns> LocalVec_T;
  typedef FMat<numFuns,numFuns> LocalMat_T;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static double constexpr volume = 4.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts =
     {
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord,
       e.pointList[3]->coord
     };
     return std::move(dofPts);
  }
};

struct RefQuadQ2
{
  typedef Quad GeoElem_T;
  typedef RefLineP2 RefFacet_T;
  static GeoElem_T const geoElem;
  static uint constexpr dim = 2U;
  static uint constexpr numFuns = 9U;
  static array<uint,4> constexpr dof_place{0,1,1,1};
  static uint constexpr dofPerFacet = 3U;
  static array<array<uint,3>,4> constexpr dofOnFacet = {
    {{0,1,4}, {1,2,5}, {2,3,6}, {3,0,7}}
  };
  typedef FVec<dim> Vec_T;
  typedef FVec<numFuns> LocalVec_T;
  typedef FMat<numFuns,numFuns> LocalMat_T;

  static array<scalarTwodFun_T,numFuns> const phiFun;
  static array<twodFun_T,numFuns> const dphiFun;
  static double constexpr volume = 4.L;

  static array<Vec3,numFuns> dofPts(GeoElem const & e)
  {
     array<Vec3,numFuns> dofPts =
     {
       e.pointList[0]->coord,
       e.pointList[1]->coord,
       e.pointList[2]->coord,
       e.pointList[3]->coord,
       0.5*(e.pointList[0]->coord+e.pointList[1]->coord),
       0.5*(e.pointList[1]->coord+e.pointList[2]->coord),
       0.5*(e.pointList[2]->coord+e.pointList[3]->coord),
       0.5*(e.pointList[3]->coord+e.pointList[0]->coord),
       e.midpoint()
     };
     return std::move(dofPts);
  }
};
