#pragma once

#include "def.hpp"
#include "dof.hpp"
#include "curfe.hpp"

template <typename Mesh,
          typename RefFE,
          typename QR,
          uint d = 1>
struct FESpace
{
  typedef Mesh Mesh_T;
  typedef RefFE RefFE_T;
  typedef QR QR_T;
  typedef DOF<Mesh, RefFE, d> DOF_T;
  typedef CurFE<RefFE, QR> CurFE_T;
  static uint const dim = d;

  explicit FESpace(std::shared_ptr<Mesh> const mesh):
    meshPtr(mesh),
    dof(*mesh)
  {}

  std::shared_ptr<Mesh> const meshPtr;
  CurFE_T curFE;
  DOF_T dof;
};

template <typename ... FESpaces>
struct FESpaceList
{
  explicit FESpaceList(FESpaces ... feSpaces):
    list{feSpaces ...}
  {}

  std::tuple<FESpaces ...> list;
};

template <typename FESpace>
void interpolateAnalyticFunction(scalarFun_T const & f,
                                   FESpace const & feSpace,
                                   Vec & v)
{
  for(auto const & e: feSpace.meshPtr->elementList)
  {
    uint i = 0;
    for(auto const & dof: feSpace.dof.elemMap[e.id])
    {
      auto const pt = FESpace::RefFE_T::dofPts(e)[i];
      v[dof] = f(pt);
      i++;
    }
  }
}
