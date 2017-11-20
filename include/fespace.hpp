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
  using Mesh_T = Mesh;
  using RefFE_T = RefFE;
  using QR_T = QR;
  using DOF_T = DOF<Mesh, RefFE, d>;
  using CurFE_T = CurFE<RefFE, QR>;
  static uint const dim = d;

  explicit FESpace(std::shared_ptr<Mesh> const mesh):
    meshPtr(mesh),
    dof(*mesh)
  {}

  std::shared_ptr<Mesh> const meshPtr;
  CurFE_T curFE;
  DOF_T dof;
};

template <typename FESpace>
void interpolateAnalyticFunction(scalarFun_T const & f,
                                 FESpace const & feSpace,
                                 Vec & v,
                                 uint const offset = 0)
{
  for(auto const & e: feSpace.meshPtr->elementList)
  {
    uint i = 0;
    for(auto const & dof: feSpace.dof.elemMap[e.id])
    {
      auto const pt = FESpace::RefFE_T::dofPts(e)[i];
      v[offset + dof] = f(pt);
      i++;
    }
  }
}
