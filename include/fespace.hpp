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
void interpolateAnalyticFunction(Fun<FESpace::dim,3> const & f,
                                 FESpace const & feSpace,
                                 Vec & v,
                                 uint const offset = 0)
{
  for(auto const & e: feSpace.meshPtr->elementList)
  {
    uint p = 0;
    for(auto const & dof: feSpace.dof.elemMap[e.id])
    {
      auto const d = p / FESpace::RefFE_T::numFuns;
      auto const pt = FESpace::RefFE_T::dofPts(e)[p % FESpace::RefFE_T::numFuns];
      v[offset + dof] = f(pt)[d];
      p++;
    }
  }
}

template <typename FESpace>
void interpolateAnalyticFunction(scalarFun_T const & f,
                                 FESpace const & feSpace,
                                 Vec & v,
                                 uint const offset = 0)
{
  interpolateAnalyticFunction([f](Vec3 const &p){return Vec1(f(p));}, feSpace, v, offset);
}
