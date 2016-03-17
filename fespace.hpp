#pragma once

#include "def.hpp"
#include "dof.hpp"
#include "curfe.hpp"

template <typename Mesh,
          typename RefFE,
          typename QR>
struct FESpace
{
  typedef Mesh Mesh_T;
  typedef RefFE RefFE_T;
  typedef QR QR_T;
  typedef DOF<Mesh, RefFE> DOF_T;
  typedef CurFE<RefFE,QR> CurFE_T;

  explicit FESpace(std::shared_ptr<Mesh> const mesh):
    meshPtr(mesh),
    dof(*mesh)
  {}

  std::shared_ptr<Mesh> const meshPtr;
  CurFE_T curFE;
  DOF_T dof;
};
