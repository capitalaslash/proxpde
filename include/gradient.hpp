#pragma once

#include "def.hpp"
#include "fe.hpp"
#include "builder.hpp"

template <typename FESpaceT>
using Grad_T =
  FESpace<typename FESpaceT::Mesh_T,
          typename FEType<typename FESpaceT::Mesh_T::Elem_T, order_v<typename FESpaceT::RefFE_T>-1>::RefFE_T,
          typename FESpaceT::QR_T,
          FESpaceT::dim * FESpaceT::Mesh_T::Elem_T::dim>;

template <typename FESpaceGrad, typename FESpaceOrig, typename Solver = LUSolver>
void computeGradient(
        Vec & grad, FESpaceGrad & feSpaceGrad,
        Vec const & u, FESpaceOrig const & feSpaceOrig)
{
  static_assert(std::is_same_v<Grad_T<FESpaceOrig>, FESpaceGrad>);

  std::tuple<> bcsGrad;
  Builder builderGrad{feSpaceGrad.dof.size * FESpaceGrad::dim};
  builderGrad.buildLhs(std::tuple{AssemblyMass{1.0, feSpaceGrad}}, bcsGrad);
  builderGrad.buildRhs(std::tuple{AssemblyGradRhs{1.0, u, feSpaceOrig, feSpaceGrad}}, bcsGrad);
  builderGrad.closeMatrix();
  Solver solverGrad;
  solverGrad.compute(builderGrad.A);
  grad = solverGrad.solve(builderGrad.b);
}
