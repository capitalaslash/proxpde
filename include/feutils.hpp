#pragma once

#include "def.hpp"
#include "fe.hpp"
#include "builder.hpp"

template <typename FESpaceTo, typename FESpaceFrom, typename Solver = LUSolver>
void l2Projection(
    Vec & to, FESpaceTo const & feSpaceTo,
    Vec const & from, FESpaceFrom const & feSpaceFrom)
{
  AssemblyMass massTo(1.0, feSpaceTo);
  AssemblyProjection projFromTo(1.0, from, feSpaceFrom, feSpaceTo);
  Builder builder{feSpaceTo.dof.size * FESpaceTo::dim};
  builder.buildLhs(std::tuple{massTo}, std::tuple{});
  builder.buildRhs(std::tuple{projFromTo}, std::tuple{});
  builder.closeMatrix();
  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  Solver solver(builder.A);
  to = solver.solve(builder.b);
}

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
  builderGrad.buildLhs(std::tuple{AssemblyScalarMass{1.0, feSpaceGrad}}, bcsGrad);
  builderGrad.buildRhs(std::tuple{AssemblyGradRhs{1.0, u, feSpaceOrig, feSpaceGrad}}, bcsGrad);
  builderGrad.closeMatrix();
  Solver solverGrad;
  solverGrad.compute(builderGrad.A);
  grad = solverGrad.solve(builderGrad.b);
}
