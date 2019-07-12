#pragma once

#include "def.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "builder.hpp"

template <typename FESpace>
struct Diagonal;
template <typename FESpace>
struct AssemblyVector;

template <typename FESpace, typename BCS, StorageType Storage = StorageType::ClmMajor>
struct Eqn
{
  using FESpace_T = FESpace;
  static uint const dim = FESpace_T::dim;

  Eqn(std::string_view const name, FESpace_T const & fe, BCS const & bcs):
    feSpace{fe},
    bcList{bcs},
    sol{name, std::vector<uint>(dim, feSpace.dof.size)},
    builder{feSpace.dof.size*dim}
  {}

  void buildLhs()
  {
    for (auto & assembly: assemblyListLhs)
    {
      builder.buildLhs(*assembly, bcList);
    }
    builder.closeMatrix();
  }

  void buildRhs()
  {
    for (auto & assembly: assemblyListRhs)
    {
      builder.buildRhs(*assembly, bcList);
    }
  }

  void build()
  {
    buildLhs();
    buildRhs();
  }

  void compute()
  {
    solver.compute(builder.A);
  }

  void solve()
  {
    sol.data = solver.solve(builder.b);
  }

  double residualNorm()
  {
    return (builder.A * sol.data - builder.b).norm();
  }

  FESpace_T const & feSpace;
  BCS const & bcList;
  BlockVar sol;
  Builder<Storage> builder;
  std::vector<std::unique_ptr<Diagonal<FESpace_T>>> assemblyListLhs;
  std::vector<std::unique_ptr<AssemblyVector<FESpace_T>>> assemblyListRhs;
  RecommendedSolverType<Storage> solver;
};
