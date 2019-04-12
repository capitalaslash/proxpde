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

template <typename FESpace>
struct Eqn
{
  using FESpace_T = FESpace;
  static uint const dim = FESpace_T::dim;

  Eqn(std::string_view const name, FESpace_T const & fe):
    feSpace{fe},
    bcList{feSpace},
    sol{name, std::vector<uint>(dim, feSpace.dof.size)},
    builder{feSpace.dof.size*dim}
  {}

  void buildLhs()
  {
    for (auto & assembly: assemblyListLhs)
    {
      builder.buildProblem(*assembly, bcList);
    }
    builder.closeMatrix();
  }

  void buildRhs()
  {
    for (auto & assembly: assemblyListRhs)
    {
      builder.buildProblem(*assembly, bcList);
    }
  }

  void build()
  {
    buildLhs();
    buildRhs();
  }

  void solve()
  {
    solver.compute(builder.A);
    sol.data = solver.solve(builder.b);
  }

  double residualNorm()
  {
    return (builder.A * sol.data - builder.b).norm();
  }

  FESpace_T feSpace;
  BCList<FESpace_T> bcList;
  BlockVar sol;
  Builder builder;
  std::vector<std::unique_ptr<Diagonal<FESpace_T>>> assemblyListLhs;
  std::vector<std::unique_ptr<AssemblyVector<FESpace_T>>> assemblyListRhs;
  LUSolver solver;
};
