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

template <typename LhsTup,
          typename RhsTup,
          typename BCS,
          StorageType Storage = StorageType::ClmMajor>
struct Eqn
{
  using FESpace_T = typename std::tuple_element_t<0, LhsTup>::FESpace_T;
  static uint const dim = FESpace_T::dim;

  Eqn(Var & s,
      LhsTup const & lhs,
      RhsTup const & rhs,
      BCS const & bcs):
    lhsAssemblies{lhs},
    rhsAssemblies{rhs},
    bcList{bcs},
    sol{s},
    builder{s.data.size()}
  {}

  void buildLhs()
  {
    builder.buildLhs(lhsAssemblies, bcList);
    builder.closeMatrix();
  }

  void buildRhs()
  {
    static_for(rhsAssemblies, [&] (auto const & /*i*/, auto & assembly)
    {
      builder.buildRhs(assembly, bcList);
    });
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

  LhsTup const & lhsAssemblies;
  RhsTup const & rhsAssemblies;
  BCS const & bcList;
  Var & sol;
  Builder<Storage> builder;
  RecommendedSolverType<Storage> solver;
};
