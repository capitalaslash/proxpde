#pragma once

#include "def.hpp"

#include "builder.hpp"

namespace proxpde
{

template <auto Storage = StorageType::ClmMajor>
struct EqnSolver
{};

template <typename FESpace>
struct Diagonal;
template <typename FESpace>
struct AssemblyVector;

template <typename LhsTup, typename RhsTup, auto Storage = StorageType::ClmMajor>
struct Eqn
{
  using FESpace_T = typename std::tuple_element_t<0, LhsTup>::FESpace_T;
  using BCList_T = std::vector<BCEss<FESpace_T>>;
  static uint const dim = FESpace_T::dim;

  Eqn(Var & s,
      LhsTup const & lhs,
      RhsTup const & rhs,
      BCList_T const & bcs = BCList_T{},
      EqnSolver<Storage> /*solver*/ = EqnSolver<Storage>{}):
      lhsAssemblies{lhs},
      rhsAssemblies{rhs},
      bcList{bcs},
      sol{s},
      builder{static_cast<uint>(s.data.size())}
  {}

  void buildLhs()
  {
    builder.buildLhs(lhsAssemblies, bcList);
    builder.closeMatrix();
  }

  void buildRhs() { builder.buildRhs(rhsAssemblies, bcList); }

  void build()
  {
    buildLhs();
    buildRhs();
  }

  void compute() { solver.compute(builder.A); }

  void solve() { sol.data = solver.solve(builder.b); }

  double residualNorm() { return (builder.A * sol.data - builder.b).norm(); }

  LhsTup const & lhsAssemblies;
  RhsTup const & rhsAssemblies;
  BCList_T const & bcList;
  Var & sol;
  Builder<Storage> builder;
  RecommendedSolverType<Storage> solver;
};

} // namespace proxpde
