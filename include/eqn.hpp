#pragma once

#include "def.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "builder.hpp"
#include "assembly.hpp"

template <
    typename Mesh,
    typename RefFE,
    typename QR,
    uint Dimension = 1>
struct Eqn
{
  using Mesh_T = Mesh;
  using RefFE_T = RefFE;
  using QR_T = QR;
  using FESpace_T = FESpace<Mesh_T, RefFE_T, QR_T, Dimension>;
  static uint const dim = Dimension;

  Eqn(std::string_view const name, Mesh const & mesh):
    feSpace{mesh},
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

