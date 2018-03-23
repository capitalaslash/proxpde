#pragma once

#include "def.hpp"
#include "fespace.hpp"

struct Var
{
  explicit Var(std::string n, uint size = 0):
    name(n),
    data(size)
  {}

  Var(std::string n, Vec const & vec):
    name(n),
    data(vec)
  {}

  Var(std::string n, Vec const & vec, uint offset, uint size):
    name(n),
    data(vec.block(offset,0,size,1))
  {}

  template <typename FESpace>
  Var(std::string n, FESpace const & feSpace):
    name(n),
    data(feSpace.dof.totalNum * feSpace.dim)
  {}

  std::string name;
  Vec data;
};
