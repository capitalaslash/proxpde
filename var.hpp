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

  std::string name;
  Vec data;
};
