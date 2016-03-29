#pragma once

#include "def.hpp"
#include "fespace.hpp"

struct Var
{
  explicit Var(std::string n, uint size = 0):
    name(n),
    data(size)
  {}

  std::string name;
  Vec data;
};
