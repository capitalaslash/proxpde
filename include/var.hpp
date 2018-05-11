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

  // template <typename FESpace>
  // Var(std::string n, FESpace const & feSpace):
  //   name(n),
  //   data(feSpace.dof.size * feSpace.dim)
  // {}

  std::string name;
  Vec data;
};

std::vector<uint> offsetInit(std::vector<uint> blocks)
{
  std::vector<uint> offsets(blocks.size(), 0U);
  std::partial_sum(blocks.begin(), blocks.end()-1, offsets.begin()+1);
  return offsets;
}

struct BlockVar: public Var
{
  BlockVar(std::string n, std::vector<uint> const & bs):
    Var{n, std::accumulate(bs.begin(), bs.end(), 0U)},
    offsets{offsetInit(bs)},
    blocks{bs}
  {}

  auto block(uint i)
  {
    return this->data.block(offsets[i], 0, blocks[i], 1);
  }

  std::vector<uint> const offsets;
  std::vector<uint> const blocks;
};
