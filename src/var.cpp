#include "var.hpp"

std::vector<uint> offsetInit(std::vector<uint> blocks)
{
  std::vector<uint> offsets(blocks.size(), 0U);
  std::partial_sum(blocks.begin(), blocks.end() - 1, offsets.begin() + 1);
  return offsets;
}

