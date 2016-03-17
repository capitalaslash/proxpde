#pragma once

#include "def.hpp"

template <typename Mesh, typename RefFE>
struct DOF
{
  explicit DOF(Mesh const& mesh)
  {
    std::cout << "_clms = " << _clms << std::endl;
    _map.resize(mesh.elementList.size());
  }

  static uint constexpr _clms = numDOFs<RefFE>();
  std::vector<std::array<id_T,_clms>> _map;
};
