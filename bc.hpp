#pragma once

#include <set>

#include "def.hpp"

typedef Eigen::Array<bool,Eigen::Dynamic,1> bool_array;

class bc_ess: public std::set<id_T>
{
public:
  void init(uint const numPts)
  {
    vec = bool_array::Constant(numPts, false);
    for(auto&id: *this)
    {
      vec[id] = true;
    }
  }

  bool is_constrained(Point const& p) const
  {
    return vec[p.id];
  }

  bool_array vec;
  scalarFun_T value;
};
