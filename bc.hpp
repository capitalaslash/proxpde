#pragma once

#include <set>

#include "def.hpp"

typedef Eigen::Array<bool,Eigen::Dynamic,1> bool_array;

class bc_set: public std::set<id_T>
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

  bool_array vec;
  scalarFun_T value;
};
