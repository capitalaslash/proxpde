#pragma once

#include <set>

#include "def.hpp"

typedef Eigen::Array<bool,Eigen::Dynamic,1> bool_array;

template <typename FESpace>
class bc_ess
{
public:
  bc_ess(FESpace const& feSpace, marker_T m, scalarFun_T v):
    marker(m),
    value(v)
  {
    for(auto& p: feSpace.meshPtr->pointList)
    {
      if(p.marker == marker)
      {
        point_set.insert(p.dof_ids[0]);
        // std::cout << "inserting " << p.id << " " << p.marker << " " << marker << std::endl;
      }
    }
  }

  void init(uint const numPts)
  {
    vec = bool_array::Constant(numPts, false);
    for(auto& id: point_set)
    {
      vec[id] = true;
    }
  }

  bool is_constrained(DOFid_T const& id) const
  {
    return vec[id];
  }

  marker_T marker;
  scalarFun_T value;
  bool_array vec;
  std::set<DOFid_T> point_set;
};

template <typename FESpace>
class bc_list: public std::vector<bc_ess<FESpace>>
{
public:
  bc_list(std::initializer_list<bc_ess<FESpace>> list):
    std::vector<bc_ess<FESpace>>(list)
    {}

  void init(uint const numPts)
  {
    for(auto& bc: *this)
    {
      bc.init(numPts);
    }
    vec = bool_array::Constant(numPts, false);
    for(uint i=0; i<numPts; ++i)
    {
      for(auto& bc: *this)
      {
        vec[i] = vec[i] || bc.vec[i];
      }
    }
  }

  bool is_constrained(Point const& p) const
  {
    return vec[p.id];
  }

  bool_array vec;
};
