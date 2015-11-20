#pragma once

#include <set>

#include "def.hpp"

typedef Eigen::Array<bool,Eigen::Dynamic,1> bool_array;

class bc_ess: public std::set<id_T>
{
public:
  bc_ess(Mesh1D const& mesh, marker_T m, scalarFun_T v):
    marker(m),
    value(v)
  {
    for(auto& p: mesh.pointList)
    {
      if(p.marker == marker)
      {
        this->insert(p.id);
        std::cout << "inserting " << p.id << " " << p.marker << " " << marker << std::endl;
      }
    }
  }

  void init(uint const numPts)
  {
    vec = bool_array::Constant(numPts, false);
    for(auto& id: *this)
    {
      vec[id] = true;
    }
  }

  bool is_constrained(Point const& p) const
  {
    return vec[p.id];
  }

  marker_T marker;
  scalarFun_T value;
  bool_array vec;
};

class bc_list: public std::vector<bc_ess>
{
public:
  bc_list(std::initializer_list<bc_ess> list):
    std::vector<bc_ess>(list)
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
