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
      }
    }
    vec = bool_array::Constant(feSpace.dof.totalNum, false);
  }

  void init()
  {
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
std::ostream & operator<<(std::ostream & out, bc_ess<FESpace> const & bc)
{
  out << "bc on marker " << bc.marker << "\n";
  out << "bool vec: ";
  for(uint i=0; i<bc.vec.size(); i++)
  {
    out << bc.vec(i) << " ";
  }
  out << "\n";
  out << "point_set: ";
  for(auto & i: bc.point_set)
  {
    out << i << " ";
  }
  out << "\n";
  return out;
}

template <typename FESpace>
class bc_list: public std::vector<bc_ess<FESpace>>
{
public:
  bc_list(FESpace const & feSpace, std::initializer_list<bc_ess<FESpace>> list):
    std::vector<bc_ess<FESpace>>(list)
  {
    vec = bool_array::Constant(feSpace.dof.totalNum, false);
  }

  void init()
  {
    for(auto& bc: *this)
    {
      bc.init();
    }
    for(uint i=0; i<vec.size(); ++i)
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

template <typename FESpace>
std::ostream & operator<<(std::ostream & out, bc_list<FESpace> const & bclist)
{
  out << "bc list\n";
  for(auto & bc: bclist)
  {
    out << bc << "\n";
  }
  out << "bool vec: ";
  for(uint i=0; i<bclist.vec.size(); i++)
  {
    out << bclist.vec(i) << " ";
  }
  out << "\n";

  return out;
}
