#pragma once

#include <set>
#include <list>

#include "def.hpp"
#include "curfe.hpp"

typedef Eigen::Array<bool,Eigen::Dynamic,1> bool_array;

template <typename FESpace>
class BCEss
{
public:
  BCEss(FESpace const& feSpace, marker_T m, scalarFun_T v):
    marker(m),
    value(v)
  {
    for(auto& f: feSpace.meshPtr->facetList)
    {
      if(f.marker == marker)
      {
        id_T const iElemId = f.facingElem[0].first->id;
        uint side = f.facingElem[0].second;
        for(uint i=0; i<FESpace::RefFE_T::dofPerFacet; i++)
        {
          DOFid_T const dof =
              feSpace.dof.elemMap[iElemId][FESpace::RefFE_T::dofOnFacet[side][i]];
          point_set.insert(dof);
        }
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
std::ostream & operator<<(std::ostream & out, BCEss<FESpace> const & bc)
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
struct BCNat
{
  typedef CurFE<
    typename FESpace::RefFE_T::RefFacet_T,
    typename FESpace::QR_T>
      CurFE_T;

  explicit BCNat(marker_T m, scalarFun_T const & f):
    marker(m),
    value(f)
  {}

  CurFE_T curFE;
  marker_T marker;
  scalarFun_T value;
};

template <typename FESpace>
class BCList: public std::vector<BCEss<FESpace>>
{
public:
  explicit BCList(FESpace const & feSpace, std::initializer_list<BCEss<FESpace>> list = {}):
    std::vector<BCEss<FESpace>>(list)
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

  void addBCNat(marker_T m, scalarFun_T const & v)
  {
    bcNat_list.emplace_back(m, v);
  }

  bool_array vec;
  std::list<BCNat<FESpace>> bcNat_list;
};

template <typename FESpace>
std::ostream & operator<<(std::ostream & out, BCList<FESpace> const & bclist)
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
