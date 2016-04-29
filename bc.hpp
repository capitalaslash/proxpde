#pragma once

#include "def.hpp"
#include "curfe.hpp"

#include <unordered_set>
#include <list>

typedef Eigen::Array<bool,Eigen::Dynamic,1> BoolArray_T;

// TODO: add static map to set components via flags

template <typename FESpace>
class BCEss
{
public:
  BCEss(FESpace const& feSpace,
        marker_T m,
        scalarFun_T v,
        std::array<uint,FESpace::dim> comp = {0}):
    marker(m),
    value(v)
  {
    for(auto& f: feSpace.meshPtr->facetList)
    {
      if(f.marker == marker)
      {
        // TODO: assert that facingElem[1] is null, so the facet is on the boundary
        id_T const iElemId = f.facingElem[0].first->id;
        uint side = f.facingElem[0].second;
        for(uint i=0; i<FESpace::RefFE_T::dofPerFacet; i++)
        {
          DOFid_T const dof =
              feSpace.dof.elemMap[iElemId][FESpace::RefFE_T::dofOnFacet[side][i]];
          for(uint d=0; d<FESpace::dim; ++d)
          {
            bool const compIsConstrained =
                std::any_of(
                  std::begin(comp),
                  std::end(comp),
                  [&d](uint i) { return i == d; }
                );
            if(compIsConstrained)
            {
              constrainedDOFset.insert(dof + d*feSpace.dof.totalNum);
            }
          }
        }
      }
    }

    boolVector = BoolArray_T::Constant(feSpace.dof.totalNum*FESpace::dim, false);
    for(auto& id: constrainedDOFset)
    {
      boolVector[id] = true;
    }
  }

  bool isConstrained(DOFid_T const& id) const
  {
    return boolVector[id];
  }

  marker_T marker;
  scalarFun_T value;
  BoolArray_T boolVector;
  std::unordered_set<DOFid_T> constrainedDOFset;
};

template <typename FESpace>
std::ostream & operator<<(std::ostream & out, BCEss<FESpace> const & bc)
{
  out << "bc on marker " << bc.marker << "\n";
  out << "boolVector: ";
  for(uint i=0; i<bc.boolVector.size(); i++)
  {
    out << bc.boolVector(i) << " ";
  }
  out << "\n";
  out << "constrainedDOFset: ";
  for(auto & i: bc.constrainedDOFset)
  {
    out << i << " ";
  }
  out << "\n";
  return out;
}

template <typename FESpace>
struct BCNat
{
  typedef typename FESpace::RefFE_T::RefFacet_T Elem_T;
  typedef typename SideQR<typename FESpace::QR_T>::Type QR_T;
  typedef CurFE<Elem_T, QR_T> CurFE_T;
  typedef typename CurFE_T::RefFE_T RefFE_T;

  explicit BCNat(marker_T m, scalarFun_T const & f):
    marker(m),
    value(f)
  {
    // TODO: create a list of constrained faces at the beginning?
  }

  CurFE_T curFE;
  marker_T marker;
  scalarFun_T value;
};

template <typename FESpace>
class BCList
{
public:
  explicit BCList(FESpace const & fe):
  feSpace(fe)
  {}

  void addEssentialBC(
      marker_T const m,
      scalarFun_T const & f,
      std::array<uint,FESpace::dim> comp = {0})
  {
    bcEss_list.emplace_back(feSpace, m, f, comp);
  }

  void addNaturalBC(marker_T const m, scalarFun_T const & f)
  {
    bcNat_list.emplace_back(m, f);
  }

  FESpace const & feSpace;
  std::list<BCEss<FESpace>> bcEss_list;
  std::list<BCNat<FESpace>> bcNat_list;
};

template <typename FESpace>
std::ostream & operator<<(std::ostream & out, BCList<FESpace> const & bclist)
{
  out << "bc list\n";
  for(auto & bc: bclist.bcEss_list)
  {
    out << bc << "\n";
  }
  return out;
}
