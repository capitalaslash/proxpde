#pragma once

#include "def.hpp"
#include "curfe.hpp"

#include <unordered_set>
#include <list>

// TODO: add static map to set components via flags

template <typename FESpace>
class BCEss
{
public:
  using DofSet_T = std::unordered_set<DOFid_T>;

  BCEss(FESpace const& feSpace,
        marker_T m,
        Fun<FESpace::dim,3> v,
        std::vector<uint> const & comp = {0}):
    marker(m),
    value(v),
    dimSize(feSpace.dof.totalNum)
  {
    fillDofSet(feSpace, comp);
    fillBoolVector(feSpace);
  }

  // use manual-provided set of DOFs for this condition
  BCEss(FESpace const& feSpace,
        DofSet_T const & dofSet,
        Fun<FESpace::dim,3> v,
        std::vector<uint> const & comp = {0}):
    marker(-1),
    value(v),
    constrainedDOFset(dofSet)
  {
    fillBoolVector(feSpace);
  }

  void fillDofSet(FESpace const& feSpace, std::vector<uint> const & comp)
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
              constrainedDOFset.insert(dof + d*dimSize);
            }
          }
        }
      }
    }
  }

  void fillBoolVector(FESpace const & feSpace)
  {
    boolVector = BoolArray_T::Constant(feSpace.dof.totalNum*FESpace::dim, false);
    for(auto& id: constrainedDOFset)
    {
      boolVector[id] = true;
    }
  }

  bool isConstrained(DOFid_T const& id, int const d = 0) const
  {
    return boolVector[id + d*dimSize];
  }

  double diag = 1.0;
  marker_T marker;
  DofSet_T constrainedDOFset;
  BoolArray_T boolVector;
  Fun<FESpace::dim,3> value;
  uint const dimSize;
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
  using Elem_T = typename FESpace::RefFE_T::RefFacet_T;
  using QR_T = typename SideQR<typename FESpace::QR_T>::Type;
  using CurFE_T = CurFE<Elem_T, QR_T>;
  using RefFE_T = typename CurFE_T::RefFE_T;

  explicit BCNat(marker_T m, Fun<FESpace::dim,3> const & f, std::vector<uint> c = {0}):
    marker(m),
    value(f),
    comp(c)
  {
    // TODO: create a list of constrained faces at the beginning?
  }

  bool hasComp(uint c)
  {
    return std::find(comp.begin(), comp.end(), c) != comp.end();
  }

  CurFE_T curFE;
  marker_T marker;
  Fun<FESpace::dim,3> value;
  std::vector<uint> comp;
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
      Fun<FESpace::dim,3> const & f,
      std::vector<uint> const & comp = {{0}})
  {
    if (checkMarkerFixed(m))
    {
      std::cerr << "the marker " << m << " has already been fixed." << std::endl;
      abort();
    }
    fixedMarkers.insert(m);
    bcEss_list.emplace_back(feSpace, m, f, comp);
  }

  void addEssentialBC(
      typename BCEss<FESpace>::DofSet_T dofSet,
      Fun<FESpace::dim,3> const & f,
      std::vector<uint> const & comp = {0})
  {
    bcEss_list.emplace_back(feSpace, dofSet, f, comp);
  }

  void addNaturalBC(
      marker_T const m,
      Fun<FESpace::dim,3> const & f,
      std::vector<uint> const & comp = {0})
  {
    if (checkMarkerFixed(m))
    {
      std::cerr << "the marker " << m << " has already been fixed." << std::endl;
      abort();
    }
    bcNat_list.emplace_back(m, f, comp);
    fixedMarkers.insert(m);
  }

  bool checkMarkerFixed(marker_T const m) const
  {
    return fixedMarkers.find(m) != fixedMarkers.end();
  }

  bool checkDofFixed(id_T const id) const
  {
    return fixedDofs.find(id) != fixedDofs.end();
  }

  FESpace const & feSpace;
  std::set<marker_T> fixedMarkers;
  std::set<id_T> fixedDofs;
  std::list<BCEss<FESpace>> bcEss_list;
  std::list<BCNat<FESpace>> bcNat_list;
};

template <typename FESpace>
std::ostream & operator<<(std::ostream & out, BCList<FESpace> const & bcList)
{
  out << "bc list with " << bcList.bcEss_list.size() << " essential bcs and "
      << bcList.bcNat_list.size() << " natural bcs" << std::endl;
  for(auto & bc: bcList.bcEss_list)
  {
    out << bc << "\n";
  }
  return out;
}
