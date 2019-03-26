#pragma once

#include "def.hpp"
#include "curfe.hpp"

#include <unordered_set>
#include <list>

// TODO: add static map to set components via flags

using DofSet_T = std::unordered_set<DOFid_T>;
using BoolArray_T = Eigen::Array<bool,Eigen::Dynamic,1>;

template <typename FESpace>
DofSet_T fillDofSet(
    FESpace const & feSpace,
    marker_T marker,
    std::vector<uint> const & comp)
{
  using RefFE_T = typename FESpace::RefFE_T;
  DofSet_T constrainedDOFset;
  for(auto const & f: feSpace.mesh.facetList)
  {
    if(f.marker == marker)
    {
      // TODO: assert that facingElem[1] is null, so the facet is on the boundary
      auto const & [elem, side] = f.facingElem[0];
      for(uint i=0; i<RefFE_T::dofPerFacet; i++)
      {
        DOFid_T const dof = feSpace.dof.getId(elem->id, RefFE_T::dofOnFacet[side][i]);
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
            constrainedDOFset.insert(dof + d*feSpace.dof.size);
          }
        }
      }
      // std::cout << "current dof set: ";
      // for (auto const & id: constrainedDOFset)
      // {
      //     std::cout << id << " ";
      // }
      // std::cout << std::endl;
    }
  }
  return constrainedDOFset;
}

template <typename FESpace>
class BCEss
{
public:
  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace_T::CurFE_T;

  // use manual-provided set of DOFs
  BCEss(FESpace_T const & feSpace,
        DofSet_T const dofSet,
        Fun<FESpace::dim,3> const v):
    curFE(feSpace.curFE),
    constrainedDOFSet(std::move(dofSet)),
    constrainedDOFVector(BoolArray_T::Constant(feSpace.dof.size*FESpace_T::dim, false)),
    dimSize(feSpace.dof.size),
    value(std::move(v))
  {
    fillConstrainedDOFVector();
    std::cout << "new bc on dofset with " << constrainedDOFSet.size() << " dofs" << std::endl;
  }

  // this can be used with scalar FESpaces only
  BCEss(FESpace_T const & feSpace,
        DofSet_T const dofSet,
        scalarFun_T const & v):
    BCEss(feSpace, std::move(dofSet), [v](Vec3 const & p){return Vec1{v(p)};})
  {}

  BCEss(FESpace_T const & feSpace,
        marker_T const m,
        Fun<FESpace_T::dim,3> const v,
        std::vector<uint> const & comp = allComp<FESpace>()):
    curFE(feSpace.curFE),
    constrainedDOFSet(fillDofSet(feSpace, m, comp)),
    constrainedDOFVector(BoolArray_T::Constant(feSpace.dof.size*FESpace_T::dim, false)),
    dimSize(feSpace.dof.size),
    value(std::move(v)),
    marker(m)
  {
    fillConstrainedDOFVector();
    std::cout << "new bc on marker " << m << " with " << constrainedDOFSet.size() << " dofs" << std::endl;
  }

  // this can be used with scalar FESpaces only
  BCEss(FESpace_T const & feSpace,
        marker_T m,
        scalarFun_T const & v,
        std::vector<uint> const & comp = allComp<FESpace>()):
    curFE(feSpace.curFE),
    constrainedDOFSet(fillDofSet(feSpace, m, comp)),
    constrainedDOFVector(BoolArray_T::Constant(feSpace.dof.size*FESpace_T::dim, false)),
    dimSize(feSpace.dof.size),
    value([v](Vec3 const & p){return Vec1{v(p)};}),
    marker(m)
  {
    fillConstrainedDOFVector();
    std::cout << "new bc on marker " << m << " with " << constrainedDOFSet.size() << " dofs" << std::endl;
  }

  bool isConstrained(DOFid_T const id, int const d = 0) const
  {
    return constrainedDOFVector[id + d*dimSize];
  }

  FVec<FESpace_T::dim> evaluate(DOFid_T const & i) const
  {
    if constexpr (Family<typename FESpace_T::RefFE_T>::value == FamilyType::LAGRANGE)
    {
      // TODO: this is a crude implementation that works only for Lagrange elements
      return value(curFE.dofPts[i]);
    }
    else if constexpr (Family<typename FESpace_T::RefFE_T>::value == FamilyType::RAVIART_THOMAS)
    {
      std::abort();
      //return value(curFE.dofPts[i]).dot(normal);
    }
    else
    {
      std::abort();
    }
  }

  friend std::ostream & operator<<(std::ostream & out, BCEss<FESpace_T> const & bc)
  {
    out << "bc constrainedDOFVector: " << bc.constrainedDOFVector.transpose() << std::endl;
    out << "constrainedDOFset: ";
    for(auto const i: bc.constrainedDOFSet)
    {
      out << i << " ";
    }
    out << std::endl;
    return out;
  }

protected:
  void fillConstrainedDOFVector()
  {
    for(auto const id: constrainedDOFSet)
    {
      constrainedDOFVector[id] = true;
    }
  }

  CurFE_T const & curFE;
  DofSet_T constrainedDOFSet;
  BoolArray_T constrainedDOFVector;
  uint const dimSize;
  Fun<FESpace::dim,3> const value;

public:
  marker_T marker = markerNotSet;
  double diag = 1.0;
};

template <typename FESpace>
struct BCNat
{
  using FESpace_T = FESpace;
  using Elem_T = typename FESpace::RefFE_T::RefFacet_T;
  using QR_T = SideQR_T<typename FESpace::QR_T>;
  using CurFE_T = CurFE<Elem_T, QR_T>;
  using RefFE_T = typename CurFE_T::RefFE_T;

  explicit BCNat(
      marker_T const m,
      Fun<FESpace::dim,3> const f,
      std::vector<uint> const c = allComp<FESpace>()):
    marker(m),
    value(std::move(f)),
    comp(std::move(c))
  {
    // TODO: create a list of constrained faces at the beginning?
  }

  bool hasComp(uint c)
  {
    return std::find(comp.begin(), comp.end(), c) != comp.end();
  }

  CurFE_T curFE;
  marker_T marker;
  Fun<FESpace::dim,3> const value;
  std::vector<uint> const comp;
};

// this class deals only with the matrix part of the mixed bc, the rhs part is
// managed as a standard natural bc
template <typename FESpace>
struct BCMixed
{
  using Elem_T = typename FESpace::RefFE_T::RefFacet_T;
  using QR_T = SideQR_T<typename FESpace::QR_T>;
  using CurFE_T = CurFE<Elem_T, QR_T>;
  using RefFE_T = typename CurFE_T::RefFE_T;

  explicit BCMixed(marker_T m, Fun<FESpace::dim,3> const f, std::vector<uint> const c = allComp<FESpace>()):
    marker(m),
    coeff(std::move(f)),
    comp(std::move(c))
  {
    // TODO: create a list of constrained faces at the beginning?
  }

  bool hasComp(uint c)
  {
    return std::find(comp.begin(), comp.end(), c) != comp.end();
  }

  CurFE_T curFE;
  marker_T marker;
  Fun<FESpace::dim,3> const coeff;
  std::vector<uint> const comp;
};

template <typename FESpace>
class BCList
{
public:
  explicit BCList(FESpace const & fe):
    feSpace(fe)
  {}

  template <typename BC>
  void addBC(BC const & bc)
  {
    static_assert(std::is_same_v<typename BC::FESpace_T, FESpace>, "this BC does not use the same FESpace");
    if (bc.marker != markerNotSet)
    {
      if (checkMarkerFixed(bc.marker))
      {
        std::cerr << "the marker " << bc.marker << " has already been fixed." << std::endl;
        abort();
      }
      fixedMarkers.insert(bc.marker);
    }

    if constexpr (std::is_same_v<BC, BCEss<FESpace>>)
    {
      bcEssList.push_back(bc);
    }
    else if constexpr (std::is_same_v<BC, BCNat<FESpace>>)
    {
      bcNatList.push_back(bc);
    }
    else if constexpr (std::is_same_v<BC, BCMixed<FESpace>>)
    {
      bcNatList.emplace_back(bc.marker, bc.b, bc.comp);
      bcMixedList.push_back(bc);
    }
    else
    {
      // this should never happen
      std::cerr << "a BC of type " << typeid(bc).name() << " has been added." << std::endl;
      std::abort();
    }
  }

  void addNaturalBC(
      marker_T const m,
      Fun<FESpace::dim,3> const f,
      std::vector<uint> const comp = allComp<FESpace>())
  {
    if (checkMarkerFixed(m))
    {
      std::cerr << "the marker " << m << " has already been fixed." << std::endl;
      abort();
    }
    bcNatList.emplace_back(m, std::move(f), std::move(comp));
    fixedMarkers.insert(m);
  }

  // this overload works only when FESpace::dim == 1
  template <typename FESpace1 = FESpace,
            std::enable_if_t<FESpace1::dim == 1, bool> = true>
  void addNaturalBC(
      marker_T const m,
      scalarFun_T const f,
      std::vector<uint> const comp = allComp<FESpace>())
  {
    addNaturalBC(m, [f] (Vec3 const & p) {return Vec1(f(p));}, std::move(comp));
  }

  // otherwise, we should fail
  template <typename FESpace1 = FESpace,
            std::enable_if_t<FESpace1::dim != 1, bool> = true>
  void addNaturalBC(
      marker_T const,
      scalarFun_T const,
      std::vector<uint> const)
  {
    std::cerr << __FILE__ << ":" << __LINE__ << " > this funcytion only works with FESpace::dim == 1" << std::endl;
    std::abort();
  }

  // mixed BC: a * u + \nabla u = b
  // - \lap u = f
  // (\nabla u, \nabla v) - <\nabla u, v> = (f, v)
  // (., .) == integration on volume
  // <., .> == integration on boundary
  // (\nabla u, \nabla v) + <a*u, v> = (f, v) + <b, v>
  // a >= \eps > 0 cohercitivity to guarantee solution
  void addMixedBC(
      marker_T const m,
      scalarFun_T const a,
      scalarFun_T const b,
      std::vector<uint> const comp = allComp<FESpace>())
  {
    addNaturalBC(m, std::move(b), comp);
    bcMixedList.emplace_back(m, [a] (Vec3 const & p) {return Vec1(a(p));}, std::move(comp));
  }

  bool checkMarkerFixed(marker_T const m) const
  {
    return fixedMarkers.find(m) != fixedMarkers.end();
  }

  // bool checkDofFixed(id_T const id) const
  // {
  //   return fixedDofs.find(id) != fixedDofs.end();
  // }

  FESpace const & feSpace;
  std::unordered_set<marker_T> fixedMarkers;
  // std::unordered_set<id_T> fixedDofs;
  std::list<BCEss<FESpace>> bcEssList;
  std::list<BCNat<FESpace>> bcNatList;
  std::list<BCMixed<FESpace>> bcMixedList;
};

template <typename FESpace>
std::ostream & operator<<(std::ostream & out, BCList<FESpace> const & bcList)
{
  out << "bc list with " << bcList.bcEssList.size() << " essential bcs and "
      << bcList.bcNatList.size() << " natural bcs" << std::endl;
  for(auto & bc: bcList.bcEssList)
  {
    out << bc << "\n";
  }
  return out;
}

template <typename FESpace>
class DOFCoordSet
{
public:
  using FESpace_T = FESpace;
  using Predicate_T = std::function<bool (Vec3 const &)>;

  DOFCoordSet(FESpace_T & fe, Predicate_T const & p):
    feSpace(fe),
    predicate(p)
  {
    for (auto const & e: feSpace.mesh.elementList)
    {
      feSpace.curFE.reinit(e);
      for (uint d=0; d<numDOFs<typename FESpace_T::RefFE_T>(); ++d)
      {
        if (predicate(feSpace.curFE.dofPts[d]))
        {
          ids.insert(feSpace.dof.getId(e.id, d));
        }
      }
    }
    std::cout << "new DOFCoordSet with " << ids.size() << " dofs" << std::endl;
  }

  FESpace_T & feSpace;
  Predicate_T const & predicate;
  DofSet_T ids;
};
