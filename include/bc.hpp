#pragma once

#include "def.hpp"

#include "qr.hpp"
#include "reffe.hpp"

template <typename RefFE, typename QR>
struct CurFE;

using DofSet_T = std::unordered_set<DOFid_T>;

// TODO: add static map to set components via flags

template <typename FESpace>
class BCEss
{
public:
  using DofMap_T = std::unordered_map<DOFid_T, id_T>;

  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace_T::CurFE_T;
  using RefFE_T = typename FESpace_T::RefFE_T;
  static constexpr uint dim = FESpace_T::dim;

  // use manual-provided set of DOFs
  // TODO: pass in a dofset anyway, the map is only needed internally
  BCEss(FESpace_T const & fe, DofMap_T const dofMap):
      feSpace(fe),
      _constrainedDofMap(std::move(dofMap)),
      data{Vec::Zero(_constrainedDofMap.size())}
  {
    std::cout << "new bc on dofset with " << _constrainedDofMap.size() << " dofs"
              << std::endl;
  }

  BCEss(
      FESpace_T const & fe,
      marker_T const m,
      std::vector<short_T> const & c = allComp<FESpace>()):
      feSpace(fe),
      marker(m),
      comp(c)
  {
    id_T counter = 0;
    for (auto const & facet: feSpace.mesh.facetList)
    {
      // TODO: store the list of facets with own marker
      if (facet.marker == marker)
      {
        facetIdList.push_back(facet.id);
        auto const & [elem, side] = facet.facingElem[0];
        if constexpr (
            order_v < RefFE_T >> 0 || family_v<RefFE_T> != FamilyType::LAGRANGE)
        {
          for (auto const dofFacet: RefFE_T::dofOnFacet[side])
          {
            for (auto const c: comp)
            {
              DOFid_T const dofId = feSpace.dof.getId(elem->id, dofFacet, c);
              [[maybe_unused]] auto const [ptr, inserted] =
                  _constrainedDofMap.insert({dofId, counter});
              if (inserted)
              {
                counter++;
              }
            }
          }
        }
        else // order_v<RefFE_T> == 0 && family_v<RefFE_T> == FamilyType::LAGRANGE
        {
          // P0 element do not have dofs on boundary, we set the element dof instead
          for (auto const c: comp)
          {
            DOFid_T const dofId = feSpace.dof.getId(elem->id, 0, c);
            // elements cannot be walked through multiple times, each insert() is
            // successful
            _constrainedDofMap.insert({dofId, counter});
            counter++;
          }
        }
      }
    }
    data = Vec::Zero(_constrainedDofMap.size());
    std::cout << "new bc on marker " << m << " with " << _constrainedDofMap.size()
              << " dofs" << std::endl;
  }

  BCEss<FESpace> & operator<<(Fun<FESpace::physicalDim(), 3> const & fun)
  {
    if (marker != markerNotSet)
    {
      for (auto const facetId: facetIdList)
      {
        auto const & facet = feSpace.mesh.facetList[facetId];
        auto const & [elem, side] = facet.facingElem[0];
        feSpace.curFE.reinit(*elem);

        for (auto const dofFacet: RefFE_T::dofOnFacet[side])
        {
          auto const value = fun(feSpace.curFE.dofPts[dofFacet]);
          if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
          {
            for (auto const c: comp)
            {
              DOFid_T const dofId = feSpace.dof.getId(elem->id, dofFacet, c);
              data[_constrainedDofMap.at(dofId)] = value[c];
            }
          }
          else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
          {
            DOFid_T const dofId = feSpace.dof.getId(elem->id, dofFacet);
            data[_constrainedDofMap.at(dofId)] =
                promote<3>(value).dot(facet.normal()) * facet.volume();
          }
        }
      }
    }
    else // no marker, work on the dof set
    {
      auto const size = numDOFs<RefFE_T>();
      for (auto const & facet: feSpace.mesh.facetList)
      {
        auto const & elem = *(facet.facingElem[0].ptr);
        for (uint d = 0; d < size * FESpace_T::dim; ++d)
        {
          auto const dof = feSpace.dof.elemMap(elem.id, d);
          if (_constrainedDofMap.count(dof) > 1)
          {
            feSpace.curFE.reinit(elem);
            auto const localDof = d % size;
            auto const c = d / size;
            auto const value = fun(feSpace.curFE.dofPts[localDof]);
            if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
            {
              data[_constrainedDofMap.at(dof)] = value[c];
            }
            else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
            {
              data[_constrainedDofMap.at(dof)] =
                  promote<3>(value).dot(facet.normal()) * facet.volume();
            }
          }
        }
      }
    }
    return *this;
  }

  BCEss<FESpace> & operator<<(scalarFun_T const & f)
  {
    static_assert(
        FESpace_T::dim == 1, "scalar functions can be used only on scalar fe spaces.");
    return operator<<([f](Vec3 const & p) { return Vec1::Constant(f(p)); });
  }

  BCEss<FESpace> & operator<<(Vec const & v)
  {
    for (auto const & [dof, id]: _constrainedDofMap)
    {
      data[id] = v[dof];
    }
    return *this;
  }

  bool isConstrained(DOFid_T const id, int const d = 0) const
  {
    return _constrainedDofMap.count(id + d * feSpace.dof.size) > 0;
  }

  double get(DOFid_T const id) const { return data[_constrainedDofMap.at(id)]; }

  void makeTangent()
  {
    static_assert(dim > 1, "cannot make tangent a scalar field.");
    for (auto const facetId: facetIdList)
    {
      auto const & facet = feSpace.mesh.facetList[facetId];
      auto const normal = narrow<dim>(facet.normal());
      auto const & [elem, side] = facet.facingElem[0];
      feSpace.curFE.reinit(*elem);
      for (auto const dofFacet: RefFE_T::dofOnFacet[side])
      {
        // get the vector value on the dof
        FVec<dim> localValue;
        std::array<DOFid_T, dim> dofIds;
        for (uint d = 0; d < dim; ++d)
        {
          dofIds[d] = feSpace.dof.getId(elem->id, dofFacet, d);
          localValue[d] = data[_constrainedDofMap[dofIds[d]]];
        }

        // take the tangential component
        localValue =
            (FMat<dim, dim>::Identity() - normal * normal.transpose()) * localValue;

        // save back
        for (uint d = 0; d < dim; ++d)
        {
          data[_constrainedDofMap[dofIds[d]]] = localValue[d];
        }
      }
    }
  }

  friend std::ostream & operator<<(std::ostream & out, BCEss<FESpace_T> const & bc)
  {
    out << "constrainedDOFset: ";
    for ([[maybe_unused]] auto const [dof, id]: bc._constrainedDofMap)
    {
      out << dof << " ";
    }
    out << std::endl;
    return out;
  }

public:
  FESpace_T const & feSpace;
  marker_T const marker = markerNotSet;
  std::vector<short_T> const comp = allComp<FESpace>();
  double diag = 1.0;

  // protected:
  DofMap_T _constrainedDofMap;
  std::vector<id_T> facetIdList;

public:
  Vec data;
};

template <typename FESpace>
struct BCNat
{
  using FESpace_T = FESpace;
  using FacetFE_T = typename FESpace_T::RefFE_T::FacetFE_T;
  using QR_T = SideQR_T<typename FESpace_T::QR_T>;
  using CurFE_T = CurFE<FacetFE_T, QR_T>;

  BCNat(
      marker_T const m,
      Fun<FESpace::dim, 3> const & f,
      std::vector<short_T> const & c = allComp<FESpace>()):
      marker(m),
      value(f),
      comp(c)
  {
    // TODO: create a list of constrained faces at the beginning?
  }

  BCNat(marker_T const m, scalarFun_T const f, std::vector<short_T> const c):
      BCNat{m, [f](Vec3 const & p) { return Vec1(f(p)); }, c}
  {
    static_assert(
        FESpace::dim == 1, "this BC constructor cannot be used on vectorial FESpaces.");
  }

  bool hasComp(uint c) { return std::find(comp.begin(), comp.end(), c) != comp.end(); }

  CurFE_T curFE;
  marker_T marker;
  Fun<FESpace::dim, 3> const value;
  // TODO: convert to std::unordered_set
  std::vector<short_T> const comp = allComp<FESpace>();
};

// this class deals only with the matrix part of the mixed bc, the rhs part is
// managed as a standard natural bc
template <typename FESpace>
struct BCMixed
{
  using FESpace_T = FESpace;
  using FacetFE_T = typename FESpace_T::RefFE_T::FacetFE_T;
  using QR_T = SideQR_T<typename FESpace_T::QR_T>;
  using CurFE_T = CurFE<FacetFE_T, QR_T>;
  using RefFE_T = typename CurFE_T::RefFE_T;

  explicit BCMixed(
      marker_T m,
      Fun<FESpace_T::dim, 3> const f,
      std::vector<short_T> const c = allComp<FESpace>()):
      marker(m),
      coef(f),
      comp(c)
  {
    // TODO: create a list of constrained faces at the beginning?
  }

  bool hasComp(uint c) { return std::find(comp.begin(), comp.end(), c) != comp.end(); }

  CurFE_T curFE;
  marker_T marker;
  Fun<FESpace_T::dim, 3> const coef;
  // TODO: convert to std::unordered_set
  std::vector<short_T> const comp = allComp<FESpace>();
};

template <typename FESpace>
class BCList
{
public:
  using FESpace_T = FESpace;

  explicit BCList(FESpace const & fe): feSpace(fe) {}

  template <typename BC>
  void addBC(BC const bc)
  {
    static_assert(
        std::is_same_v<typename BC::FESpace_T, FESpace>,
        "this BC does not use the same FESpace");
    if (bc.marker != markerNotSet)
    {
      if (fixedMarkers.contains(bc.marker))
      {
        std::cerr << "the marker " << bc.marker << " has already been fixed."
                  << std::endl;
        abort();
      }
      fixedMarkers.insert(bc.marker);
    }

    if constexpr (std::is_same_v<typename BC::FESpace_T, FESpace>)
    {
      bcEssList.push_back(std::move(bc));
    }
    else
    {
      // this should never happen
      std::cerr << "a BC of type " << typeid(bc).name() << " has been added."
                << std::endl;
      std::abort();
    }
  }

  FESpace const & feSpace;
  std::unordered_set<marker_T> fixedMarkers;
  // std::unordered_set<id_T> fixedDofs;
  std::list<BCEss<FESpace>> bcEssList;
};

template <typename FESpace>
std::ostream & operator<<(std::ostream & out, BCList<FESpace> const & bcList)
{
  out << "bc list with " << bcList.bcEssList.size() << " essential bcs" << std::endl;
  for (auto & bc: bcList.bcEssList)
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
  using Predicate_T = std::function<bool(Vec3 const &)>;

  DOFCoordSet(FESpace_T & fe, Predicate_T const & p): feSpace(fe), predicate(p)
  {
    id_T counter = 0;
    for (auto const & e: feSpace.mesh.elementList)
    {
      feSpace.curFE.reinit(e);
      for (uint d = 0; d < numDOFs<typename FESpace_T::RefFE_T>(); ++d)
      {
        if (predicate(feSpace.curFE.dofPts[d]))
        {
          ids.insert({feSpace.dof.getId(e.id, d), counter});
          counter++;
        }
      }
    }
    std::cout << "new DOFCoordSet with " << ids.size() << " dofs" << std::endl;
  }

  FESpace_T & feSpace;
  Predicate_T const & predicate;
  typename BCEss<FESpace_T>::DofMap_T ids;
};

// apply a list of bcs to an already-assembled vector
template <typename BCList>
void applyBCs(Vec & v, BCList bcs)
{
  static_for(
      bcs,
      [&v]([[maybe_unused]] auto const i, auto const & bc)
      {
        for (auto const & [dof, counter]: bc._constrainedDofMap)
        {
          v[dof] = bc.data[counter];
        }
      });
}
