#pragma once

#include "def.hpp"

#include "assembly_lhs.hpp"
#include "assembly_rhs.hpp"

namespace proxpde
{

template <typename FESpace>
struct AssemblyBCNaturalAnalytic: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  using FEFacet_T = typename FESpace_T::RefFE_T::FEFacet_T;
  using QR_T = SideQR_T<typename FESpace_T::QR_T>;
  using FacetCurFE_T = typename CurFETraits<FEFacet_T, QR_T>::type;

  AssemblyBCNaturalAnalytic(
      Fun<FESpace_T::dim, 3> const r,
      marker_T const m,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      rhs(std::move(r)),
      marker(m)
  {}

  AssemblyBCNaturalAnalytic(
      scalarFun_T const r,
      marker_T const m,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyBCNaturalAnalytic<FESpace_T>(
          [r](Vec3 const & p) { return Vec1::Constant(r(p)); }, m, fe, cmp)
  {
    static_assert(FESpace_T::dim == 1);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    auto const & mesh = *this->feSpace->mesh;
    auto const & e = *this->feSpace->curFE.elem;
    uint facetCounter = 0;
    for (auto const facetId: mesh.elemToFacet[e.id])
    {
      if (facetId != dofIdNotSet && mesh.facetList[facetId].marker == marker)
      {
        auto const & facet = mesh.facetList[facetId];
        facetCurFE.reinit(facet);
        for (uint q = 0; q < QR_T::numPts; ++q)
        {
          auto const value = rhs(facetCurFE.qpoint[q]);
          for (uint const d: this->comp)
          {
            for (uint i = 0; i < FEFacet_T::numDOFs; ++i)
            {
              auto const id =
                  CurFE_T::RefFE_T::dofOnFacet[facetCounter][i] + d * CurFE_T::numDOFs;
              Fe[id] += facetCurFE.JxW[q] * facetCurFE.phi[q](i) * value[d];
            }
          }
        }
      }
      facetCounter++;
    }
  }

  Fun<FESpace_T::dim, 3> rhs;
  marker_T const marker;
  FacetCurFE_T mutable facetCurFE;
};

template <typename FESpace, typename FESpaceData>
struct AssemblyBCNaturalFE: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceData_T = FESpaceData;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  using FEFacet_T = typename FESpace_T::RefFE_T::FEFacet_T;
  using QR_T = SideQR_T<typename FESpace_T::QR_T>;
  using FacetCurFE_T = typename CurFETraits<FEFacet_T, QR_T>::type;

  static constexpr auto numPtsQR = SideQR_T<typename FESpace_T::QR_T>::numPts;

  AssemblyBCNaturalFE(
      FEVar<FESpaceData_T> & d,
      std::vector<uint> const components,
      marker_T const m,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>{fe, cmp},
      marker{m},
      data{&d},
      dataComponents{components}
  {
    assert(dataComponents.size() == cmp.size());
    // TODO: we should allow different fespaces, as long as they are on sides
    static_assert(std::is_same_v<
                  typename FESpaceData_T::QR_T,
                  SideGaussQR<typename FESpace_T::Mesh_T::Elem_T, numPtsQR>>);
    // TODO: the data vector is volumetric, it is not sized on the boundary
    // auto count = 0u;
    // for (auto const & facet: data->feSpace->mesh->facetList)
    //   if (facet.marker == marker)
    //     count++;
    // // check that data is coherent with marked facets
    // assert(data->data.size() == count);
  }

  AssemblyBCNaturalFE(FEVar<FESpaceData_T> & d, marker_T const m, FESpace_T const & fe):
      AssemblyVector<FESpace_T>{fe, allComp<FESpace_T>()},
      marker{m},
      data{&d},
      dataComponents(allComp<FESpaceData_T>())
  {
    static_assert(FESpaceData_T::dim == FESpace_T::dim);
    // TODO: we should allow different fespaces, as long as they are on sides
    static_assert(std::is_same_v<
                  typename FESpaceData_T::QR_T,
                  SideGaussQR<typename FESpace_T::Mesh_T::Elem_T, numPtsQR>>);
    // TODO: the data vector is volumetric, it is not sized on the boundary
    // auto count = 0u;
    // for (auto const & facet: data->feSpace->mesh->facetList)
    //   if (facet.marker == marker)
    //     count++;
    // // check that data is coherent with marked facets
    // assert(data->data.size() == count);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    auto const & mesh = *this->feSpace->mesh;
    auto const & e = *this->feSpace->curFE.elem;
    for (auto f = 0u; f < FESpace_T::Mesh_T::Elem_T::numFacets; f++)
    {
      auto const facetId = mesh.elemToFacet[e.id][f];
      if (facetId != dofIdNotSet && mesh.facetList[facetId].marker == marker)
      {
        data->reinit(e);
        auto const & facet = mesh.facetList[facetId];
        auto const [ePtr, side] = facet.facingElem[0];
        assert(ePtr->id == e.id);
        facetCurFE.reinit(facet);
        for (uint q = 0u; q < QR_T::numPts; ++q)
        {
          auto const dataQ = data->evaluate(numPtsQR * side + q);
          for (uint k = 0u; k < this->comp.size(); k++)
          {
            auto const comp = this->comp[k];
            auto const dataComp = dataComponents[k];
            for (uint i = 0u; i < FEFacet_T::numDOFs; ++i)
            {
              auto const id =
                  CurFE_T::RefFE_T::dofOnFacet[f][i] + comp * CurFE_T::numDOFs;
              Fe[id] += facetCurFE.JxW[q] * facetCurFE.phi[q](i) * dataQ[dataComp];
            }
          }
        }
      }
    }
  }

  marker_T const marker;
  FEVar<FESpaceData_T> * data;
  std::vector<uint> dataComponents;
  FacetCurFE_T mutable facetCurFE;
};

template <typename FESpace>
struct AssemblyBCNormal: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  using FEFacet_T = typename FESpace::RefFE_T::FEFacet_T;
  using QR_T = SideQR_T<typename FESpace::QR_T>;
  using FacetCurFE_T = typename CurFETraits<FEFacet_T, QR_T>::type;

  AssemblyBCNormal(
      scalarFun_T const & r,
      marker_T const m,
      FESpace & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      AssemblyVector<FESpace>{fe, cmp},
      rhs{r},
      marker{m}
  {}

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    auto const & mesh = *this->feSpace->mesh;
    auto const & e = *this->feSpace->curFE.elem;
    uint facetCounter = 0;
    for (auto const facetId: mesh.elemToFacet[e.id])
    {
      if (facetId != dofIdNotSet && mesh.facetList[facetId].marker == marker)
      {
        auto const & facet = mesh.facetList[facetId];
        auto const normal = facet.normal();
        facetCurFE.reinit(facet);
        for (uint q = 0; q < QR_T::numPts; ++q)
        {
          for (auto const d: this->comp)
          {
            auto const value = rhs(facetCurFE.qpoint[q]);
            for (uint i = 0; i < FEFacet_T::numDOFs; ++i)
            {
              auto const id =
                  CurFE_T::RefFE_T::dofOnFacet[facetCounter][i] + d * CurFE_T::numDOFs;
              Fe(id) += facetCurFE.JxW[q] * facetCurFE.phi[q](i) * normal[d] * value;
            }
          }
        }
      }
      facetCounter++;
    }
  }

  scalarFun_T rhs;
  marker_T marker = markerNotSet;
  FacetCurFE_T mutable facetCurFE;
};

template <typename FESpace>
struct AssemblyBCMixed: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  using FEFacet_T = typename FESpace::RefFE_T::FEFacet_T;
  using QR_T = SideQR_T<typename FESpace::QR_T>;
  using FacetCurFE_T = typename CurFETraits<FEFacet_T, QR_T>::type;

  AssemblyBCMixed(
      scalarFun_T const c,
      marker_T const m,
      FESpace & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      coef(std::move(c)),
      marker(m)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    auto const & mesh = *this->feSpace->mesh;
    uint facetCounter = 0;
    for (auto const facetId: mesh.elemToFacet[this->feSpace->curFE.elem->id])
    {
      if (facetId != dofIdNotSet && mesh.facetList[facetId].marker == marker)
      {
        auto const & facet = mesh.facetList[facetId];
        facetCurFE.reinit(facet);
        for (uint q = 0; q < QR_T::numPts; ++q)
        {
          auto const localCoef = coef(facetCurFE.qpoint[q]);
          for (uint const d: this->comp)
          {
            for (uint i = 0; i < FEFacet_T::numDOFs; ++i)
            {
              auto const idI =
                  CurFE_T::RefFE_T::dofOnFacet[facetCounter][i] + d * CurFE_T::numDOFs;
              for (uint j = 0; j < FEFacet_T::numDOFs; ++j)
              {
                auto const idJ = CurFE_T::RefFE_T::dofOnFacet[facetCounter][j] +
                                 d * CurFE_T::numDOFs;
                Ke(idI, idJ) += facetCurFE.JxW[q] * localCoef * facetCurFE.phi[q](i) *
                                facetCurFE.phi[q](j);
              }
            }
          }
        }
      }
      facetCounter++;
    }
  }

  scalarFun_T const coef;
  marker_T const marker;
  FacetCurFE_T mutable facetCurFE;
};

} // namespace proxpde
