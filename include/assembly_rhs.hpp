#pragma once

#include "def.hpp"

#include "assembly_base.hpp"
#include "reffe.hpp"

namespace proxpde
{

template <typename FESpace>
struct AssemblyVector: public AssemblyBase
{
  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace::CurFE_T;
  using LMat_T =
      FMat<FESpace_T::dim * CurFE_T::numDOFs, FESpace_T::dim * CurFE_T::numDOFs>;
  using LVec_T = FVec<FESpace_T::dim * CurFE_T::numDOFs>;

  AssemblyVector(FESpace_T const & fe, CompList const & cmp):
      AssemblyBase{cmp},
      feSpace{&fe}
  {}

  virtual ~AssemblyVector() = default;

  virtual void build(LVec_T & Fe) const = 0;

  virtual void reinit(GeoElem const &) const {}

  FESpace_T const * feSpace;
};

template <typename FESpace>
struct AssemblyRhsAnalytic: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyRhsAnalytic(
      Fun<FESpace::physicalDim(), 3> const r,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      rhs(std::move(r))
  {}

  AssemblyRhsAnalytic(
      scalarFun_T const r,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyRhsAnalytic<FESpace_T>(
          [r](Vec3 const & p) { return Vec1(r(p)); }, fe, cmp)
  {
    static_assert(
        FESpace_T::dim == 1, "this ctor is available only for scalar fe spaces.");
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      FVec<FESpace::physicalDim()> const localRhs = rhs(this->feSpace->curFE.qpoint[q]);

      if constexpr (family_v<typename FESpace_T::RefFE_T> == FamilyType::LAGRANGE)
      {
        for (uint const c: this->comp)
        {
          Fe.template block<CurFE_T::numDOFs, 1>(c * CurFE_T::numDOFs, 0) +=
              this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] * localRhs[c];
        }
      }
      else if constexpr (
          family_v<typename FESpace_T::RefFE_T> == FamilyType::RAVIART_THOMAS)
      {
        auto const localRhs3 = promote<3u>(localRhs);
        Fe += this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phiVect[q] * localRhs3;
      }
      else
      {
        std::cerr << "only Lagrange and Raviart-Thomas are cyrrently implemented."
                  << std::endl;
        std::exit(PROXPDE_NOT_IMPLEMENTED);
      }
    }
  }

  Fun<FESpace::physicalDim(), 3> const rhs;
};

template <typename FESpace, typename FESpaceRhs = FESpace>
struct AssemblyS2SProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyS2SProjection(
      double const c,
      Vec const & r,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      rhs(&r),
      feSpaceRhs(&fe)
  {}

  AssemblyS2SProjection(
      double const c,
      Vec const & r,
      FESpaceRhs_T const & feRhs,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      rhs(&r),
      feSpaceRhs(&feRhs)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace_T::CurFE_T::QR_T,
            typename FESpaceRhs_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    // if constexpr (!std::is_same_v<FESpace_T,FESpaceRhs_T>)
    feSpaceRhs->curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFERhs_T = typename FESpaceRhs_T::CurFE_T;
    uint d = 0;
    for (auto const c: this->comp)
    {
      FVec<CurFERhs_T::numDOFs> rhsLocal;
      for (uint n = 0; n < CurFERhs_T::RefFE_T::numDOFs; ++n)
      {
        id_T const dofId = feSpaceRhs->dof.getId(feSpaceRhs->curFE.elem->id, n, c);
        rhsLocal[n] = (*rhs)[dofId];
      }
      auto Fec = Fe.template block<CurFE_T::numDOFs, 1>(d * CurFE_T::numDOFs, 0);
      for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
      {
        auto const rhsQ = feSpaceRhs->curFE.phi[q].dot(rhsLocal);
        Fec += coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] * rhsQ;
      }
      d++;
    }
  }

  double coef;
  Vec const * rhs;
  FESpaceRhs_T const * feSpaceRhs;
};

template <typename FESpace, typename FESpaceRhs>
struct AssemblyS2VProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyS2VProjection(
      double const c, Vec const & r, FESpaceRhs_T const & feRhs, FESpace_T const & fe):
      AssemblyVector<FESpace_T>(fe, {0}),
      coef(c),
      rhs(&r),
      feSpaceRhs(&feRhs)
  {
    static_assert(
        FESpace_T::RefFE_T::dim == FESpaceRhs_T::dim,
        "the two fespaces are not of the same dimension");
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace_T::CurFE_T::QR_T,
            typename FESpaceRhs_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override { feSpaceRhs->curFE.reinit(elem); }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFERhs_T = typename FESpaceRhs_T::CurFE_T;
    for (uint d = 0; d < CurFE_T::RefFE_T::dim; ++d)
    {
      FVec<CurFERhs_T::numDOFs> localRhs;
      for (uint n = 0; n < CurFERhs_T::numDOFs; ++n)
      {
        id_T const dofId = feSpaceRhs->dof.getId(feSpaceRhs->curFE.elem->id, n, d);
        localRhs[n] = (*rhs)[dofId];
      }
      for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
      {
        Fe += coef * this->feSpace->curFE.JxW[q] *
              this->feSpace->curFE.phiVect[q].col(d) *
              (feSpaceRhs->curFE.phi[q].dot(localRhs)); // u*[q] = sum_k u*_k phi_k[q]
      }
    }
  }

  double const coef;
  Vec const * rhs;
  FESpaceRhs_T const * feSpaceRhs;
};

template <typename FESpace, typename FESpaceRhs>
struct AssemblyV2SProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyV2SProjection(
      double const c,
      Vec const & r,
      FESpaceRhs_T const & feRhs,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      rhs(&r),
      feSpaceRhs(&feRhs)
  {
    // TODO: when comp is specified, this is not guaranteed anymore
    // static_assert(FESpace_T::dim == FESpaceRhs_T::RefFE_T::dim,
    //               "the two fespaces are not of the same dimension");
    assert(cmp.size() == FESpace_T::dim);

    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace_T::CurFE_T::QR_T,
            typename FESpaceRhs_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override { feSpaceRhs->curFE.reinit(elem); }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFERhs_T = typename FESpaceRhs_T::CurFE_T;

    FVec<CurFERhs_T::numDOFs> localRhs;
    for (uint n = 0; n < CurFERhs_T::RefFE_T::numDOFs; ++n)
    {
      id_T const dofId = feSpaceRhs->dof.getId(feSpaceRhs->curFE.elem->id, n);
      localRhs[n] = (*rhs)[dofId];
    }

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d = 0; d < this->comp.size(); ++d)
      {
        auto const c = this->comp[d];
        // u*[q] = sum_k u*_k vec{psi}_k[q]
        double rhsQ = feSpaceRhs->curFE.phiVect[q].col(c).dot(localRhs);

        Fe.template block<CurFE_T::numDOFs, 1>(d * CurFE_T::numDOFs, 0) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] * rhsQ;
      }
    }
  }

  double const coef;
  Vec const * rhs;
  FESpaceRhs_T const * feSpaceRhs;
};

template <typename FESpace, typename FESpaceRhs = FESpace>
struct AssemblyProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyProjection(
      double const c,
      Vec const & r,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp)
  {
    if constexpr (fedim_v<typename FESpace_T::RefFE_T> == FEDimType::SCALAR)
    {
      assembly = std::make_unique<AssemblyS2SProjection<FESpace_T>>(
          AssemblyS2SProjection<FESpace_T>{c, r, fe, cmp});
    }
    else
    {
      // no support for vector -> vector projection
      static_assert(
          dependent_false_v<FESpace>,
          "vector -> vector projection not yet implemented");
    }
  }

  AssemblyProjection(
      double const c,
      Vec const & r,
      FESpaceRhs_T const & feRhs,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace_T::CurFE_T::QR_T,
            typename FESpaceRhs_T::CurFE_T::QR_T>,
        "the two quad rules are not the same");

    if constexpr (
        fedim_v<typename FESpace_T::RefFE_T> == FEDimType::SCALAR &&
        fedim_v<typename FESpaceRhs_T::RefFE_T> == FEDimType::SCALAR)
    {
      assembly = std::make_shared<AssemblyS2SProjection<FESpace_T, FESpaceRhs_T>>(
          c, r, feRhs, fe, cmp);
    }
    else if constexpr (
        fedim_v<typename FESpace_T::RefFE_T> == FEDimType::VECTOR &&
        fedim_v<typename FESpaceRhs_T::RefFE_T> == FEDimType::SCALAR)
    {
      assembly = std::make_shared<AssemblyS2VProjection<FESpace_T, FESpaceRhs_T>>(
          c, r, feRhs, fe);
    }
    else if constexpr (
        fedim_v<typename FESpace_T::RefFE_T> == FEDimType::SCALAR &&
        fedim_v<typename FESpaceRhs_T::RefFE_T> == FEDimType::VECTOR)
    {
      // static_assert (dependent_false_v<FESpace>, "not yet implemented");
      assembly = std::make_shared<AssemblyV2SProjection<FESpace_T, FESpaceRhs_T>>(
          c, r, feRhs, fe, cmp);
    }
    else if constexpr (
        fedim_v<typename FESpace_T::RefFE_T> == FEDimType::VECTOR &&
        fedim_v<typename FESpaceRhs_T::RefFE_T> == FEDimType::VECTOR)
    {
      // no support for vector -> vector projection
      static_assert(
          dependent_false_v<FESpace>,
          "vector -> vector projection not yet implemented");
    }
  }

  void reinit(GeoElem const & elem) const override { assembly->reinit(elem); }

  void build(LVec_T & Fe) const override { assembly->build(Fe); }

  // LVec_T build() const
  // {
  //   return LVec_T{};
  // }

  // TODO: check own to avoid copies of the class
  // that would be implicitly deleted by the use of a unique_ptr
  std::shared_ptr<AssemblyVector<FESpace_T>> assembly;
};

template <typename FESpace, typename FESpaceData>
struct AssemblyDivRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceData_T = FESpaceData;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyDivRhs(
      double const c,
      Vec const & vec,
      FESpaceData_T const & feData,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = {}):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      data(&vec),
      feSpaceData(&feData)
  {
    static_assert(FESpace_T::dim == 1, "this assembly wors only on scalar fe spaces");
  }

  void reinit(GeoElem const & elem) const override { feSpaceData->curFE.reinit(elem); }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFEData_T = typename FESpaceData_T::CurFE_T;
    uint constexpr sizeData = CurFEData_T::numDOFs;

    FMat<FESpaceData_T::dim, sizeData> localData;
    for (uint d2 = 0; d2 < FESpaceData_T::dim; ++d2)
    {
      for (uint n = 0; n < sizeData; ++n)
      {
        id_T const dofId = feSpaceData->dof.getId(feSpaceData->curFE.elem->id, n, d2);
        localData(d2, n) = (*data)[dofId];
      }
    }

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d2 = 0; d2 < FESpaceData_T::dim; ++d2)
      {
        double const localDataGradQ =
            localData.row(d2) * feSpaceData->curFE.dphi[q].col(d2);
        Fe += coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] *
              localDataGradQ;
      }
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coef;
  Vec const * data;
  FESpaceData const * feSpaceData;
};

template <typename FESpace, typename FESpaceData>
struct AssemblyGradRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceData_T = FESpaceData;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyGradRhs(
      double const c,
      Vec const & vec,
      FESpaceData_T const & feData,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      data(&vec),
      feSpaceData(&feData)
  {}

  void reinit(GeoElem const & elem) const override { feSpaceData->curFE.reinit(elem); }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFEData_T = typename FESpaceData_T::CurFE_T;
    uint constexpr size = CurFE_T::numDOFs;
    uint constexpr sizeData = CurFEData_T::numDOFs;
    FMat<FESpaceData_T::dim, sizeData> localData;

    uint const elemId = feSpaceData->curFE.elem->id;
    for (uint d2 = 0; d2 < FESpaceData_T::dim; ++d2)
    {
      for (uint n = 0; n < sizeData; ++n)
      {
        id_T const dofId = feSpaceData->dof.getId(elemId, n, d2);
        localData(d2, n) = (*data)[dofId];
      }
    }

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      FMat<FESpaceData_T::dim, 3> const localDataGradQ =
          localData * feSpaceData->curFE.dphi[q];
      // FVec<FESpaceData_T::dim> const localDataQ =
      //     localData * feSpaceData.curFE.phi[q];

      uint counter = 0;
      for (auto const c: this->comp)
      {
        uint const d1 = c % FESpace_T::Mesh_T::Elem_T::dim;
        uint const d2 = c / FESpace_T::Mesh_T::Elem_T::dim;

        // d u_d2 / d x_d1
        Fe.template block<size, 1>(counter * size, 0) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] *
            localDataGradQ(d2, d1);
        // // u_d2
        // Fe.template block<size, 1>(counter * size, 0) +=
        //     coef * this->feSpace->curFE.JxW[q] *
        //     this->feSpace->curFE.dphi[q].col(d1) *
        //     localDataQ[d2];
        counter++;
      }
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coef;
  Vec const * data;
  FESpaceData_T const * feSpaceData;
};

template <typename FESpace, typename FESpaceData>
struct AssemblyGradRhs2: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceData_T = FESpaceData;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyGradRhs2(
      double const c,
      Vec const & vec,
      FESpaceData_T const & feData,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      feSpaceData(&feData),
      data(&vec)
  {}

  void reinit(GeoElem const & elem) const override { feSpaceData->curFE.reinit(elem); }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFEData_T = typename FESpaceData_T::CurFE_T;
    uint constexpr size = CurFE_T::numDOFs;
    uint constexpr sizeData = CurFEData_T::numDOFs;
    FMat<FESpaceData_T::dim, sizeData> localData;

    uint const elemId = feSpaceData->curFE.elem->id;
    for (uint d2 = 0; d2 < FESpaceData_T::dim; ++d2)
    {
      for (uint n = 0; n < sizeData; ++n)
      {
        id_T const dofId = feSpaceData->dof.getId(elemId, n, d2);
        localData(d2, n) = (*data)[dofId];
      }
    }

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      // FMat<FESpaceData_T::dim, 3> const localdataGradQ =
      //     localData * feSpaceData.curFE.dphi[q];
      FVec<FESpaceData_T::dim> const dataQ = localData * feSpaceData->curFE.phi[q];

      uint counter = 0;
      for (auto const c: this->comp)
      {
        uint const d1 = c % FESpace_T::Mesh_T::Elem_T::dim;
        uint const d2 = c / FESpace_T::Mesh_T::Elem_T::dim;

        // // d u_d2 / d x_d1
        // Fe.template block<size, 1>(counter * size, 0) +=
        //     coef * this->feSpace->curFE.JxW[q] *
        //     this->feSpace->curFE.phi[q] *
        //     localdataGradQ(d2, d1);
        // u_d2
        Fe.template block<size, 1>(counter * size, 0) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.dphi[q].col(d1) *
            dataQ[d2];
        counter++;
      }
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coef;
  FESpaceData_T const * feSpaceData;
  Vec const * data;
};

template <typename FESpace, typename FESpaceData = FESpace>
struct AssemblyStiffnessRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceData_T = FESpaceData;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyStiffnessRhs(
      double const c,
      Vec const & vec,
      FESpaceData_T & feData,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      data(&vec),
      feSpaceData(&feData)
  {}

  void reinit(GeoElem const & elem) const override { feSpaceData->curFE.reinit(elem); }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFEData_T = typename FESpaceData_T::CurFE_T;
    uint constexpr size = CurFE_T::numDOFs;
    uint constexpr sizeData = CurFEData_T::numDOFs;
    FMat<FESpaceData_T::dim, sizeData> localData;

    uint const elemId = feSpaceData->curFE.elem->id;
    for (uint d2 = 0; d2 < FESpaceData_T::dim; ++d2)
    {
      for (uint n = 0; n < sizeData; ++n)
      {
        id_T const dofId = feSpaceData->dof.getId(elemId, n, d2);
        localData(d2, n) = (*data)[dofId];
      }
    }

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      FMat<FESpaceData_T::dim, 3> const localDataGradQ =
          localData * feSpaceData->curFE.dphi[q];

      uint counter = 0;
      for (auto const c: this->comp)
      {
        Fe.template block<size, 1>(counter * size, 0) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.dphi[q].col(c) *
            localDataGradQ(c, c);
        counter++;
      }
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coef;
  Vec const * data;
  FESpaceData_T * feSpaceData;
};

template <typename FESpace>
struct AssemblyTensorStiffnessRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyTensorStiffnessRhs(
      double const c,
      Vec const & uOld,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Super_T{fe, cmp},
      coef{c},
      velOld{&uOld}
  {}

  void build(LVec_T & /*Fe*/) const override
  {
    // using CurFE_T = typename FESpace_T::CurFE_T;
    // for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    // {
    //   for (uint di=0; di<FESpace_T::dim; ++di)
    //   {
    //     // d u_i / d x_j * d \phi_i / d x_j
    //     Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(di*CurFE_T::numDOFs,
    //     di*CurFE_T::numDOFs) +=
    //         coef * this->feSpace->curFE.JxW[q] *
    //         this->feSpace->curFE.dphi[q] *
    //         this->feSpace->curFE.dphi[q].transpose();
    //     for (uint dj=0; dj<FESpace_T::dim; ++dj)
    //     {
    //       // d u_j / d x_i * d \phi_i / d x_j
    //       Ke.template
    //       block<CurFE_T::numDOFs,CurFE_T::numDOFs>(di*CurFE_T::numDOFs,
    //       dj*CurFE_T::numDOFs) +=
    //           coef * this->feSpace->curFE.JxW[q] *
    //           this->feSpace->curFE.dphi[q].col(di) *
    //           this->feSpace->curFE.dphi[q].col(dj).transpose();
    //       // for (uint i=0; i<CurFE_T::numDOFs; ++i)
    //       //   for (uint j=0; j<CurFE_T::numDOFs; ++j)
    //       //   {
    //       //     // d u_j / d x_i * d \phi_i / d x_j
    //       //     Ke(i+di*CurFE_T::numDOFs, j+dj*CurFE_T::numDOFs) +=
    //       //         coef * this->feSpace->curFE.JxW[q] *
    //       //         this->feSpace->curFE.dphi[q].row(j)(di) *
    //       //         this->feSpace->curFE.dphi[q].row(i)(dj);
    //       //   }
    //          }
    //        }
    //      }
  }

  double coef;
  Vec const * velOld;
};

template <typename FESpace, typename FESpaceVel, typename FESpaceData = FESpace>
struct AssemblyAdvectionRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceVel_T = FESpaceVel;
  using FESpaceData_T = FESpaceData;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyAdvectionRhs(
      double const c,
      Vec const & u,
      FESpaceVel_T const & feVel,
      Vec const & d,
      FESpaceData_T const & feData,
      FESpace_T & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      vel(&u),
      feSpaceVel(&feVel),
      data(&d),
      feSpaceData(&feData)
  {}

  AssemblyAdvectionRhs(
      double const c,
      Vec const & u,
      FESpaceVel_T const & feVel,
      Vec const & vec,
      FESpace_T & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      vel(&u),
      feSpaceVel(&feVel),
      data(&vec),
      feSpaceData(&fe)
  {}

  void reinit(GeoElem const & elem) const override
  {
    feSpaceData->curFE.reinit(elem);
    feSpaceVel->curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFEData_T = typename FESpaceData_T::CurFE_T;
    using CurFEVel_T = typename FESpaceVel_T::CurFE_T;

    FMat<FESpaceData_T::dim, CurFEData_T::numDOFs> localData;
    for (uint n = 0; n < CurFEData_T::RefFE_T::numDOFs; ++n)
    {
      for (uint d2 = 0; d2 < FESpaceData_T::dim; ++d2)
      {
        id_T const dofId = feSpaceData->dof.getId(feSpaceData->curFE.elem->id, n, d2);
        localData(d2, n) = (*data)[dofId];
      }
    }

    FMat<3, CurFEVel_T::numDOFs> localVel = FVec<3>::Zero();
    for (uint n = 0; n < CurFEVel_T::numDOFs; ++n)
    {
      for (uint d3 = 0; d3 < FESpaceVel_T::dim; ++d3)
      {
        id_T const dofId = feSpaceVel->dof.getId(feSpaceVel->curFE.elem->id, n, d3);
        localVel(d3, n) = (*vel)[dofId];
      }
    }

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      FVec<3> localVelQ = promote<3>(localVel * feSpaceVel->curFE.phi[q]);

      FMat<FESpaceData_T::dim, 3> const localDataGradQ =
          localData * feSpaceData->curFE.dphi[q];

      uint counter = 0;
      for (auto const c: this->comp)
      {
        Fe.template block<CurFE_T::numDOFs, 1>(counter * CurFE_T::numDOFs, 0) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] *
            (localDataGradQ * localVelQ)[c];
        counter++;
      }
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coef;
  Vec const * vel;
  FESpaceVel_T const * feSpaceVel;
  Vec const * data;
  FESpaceData_T const * feSpaceData;
};

} // namespace proxpde
