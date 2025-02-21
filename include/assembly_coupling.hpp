#pragma once

#include "def.hpp"

#include "assembly_base.hpp"

namespace proxpde
{

template <typename FESpace1, typename FESpace2>
struct Coupling: public AssemblyBase
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using CurFE1_T = typename FESpace1_T::CurFE_T;
  using CurFE2_T = typename FESpace2_T::CurFE_T;
  using LMat_T =
      FMat<FESpace1_T::dim * CurFE1_T::numDOFs, FESpace2_T::dim * CurFE2_T::numDOFs>;
  using LVec_T = FVec<FESpace1_T::dim * CurFE1_T::numDOFs>;

  Coupling(FESpace1_T const & fe1, FESpace2_T const & fe2, CompList const & cmp):
      AssemblyBase{cmp},
      feSpace1{&fe1},
      feSpace2{&fe2}
  {}

  virtual ~Coupling() = default;

  virtual void build(LMat_T & Ke) const = 0;

  virtual void reinit(GeoElem const & elem) const final
  {
    feSpace1->curFE.reinit(elem);
    feSpace2->curFE.reinit(elem);
  }

  FESpace1_T const * feSpace1;
  FESpace2_T const * feSpace2;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyGrad: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyGrad(
      double const c,
      FESpace1_T & fe1,
      FESpace2_T & fe2,
      AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
      Coupling<FESpace1_T, FESpace2_T>(fe1, fe2, cmp),
      coef(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace1_T::CurFE_T::QR_T,
            typename FESpace2_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    uint d = 0;
    for (auto const c: this->comp)
    {
      auto Kec = Ke.template block<CurFE1_T::numDOFs, CurFE2_T::numDOFs>(
          d * CurFE1_T::numDOFs, 0);
      for (uint q = 0; q < FESpace1_T::CurFE_T::QR_T::numPts; ++q)
      {
        Kec += coef * this->feSpace1->curFE.JxW[q] *
               this->feSpace1->curFE.dphi[q].col(c) *
               this->feSpace2->curFE.phi[q].transpose();
      }
      d++;
    }
  }

  double const coef;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyDiv: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyDiv(
      double const c,
      FESpace1_T & fe1,
      FESpace2_T & fe2,
      AssemblyBase::CompList const & cmp = allComp<FESpace2_T>()):
      Coupling<FESpace1_T, FESpace2_T>(fe1, fe2, cmp),
      coef(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace1_T::CurFE_T::QR_T,
            typename FESpace2_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    uint d = 0;
    for (auto const c: this->comp)
    {
      auto Kec = Ke.template block<CurFE1_T::numDOFs, CurFE2_T::numDOFs>(
          0, d * CurFE2_T::numDOFs);
      for (uint q = 0; q < CurFE1_T::QR_T::numPts; ++q)
      {
        Kec += coef * this->feSpace1->curFE.JxW[q] * this->feSpace1->curFE.phi[q] *
               this->feSpace2->curFE.dphi[q].col(c).transpose();
      }
      d++;
    }
  }

  double const coef;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyVectorGrad: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyVectorGrad(
      double const c,
      FESpace1_T & fe1,
      FESpace2_T & fe2,
      AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
      Coupling<FESpace1_T, FESpace2_T>(fe1, fe2, cmp),
      coef(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace1_T::CurFE_T::QR_T,
            typename FESpace2_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    for (uint q = 0; q < FESpace1_T::QR_T::numPts; ++q)
    {
      Ke += coef * this->feSpace1->curFE.JxW[q] * this->feSpace1->curFE.divphi[q] *
            this->feSpace2->curFE.phi[q].transpose();
    }
  }

  double const coef;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyVectorGRAD: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyVectorGRAD(
      double const c,
      FESpace1_T & fe1,
      FESpace2_T & fe2,
      AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
      Coupling<FESpace1_T, FESpace2_T>(fe1, fe2, cmp),
      coef(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace1_T::CurFE_T::QR_T,
            typename FESpace2_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    for (uint q = 0; q < FESpace1_T::QR_T::numPts; ++q)
    {
      Ke += coef * this->feSpace1->curFE.JxW[q] * this->feSpace1->curFE.divphi0[q] *
            this->feSpace2->curFE.phi[q].transpose();
    }
  }

  double const coef;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyVectorDiv: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyVectorDiv(
      double const c,
      FESpace1_T & fe1,
      FESpace2_T & fe2,
      AssemblyBase::CompList const & cmp = allComp<FESpace2_T>()):
      Coupling<FESpace1_T, FESpace2_T>(fe1, fe2, cmp),
      coef(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace1_T::CurFE_T::QR_T,
            typename FESpace2_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    for (uint q = 0; q < FESpace1_T::QR_T::numPts; ++q)
    {
      Ke += coef * this->feSpace1->curFE.JxW[q] * this->feSpace1->curFE.phi[q] *
            this->feSpace2->curFE.divphi[q].transpose();
    }
  }

  double const coef;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyVectorDIV: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyVectorDIV(
      double const c,
      FESpace1_T & fe1,
      FESpace2_T & fe2,
      AssemblyBase::CompList const & cmp = allComp<FESpace2_T>()):
      Coupling<FESpace1_T, FESpace2_T>(fe1, fe2, cmp),
      coef(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
        std::is_same_v<
            typename FESpace1_T::CurFE_T::QR_T,
            typename FESpace2_T::CurFE_T::QR_T>,
        "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    for (uint q = 0; q < FESpace1_T::QR_T::numPts; ++q)
    {
      Ke += coef * this->feSpace1->curFE.JxW[q] * this->feSpace1->curFE.phi[q] *
            this->feSpace2->curFE.divphi0[q].transpose();
    }
  }

  double const coef;
};

} // namespace proxpde
