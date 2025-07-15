#pragma once

#include "def.hpp"

#include "assembly_base.hpp"
#include "reffe.hpp"
#include "var.hpp"

namespace proxpde
{

template <typename FESpace>
struct Diagonal: public AssemblyBase
{
  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace_T::CurFE_T;
  using LMat_T =
      FMat<FESpace_T::dim * CurFE_T::numDOFs, FESpace_T::dim * CurFE_T::numDOFs>;
  using LVec_T = FVec<FESpace_T::dim * CurFE_T::numDOFs>;

  Diagonal(FESpace_T const & fe, CompList const & cmp): AssemblyBase{cmp}, feSpace{&fe}
  {}

  virtual ~Diagonal() = default;

  virtual void build(LMat_T & Ke) const = 0;
  // virtual LMat_T build(/*LMat_T & Ke*/) const = 0;

  virtual void reinit(GeoElem const &) const {}

  FESpace_T const * feSpace;
  // LMat_T mat;
  // LVec_T vec;
};

template <typename FESpace>
struct AssemblyStiffness: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyStiffness(
      double const c,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace_T>(fe, cmp),
      coef(c)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (auto q = 0u; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (auto d = 0u; d < FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
            d * CurFE_T::numDOFs, d * CurFE_T::numDOFs) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.dphi[q] *
            this->feSpace->curFE.dphi[q].transpose();
      }
    }
  }

  // LMat_T build(/*LMat_T & Ke*/) const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double coef;
};

template <typename FESpace, typename Coef>
struct AssemblyStiffnessFE: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyStiffnessFE(
      Coef & c,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace_T>(fe, cmp),
      coef(&c)
  {}

  void reinit(GeoElem const & elem) const override { coef->reinit(elem); }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (auto q = 0u; q < CurFE_T::QR_T::numPts; ++q)
    {
      double const coefQPoint = coef->evaluate(q)[0];
      for (auto d = 0u; d < FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
            d * CurFE_T::numDOFs, d * CurFE_T::numDOFs) +=
            coefQPoint * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.dphi[q] *
            this->feSpace->curFE.dphi[q].transpose();
      }
    }
  }

private:
  Coef * coef;
};

template <typename FESpace>
struct AssemblyTensorStiffness: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyTensorStiffness(
      double const c,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace_T>(fe, cmp),
      coef(c)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (auto q = 0u; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (auto di = 0u; di < FESpace_T::dim; ++di)
      {
        // d u_i / d x_j * d \phi_i / d x_j
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
            di * CurFE_T::numDOFs, di * CurFE_T::numDOFs) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.dphi[q] *
            this->feSpace->curFE.dphi[q].transpose();
        for (auto dj = 0u; dj < FESpace_T::dim; ++dj)
        {
          // d u_j / d x_i * d \phi_i / d x_j
          Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
              di * CurFE_T::numDOFs, dj * CurFE_T::numDOFs) +=
              coef * this->feSpace->curFE.JxW[q] *
              this->feSpace->curFE.dphi[q].col(di) *
              this->feSpace->curFE.dphi[q].col(dj).transpose();
          // for (auto i = 0u; i < CurFE_T::numDOFs; ++i)
          //   for (auto j = 0u; j < CurFE_T::numDOFs; ++j)
          //   {
          //     // d u_j / d x_i * d \phi_i / d x_j
          //     Ke(i+di*CurFE_T::numDOFs, j+dj*CurFE_T::numDOFs) +=
          //         coef * this->feSpace->curFE.JxW[q] *
          //         this->feSpace->curFE.dphi[q].row(j)(di) *
          //         this->feSpace->curFE.dphi[q].row(i)(dj);
          //   }
        }
      }
    }
  }

  // LMat_T build(/*LMat_T & Ke*/) const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double coef;
};

template <typename FESpace>
struct AssemblyDummy: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyDummy(
      FESpace_T const & fe, AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp)
  {}

  void build(LMat_T &) const override {}
};

template <typename FESpace>
struct AssemblyScalarMass: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyScalarMass(
      double const & c,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace_T>(fe, cmp),
      coef(c)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for (auto q = 0u; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (auto d = 0u; d < FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
            d * CurFE_T::numDOFs, d * CurFE_T::numDOFs) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] *
            this->feSpace->curFE.phi[q].transpose();
      }
    }
  }

  // LMat_T build(/*LMat_T & Ke*/) const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double coef;
};

template <typename FESpace, typename Coef>
struct AssemblyMassFE: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyMassFE(
      double const c,
      Coef & cFE,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      coef{c},
      coefFE{&cFE}
  {}

  AssemblyMassFE(
      Coef & cFE,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      coef{1.0},
      coefFE{&cFE}
  {}

  void reinit(GeoElem const & elem) const override { coefFE->reinit(elem); }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for (auto q = 0u; q < CurFE_T::QR_T::numPts; ++q)
    {
      double const coefQPoint = coefFE->evaluate(q)[0];
      for (auto const c: this->comp)
      {
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
            c * CurFE_T::numDOFs, c * CurFE_T::numDOFs) +=
            coef * coefQPoint * this->feSpace->curFE.JxW[q] *
            this->feSpace->curFE.phi[q] * this->feSpace->curFE.phi[q].transpose();
      }
    }
  }

  double const coef;

private:
  Coef * coefFE;
};

template <typename FESpace>
struct AssemblyVectorMass: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyVectorMass(
      double const & c,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      coef(c)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for (auto q = 0u; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d = 0; d < FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
            d * CurFE_T::numDOFs, d * CurFE_T::numDOFs) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phiVect[q] *
            this->feSpace->curFE.phiVect[q].transpose();
      }
    }
  }

  // LMat_T build() const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double coef;
};

template <typename FESpace, FEDimType feDim>
struct AssemblyMassSelector
{};

template <typename FESpace>
struct AssemblyMassSelector<FESpace, FEDimType::SCALAR>
{
  using type = AssemblyScalarMass<FESpace>;
};

template <typename FESpace>
struct AssemblyMassSelector<FESpace, FEDimType::VECTOR>
{
  using type = AssemblyVectorMass<FESpace>;
};

template <typename FESpace>
using AssemblyMassSelector_T =
    typename AssemblyMassSelector<FESpace, fedim_v<typename FESpace::RefFE_T>>::type;

template <typename FESpace>
struct AssemblyMass: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;

  AssemblyMass(
      double const & c,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      assembly{c, fe, cmp}
  {}

  void build(LMat_T & Ke) const { assembly.build(Ke); }

  AssemblyMassSelector_T<FESpace> assembly;
};

// TODO: implement 2nd order formulation with
// \vec{v} \cdot \nabla \vec{v} = \vec{v}^{n+1} \cdot \nabla \vec{v}^{n} +
// \vec{v}^{n} \cdot \nabla \vec{v}^{n+1} - \vec{v}^{n} \cdot \nabla \vec{v}^{n}
// that inserts values also at rhs
template <typename FESpace, typename FESpaceVel>
struct AssemblyAdvection: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceVel_T = FESpaceVel;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;

  AssemblyAdvection(
      double const c,
      Vec const & u,
      FESpaceVel_T const & feVel,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>{fe, cmp},
      coef{c},
      vel{&u},
      feSpaceVel{&feVel}
  {}

  AssemblyAdvection(
      double const c,
      FEVar<FESpaceVel_T> const & velFE,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      AssemblyAdvection{c, velFE.data, *velFE.feSpace, fe, cmp}
  {}

  AssemblyAdvection(
      Vec const & u,
      FESpaceVel_T const & feVel,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      AssemblyAdvection{1.0, u, feVel, fe, cmp}
  {}

  AssemblyAdvection(
      FEVar<FESpaceVel_T> const & velFE,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      AssemblyAdvection{1.0, velFE.data, *velFE.feSpace, fe, cmp}
  {}

  AssemblyAdvection(
      FESpaceVel_T const & feVel,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      AssemblyAdvection{1.0, nullptr, feVel, fe, cmp}
  {}

  void reinit(GeoElem const & elem) const override { feSpaceVel->curFE.reinit(elem); }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (auto q = 0u; q < CurFE_T::QR_T::numPts; ++q)
    {
      FVec<3> localVel = FVec<3>::Zero();
      if constexpr (family_v<typename FESpaceVel_T::RefFE_T> == FamilyType::LAGRANGE)
      {
        for (auto n = 0u; n < FESpaceVel_T::RefFE_T::numDOFs; ++n)
        {
          for (auto d = 0u; d < FESpaceVel_T::dim; ++d)
          {
            id_T const dofId =
                feSpaceVel->dof.getId(this->feSpace->curFE.elem->id, n, d);
            localVel[d] += (*vel)[dofId] * this->feSpace->curFE.phi[q](n);
          }
        }
      }
      else if constexpr (
          family_v<typename FESpaceVel_T::RefFE_T> == FamilyType::RAVIART_THOMAS)
      {
        for (auto n = 0u; n < FESpaceVel_T::RefFE_T::numDOFs; ++n)
        {
          id_T const dofId = feSpaceVel->dof.getId(this->feSpace->curFE.elem->id, n);
          localVel +=
              (*vel)[dofId] * this->feSpaceVel->curFE.phiVect[q].row(n).transpose();
        }
      }
      else
      {
        std::cerr << "only Lagrange and Raviart-Thomas are cyrrently implemented."
                  << std::endl;
        std::exit(PROXPDE_NOT_IMPLEMENTED);
      }

      // TODO: the block is independent from d, so it could be computed once and then
      // copied over to the others
      for (auto d = 0u; d < FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
            d * CurFE_T::numDOFs, d * CurFE_T::numDOFs) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] *
            (this->feSpace->curFE.dphi[q] * localVel).transpose();
      }
    }
  }

  // LMat_T build(/*LMat_T & Ke*/) const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double coef = 1.0;
  Vec const * vel;
  FESpaceVel_T const * feSpaceVel;
};

template <typename FESpace, typename Vel>
struct AssemblyAdvectionFE: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyAdvectionFE(
      double const c,
      Vel & u,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      coef{c},
      vel{&u}
  {}

  AssemblyAdvectionFE(
      Vel & u,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      coef{1.0},
      vel{&u}
  {}

  void reinit(GeoElem const & elem) const override { vel->reinit(elem); }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    uint constexpr size = CurFE_T::numDOFs;
    for (auto q = 0u; q < CurFE_T::QR_T::numPts; ++q)
    {
      auto const localV = vel->evaluate(q);
      FVec<3> const localVel = promote<3>(localV);
      for (auto d = 0u; d < FESpace_T::dim; ++d)
      {
        Ke.template block<size, size>(d * size, d * size) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phi[q] *
            (this->feSpace->curFE.dphi[q] * localVel).transpose();
      }
    }
  }

  double const coef;
  Vel * vel;
};

} // namespace proxpde
