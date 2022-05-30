#pragma once

#include "def.hpp"

#include "bc.hpp"
#include "blockmatrix.hpp"
#include "fespace.hpp"
#include "var.hpp"

struct ScalarCoef
{
  explicit ScalarCoef(double const c): coef{c} {}

  void reinit(GeoElem const & /*elem*/) {}

  double evaluate(short_T const /*q*/) const { return coef; }

  double coef;
};

struct AssemblyBase
{
  // TODO: convert to std::unordered_set
  using CompList = std::vector<short_T>;

  CompList const comp = {};

  bool hasComp(short_T c) const
  {
    return std::find(comp.begin(), comp.end(), c) != comp.end();
  }
};

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
    for (short_T q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (short_T d = 0; d < FESpace_T::dim; ++d)
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
    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      double const coefQPoint = coef->evaluate(q)[0];
      for (uint d = 0; d < FESpace_T::dim; ++d)
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
    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (uint di = 0; di < FESpace_T::dim; ++di)
      {
        // d u_i / d x_j * d \phi_i / d x_j
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
            di * CurFE_T::numDOFs, di * CurFE_T::numDOFs) +=
            coef * this->feSpace->curFE.JxW[q] * this->feSpace->curFE.dphi[q] *
            this->feSpace->curFE.dphi[q].transpose();
        for (uint dj = 0; dj < FESpace_T::dim; ++dj)
        {
          // d u_j / d x_i * d \phi_i / d x_j
          Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>(
              di * CurFE_T::numDOFs, dj * CurFE_T::numDOFs) +=
              coef * this->feSpace->curFE.JxW[q] *
              this->feSpace->curFE.dphi[q].col(di) *
              this->feSpace->curFE.dphi[q].col(dj).transpose();
          // for (uint i=0; i<CurFE_T::numDOFs; ++i)
          //   for (uint j=0; j<CurFE_T::numDOFs; ++j)
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
struct AssemblyTensorStiffnessRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyTensorStiffnessRhs(
      double const c,
      Vec const & uOld,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      coef(c),
      velOld(&uOld)
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

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d = 0; d < FESpace_T::dim; ++d)
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
      coef(c),
      coefFE(&cFE)
  {}

  void reinit(GeoElem const & elem) const override { coefFE->reinit(elem); }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      double const coefQPoint = coefFE->evaluate(q);
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

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
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

template <typename FESpace>
struct AssemblyAnalyticRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyAnalyticRhs(
      Fun<FESpace::physicalDim(), 3> const r,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      rhs(std::move(r))
  {}

  AssemblyAnalyticRhs(
      scalarFun_T const r,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyAnalyticRhs<FESpace_T>(
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
        Vec3 localRhs3 = promote<3>(localRhs);
        Fe += this->feSpace->curFE.JxW[q] * this->feSpace->curFE.phiVect[q] * localRhs3;
      }
      else
      {
        std::abort();
      }
    }
  }

  Fun<FESpace::physicalDim(), 3> const rhs;
};

template <typename FESpace, typename FESpaceVel>
struct AssemblyAdvection: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceVel_T = FESpaceVel;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;

  AssemblyAdvection(
      FESpaceVel_T const & feVel,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>{fe, cmp},
      feSpaceVel{feVel}
  {}

  AssemblyAdvection(
      double const c,
      Vec const & u,
      FESpaceVel_T const & feVel,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      coef(c),
      vel(&u),
      feSpaceVel(&feVel)
  {}

  void reinit(GeoElem const & elem) const override { feSpaceVel->curFE.reinit(elem); }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      FVec<3> localVel = FVec<3>::Zero();
      if constexpr (family_v<typename FESpaceVel_T::RefFE_T> == FamilyType::LAGRANGE)
      {
        for (uint n = 0; n < FESpaceVel_T::RefFE_T::numDOFs; ++n)
        {
          for (uint d = 0; d < FESpaceVel_T::dim; ++d)
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
        for (uint n = 0; n < FESpaceVel_T::RefFE_T::numDOFs; ++n)
        {
          id_T const dofId = feSpaceVel->dof.getId(this->feSpace->curFE.elem->id, n);
          localVel +=
              (*vel)[dofId] * this->feSpace->curFE.phiVect[q].row(n).transpose();
        }
      }
      else
      {
        std::abort();
      }

      // TODO: the block is independent from d, so it could be computed once and then
      // copied over to the others
      for (uint d = 0; d < FESpace_T::dim; ++d)
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

  double coef;
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
      coef(c),
      vel(&u)
  {}

  void reinit(GeoElem const & elem) const override { vel->reinit(elem); }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    uint constexpr size = CurFE_T::numDOFs;
    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      auto const localV = vel->evaluate(q);
      FVec<3> const localVel = promote<3>(localV);
      for (uint d = 0; d < FESpace_T::dim; ++d)
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
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    // using CurFE2_T = typename FESpace2_T::CurFE_T;
    for (uint q = 0; q < CurFE1_T::QR_T::numPts; ++q)
    {
      Ke += coef * this->feSpace1->curFE.JxW[q] * this->feSpace1->curFE.divphi[q] *
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
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    // using CurFE2_T = typename FESpace2_T::CurFE_T;
    for (uint q = 0; q < CurFE1_T::QR_T::numPts; ++q)
    {
      Ke += coef * this->feSpace1->curFE.JxW[q] * this->feSpace1->curFE.phi[q] *
            this->feSpace2->curFE.divphi[q].transpose();
    }
  }

  double const coef;
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
      assembly = std::make_unique<AssemblyS2SProjection<FESpace_T, FESpaceRhs_T>>(
          AssemblyS2SProjection<FESpace_T, FESpaceRhs_T>{
              c, r, feRhs, fe, allComp<FESpace_T>()});
    }
    else if constexpr (
        fedim_v<typename FESpace_T::RefFE_T> == FEDimType::VECTOR &&
        fedim_v<typename FESpaceRhs_T::RefFE_T> == FEDimType::SCALAR)
    {
      assembly = std::make_unique<AssemblyS2VProjection<FESpace_T, FESpaceRhs_T>>(
          AssemblyS2VProjection<FESpace_T, FESpaceRhs_T>{c, r, feRhs, fe});
    }
    else if constexpr (
        fedim_v<typename FESpace_T::RefFE_T> == FEDimType::SCALAR &&
        fedim_v<typename FESpaceRhs_T::RefFE_T> == FEDimType::VECTOR)
    {
      // static_assert (dependent_false_v<FESpace>, "not yet implemented");
      assembly = std::make_unique<AssemblyV2SProjection<FESpace_T, FESpaceRhs_T>>(
          AssemblyV2SProjection<FESpace_T, FESpaceRhs_T>{c, r, feRhs, fe, cmp});
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

template <typename FESpace>
struct AssemblyBCNatural: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  using FacetFE_T = typename FESpace_T::RefFE_T::FacetFE_T;
  using QR_T = SideQR_T<typename FESpace_T::QR_T>;
  using FacetCurFE_T = CurFE<FacetFE_T, QR_T>;

  AssemblyBCNatural(
      Fun<FESpace_T::dim, 3> const r,
      marker_T const m,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      rhs(std::move(r)),
      marker(m)
  {}

  AssemblyBCNatural(
      scalarFun_T const r,
      marker_T const m,
      FESpace_T const & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyBCNatural<FESpace_T>(
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
            for (uint i = 0; i < FacetFE_T::numDOFs; ++i)
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

  Fun<FESpace_T::dim, 3> const rhs;
  marker_T const marker;
  FacetCurFE_T mutable facetCurFE;
};

template <typename FESpace>
struct AssemblyBCNormal: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  using FacetFE_T = typename FESpace::RefFE_T::FacetFE_T;
  using QR_T = SideQR_T<typename FESpace::QR_T>;
  using FacetCurFE_T = CurFE<FacetFE_T, QR_T>;

  AssemblyBCNormal(
      scalarFun_T const r,
      marker_T const m,
      FESpace & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      AssemblyVector<FESpace>(fe, cmp),
      rhs(std::move(r)),
      marker(m)
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
            for (uint i = 0; i < FacetFE_T::numDOFs; ++i)
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

  scalarFun_T const rhs;
  marker_T const marker;
  FacetCurFE_T mutable facetCurFE;
};

template <typename FESpace>
struct AssemblyBCMixed: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  using FacetFE_T = typename FESpace::RefFE_T::FacetFE_T;
  using QR_T = SideQR_T<typename FESpace::QR_T>;
  using FacetCurFE_T = CurFE<FacetFE_T, QR_T>;

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
            for (uint i = 0; i < FacetFE_T::numDOFs; ++i)
            {
              auto const idI =
                  CurFE_T::RefFE_T::dofOnFacet[facetCounter][i] + d * CurFE_T::numDOFs;
              for (uint j = 0; j < FacetFE_T::numDOFs; ++j)
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
