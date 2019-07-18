#pragma once

#include "def.hpp"
#include "blockmatrix.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"

struct ScalarCoef
{
  explicit ScalarCoef(double const c):
    coef{c}
  {}

  void reinit(GeoElem const & /*elem*/) {}

  double evaluate(uint const /*q*/) const
  {
    return coef;
  }

  double coef;
};

struct AssemblyBase
{
  using CompList = std::vector<uint>;

  CompList const comp = {};

  bool hasComp(uint c) const
  {
    return std::find(comp.begin(), comp.end(), c) != comp.end();
  }
};

template <typename FESpace>
struct Diagonal: public AssemblyBase
{
  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace_T::CurFE_T;
  using LMat_T = FMat<FESpace_T::dim*CurFE_T::numDOFs, FESpace_T::dim*CurFE_T::numDOFs>;
  using LVec_T = FVec<FESpace_T::dim*CurFE_T::numDOFs>;

  explicit Diagonal(FESpace_T & fe, CompList const & cmp):
    AssemblyBase{cmp},
    feSpace(fe)
  {}

  virtual ~Diagonal() = default;

  virtual void build(LMat_T & Ke) const = 0;
  // virtual LMat_T build(/*LMat_T & Ke*/) const = 0;

  virtual void reinit(GeoElem const & elem) const
  {
    feSpace.curFE.reinit(elem);
  }

  FESpace_T & feSpace;
  LMat_T mat;
  LVec_T vec;
};

template <typename FESpace1, typename FESpace2>
struct Coupling: public AssemblyBase
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using CurFE1_T = typename FESpace1_T::CurFE_T;
  using CurFE2_T = typename FESpace2_T::CurFE_T;
  using LMat_T = FMat<FESpace1_T::dim*CurFE1_T::numDOFs, FESpace2_T::dim*CurFE2_T::numDOFs>;
  using LVec_T = FVec<FESpace1_T::dim*CurFE1_T::numDOFs>;

  explicit Coupling(FESpace1_T & fe1,
                    FESpace2 & fe2,
                    CompList const & cmp):
    AssemblyBase{cmp},
    feSpace1(fe1),
    feSpace2(fe2)
  {}

  virtual ~Coupling() = default;

  virtual void build(LMat_T & Ke) const = 0;

  virtual void reinit(GeoElem const & elem) const final
  {
    feSpace1.curFE.reinit(elem);
    feSpace2.curFE.reinit(elem);
  }

  FESpace1_T & feSpace1;
  FESpace2_T & feSpace2;
};

template <typename FESpace>
struct AssemblyVector: public AssemblyBase
{
  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace::CurFE_T;
  using LMat_T = FMat<FESpace_T::dim*CurFE_T::numDOFs, FESpace_T::dim*CurFE_T::numDOFs>;
  using LVec_T = FVec<FESpace_T::dim*CurFE_T::numDOFs>;

  explicit AssemblyVector(FESpace_T & fe, CompList const & cmp):
    AssemblyBase{cmp},
    feSpace(fe)
  {}

  virtual ~AssemblyVector() = default;

  virtual void build(LVec_T & Fe) const = 0;

  virtual void reinit(GeoElem const & elem) const
  {
    feSpace.curFE.reinit(elem);
  }

  FESpace_T & feSpace;
};

template <typename FESpace>
struct AssemblyStiffness: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyStiffness(double const c,
                    FESpace_T & fe,
                    AssemblyBase::CompList const & cmp = allComp<FESpace>()):
    Diagonal<FESpace_T>(fe, cmp),
    coeff(c)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(d*CurFE_T::numDOFs, d*CurFE_T::numDOFs) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.dphi[q] *
            this->feSpace.curFE.dphi[q].transpose();
      }
    }
  }

  // LMat_T build(/*LMat_T & Ke*/) const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double coeff;
};

template <typename FESpace, typename Coef>
struct AssemblyStiffnessFE: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyStiffnessFE(Coef & c,
                      FESpace_T & fe,
                      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
    Diagonal<FESpace_T>(fe, cmp),
    coef(c)
  {}

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    coef.reinit(elem);
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double const coefQPoint = coef.evaluate(q);
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>
            (d*CurFE_T::numDOFs, d*CurFE_T::numDOFs) +=
            coefQPoint * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.dphi[q] *
            this->feSpace.curFE.dphi[q].transpose();
      }
    }
  }

private:
  Coef & coef;
};

template <typename FESpace>
struct AssemblyTensorStiffness: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyTensorStiffness(double const c,
                    FESpace_T & fe,
                    AssemblyBase::CompList const & cmp = allComp<FESpace>()):
    Diagonal<FESpace_T>(fe, cmp),
    coeff(c)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint di=0; di<FESpace_T::dim; ++di)
      {
        // d u_i / d x_j * d \phi_i / d x_j
        Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(di*CurFE_T::numDOFs, di*CurFE_T::numDOFs) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.dphi[q] *
            this->feSpace.curFE.dphi[q].transpose();
        for (uint dj=0; dj<FESpace_T::dim; ++dj)
        {
          // d u_j / d x_i * d \phi_i / d x_j
          Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(di*CurFE_T::numDOFs, dj*CurFE_T::numDOFs) +=
              coeff * this->feSpace.curFE.JxW[q] *
              this->feSpace.curFE.dphi[q].col(di) *
              this->feSpace.curFE.dphi[q].col(dj).transpose();
//          for (uint i=0; i<CurFE_T::numDOFs; ++i)
//            for (uint j=0; j<CurFE_T::numDOFs; ++j)
//            {
//              // d u_j / d x_i * d \phi_i / d x_j
//              Ke(i+di*CurFE_T::numDOFs, j+dj*CurFE_T::numDOFs) +=
//                  coeff * this->feSpace.curFE.JxW[q] *
//                  this->feSpace.curFE.dphi[q].row(j)(di) *
//                  this->feSpace.curFE.dphi[q].row(i)(dj);
//            }
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

  double coeff;
};

template <typename FESpace>
struct AssemblyTensorStiffnessRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyTensorStiffnessRhs(double const c,
                             Vec const & uOld,
                             FESpace_T & fe,
                             AssemblyBase::CompList const & cmp = allComp<FESpace>()):
    AssemblyVector<FESpace_T>(fe, cmp),
    coeff(c),
    velOld(uOld)
  {}

  void build(LVec_T & Fe) const override
  {
//     using CurFE_T = typename FESpace_T::CurFE_T;
//     for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
//     {
//       for (uint di=0; di<FESpace_T::dim; ++di)
//       {
//         // d u_i / d x_j * d \phi_i / d x_j
//         Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(di*CurFE_T::numDOFs, di*CurFE_T::numDOFs) +=
//             coeff * this->feSpace.curFE.JxW[q] *
//             this->feSpace.curFE.dphi[q] *
//             this->feSpace.curFE.dphi[q].transpose();
//         for (uint dj=0; dj<FESpace_T::dim; ++dj)
//         {
//           // d u_j / d x_i * d \phi_i / d x_j
//           Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(di*CurFE_T::numDOFs, dj*CurFE_T::numDOFs) +=
//               coeff * this->feSpace.curFE.JxW[q] *
//               this->feSpace.curFE.dphi[q].col(di) *
//               this->feSpace.curFE.dphi[q].col(dj).transpose();
// //          for (uint i=0; i<CurFE_T::numDOFs; ++i)
// //            for (uint j=0; j<CurFE_T::numDOFs; ++j)
// //            {
// //              // d u_j / d x_i * d \phi_i / d x_j
// //              Ke(i+di*CurFE_T::numDOFs, j+dj*CurFE_T::numDOFs) +=
// //                  coeff * this->feSpace.curFE.JxW[q] *
// //                  this->feSpace.curFE.dphi[q].row(j)(di) *
// //                  this->feSpace.curFE.dphi[q].row(i)(dj);
// //            }
//         }
//       }
//     }
  }

  double coeff;
  Vec const & velOld;
};

template <typename FESpace>
struct AssemblyMass: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyMass(double const & c,
                        FESpace_T & fe,
                        AssemblyBase::CompList const & cmp = allComp<FESpace>()):
     Diagonal<FESpace>(fe, cmp),
     coeff(c)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(d*CurFE_T::numDOFs, d*CurFE_T::numDOFs) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            this->feSpace.curFE.phi[q].transpose();
      }
    }
  }

  // LMat_T build(/*LMat_T & Ke*/) const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double coeff;
};

template <typename FESpace, typename Coef>
struct AssemblyMassFE: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyMassFE(double const c,
                          Coef & cFE,
                          FESpace_T & fe,
                          AssemblyBase::CompList const & cmp = allComp<FESpace>()):
     Diagonal<FESpace>(fe, cmp),
     coef(c),
     coefFE(cFE)
  {}

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    coefFE.reinit(elem);
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double const coefQPoint = coefFE.evaluate(q);
      for (auto const c: this->comp)
      {
        Ke.template block<CurFE_T::numDOFs, CurFE_T::numDOFs>
            (c*CurFE_T::numDOFs, c*CurFE_T::numDOFs) +=
            coef * coefQPoint *
            this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            this->feSpace.curFE.phi[q].transpose();
      }
    }
  }

  double const coef;

private:
  Coef & coefFE;
};

template <typename FESpace>
struct AssemblyVectorMass: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyVectorMass(double const & c,
                              FESpace_T & fe,
                              AssemblyBase::CompList const & cmp = allComp<FESpace>()):
     Diagonal<FESpace>(fe, cmp),
     coeff(c)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(d*CurFE_T::numDOFs, d*CurFE_T::numDOFs) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phiVect[q] *
            this->feSpace.curFE.phiVect[q].transpose();
      }
    }
  }

  // LMat_T build() const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double coeff;
};

template <typename FESpace>
struct AssemblyAnalyticRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyAnalyticRhs(Fun<FESpace_T::dim,3> const r,
                      FESpace & fe,
                      AssemblyBase::CompList const cmp = allComp<FESpace>()):
    AssemblyVector<FESpace>(fe, std::move(cmp)),
    rhs(std::move(r))
  {}

  // this constructor is available only when FESpace::dim == 1
  template <typename FESpace1 = FESpace,
            std::enable_if_t<FESpace1::dim == 1, bool> = true>
  AssemblyAnalyticRhs(scalarFun_T const r,
                      FESpace & fe,
                      AssemblyBase::CompList const cmp = allComp<FESpace>()):
    AssemblyAnalyticRhs<FESpace>([r](Vec3 const & p) {return Vec1(r(p));}, fe, std::move(cmp))
  {}

  // otherwise, we fail
  template <typename FESpace1 = FESpace,
            std::enable_if_t<FESpace1::dim != 1, bool> = true>
  AssemblyAnalyticRhs(scalarFun_T const &,
                      FESpace & fe,
                      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
    AssemblyVector<FESpace>(fe, cmp)
  {
    std::abort();
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Fe.template block<CurFE_T::numDOFs, 1>(d*CurFE_T::numDOFs, 0) +=
            this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            rhs(this->feSpace.curFE.qpoint[q])(d);
      }
    }
  }

  Fun<FESpace_T::dim,3> const rhs;
};

template <typename FESpace, typename FESpaceVel>
struct AssemblyAdvection: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using VelVec_T = Eigen::Matrix<double, Eigen::Dynamic, FESpace_T::RefFE_T::dim>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyAdvection(double const c,
                             Vec const & u,
                             FESpaceVel & feVel,
                             FESpace_T & fe,
                             AssemblyBase::CompList const & cmp = allComp<FESpace>()):
    Diagonal<FESpace>(fe, cmp),
    coeff(c),
    vel(u),
    feSpaceVel(feVel)
  {}

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    // if (&feSpaceVel != &this->feSpace)
    {
      feSpaceVel.curFE.reinit(elem);
    }
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      //TODO: the velocity should have its own fe space
      FVec<3> localVel = FVec<3>::Zero();
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        for (uint d=0; d<FESpaceVel::RefFE_T::dim; ++d)
        {
          id_T const dofId = feSpaceVel.dof.getId(this->feSpace.curFE.elem->id, n, d);
          localVel[d] += vel[dofId] * this->feSpace.curFE.phi[q](n);
        }
      }
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::numDOFs,CurFE_T::numDOFs>(d * CurFE_T::numDOFs, d * CurFE_T::numDOFs) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            (this->feSpace.curFE.dphi[q] * localVel).transpose();
      }
    }
  }

  // LMat_T build(/*LMat_T & Ke*/) const override
  // {
  //   LMat_T Ke;
  //   this->build(Ke);
  //   return Ke;
  // }

  double const coeff;
  Vec const & vel;
  FESpaceVel & feSpaceVel;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyGrad: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyGrad(double const c,
                        FESpace1_T & fe1,
                        FESpace2_T & fe2,
                        AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, cmp),
    coeff(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace1_T::CurFE_T::QR_T,
          typename FESpace2_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    uint d=0;
    for (auto const c: this->comp)
    {
      auto Kec = Ke.template block<CurFE1_T::numDOFs,CurFE2_T::numDOFs>(d*CurFE1_T::numDOFs, 0);
      for(uint q=0; q<FESpace1_T::CurFE_T::QR_T::numPts; ++q)
      {
        Kec += coeff * this->feSpace1.curFE.JxW[q] *
               this->feSpace1.curFE.dphi[q].col(c)*
               this->feSpace2.curFE.phi[q].transpose();
      }
      d++;
    }
  }

  double const coeff;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyDiv: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyDiv(double const c,
                       FESpace1_T & fe1,
                       FESpace2_T & fe2,
                       AssemblyBase::CompList const & cmp = allComp<FESpace2_T>()):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, cmp),
    coeff(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace1_T::CurFE_T::QR_T,
          typename FESpace2_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    uint d = 0;
    for (auto const c: this->comp)
    {
      auto Kec = Ke.template block<CurFE1_T::numDOFs,CurFE2_T::numDOFs>(0, d*CurFE2_T::numDOFs);
      for (uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
      {
        Kec += coeff * this->feSpace1.curFE.JxW[q] *
               this->feSpace1.curFE.phi[q] *
               this->feSpace2.curFE.dphi[q].col(c).transpose();
      }
      d++;
    }
  }

  double const coeff;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyVectorGrad: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyVectorGrad(double const c,
                              FESpace1_T & fe1,
                              FESpace2_T & fe2,
                              AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, cmp),
    coeff(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace1_T::CurFE_T::QR_T,
          typename FESpace2_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    // using CurFE2_T = typename FESpace2_T::CurFE_T;
    for(uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
    {
      Ke += coeff * this->feSpace1.curFE.JxW[q] *
          this->feSpace1.curFE.divphi[q] *
          this->feSpace2.curFE.phi[q].transpose();
    }
  }

  double const coeff;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyVectorDiv: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyVectorDiv(double const c,
                             FESpace1_T & fe1,
                             FESpace2_T & fe2,
                             AssemblyBase::CompList const & cmp = allComp<FESpace2_T>()):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, cmp),
    coeff(c)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace1_T::CurFE_T::QR_T,
          typename FESpace2_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    // using CurFE2_T = typename FESpace2_T::CurFE_T;
    for (uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
    {
      Ke += coeff * this->feSpace1.curFE.JxW[q] *
          this->feSpace1.curFE.phi[q] *
          this->feSpace2.curFE.divphi[q].transpose();
    }
  }

  double const coeff;
};

template <typename FESpace, typename FESpaceRhs = FESpace>
struct AssemblyS2SProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyS2SProjection(double const c,
                              Vec const & r,
                              FESpace_T & fe,
                              AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
    AssemblyVector<FESpace_T>(fe, cmp),
    coef(c),
    rhs(r),
    feSpaceRhs(fe)
  {}

  explicit AssemblyS2SProjection(double const c,
                                 Vec const & r,
                                 FESpace_T & fe,
                                 FESpaceRhs_T & feRhs,
                                 AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
    AssemblyVector<FESpace_T>(fe, cmp),
    coef(c),
    rhs(r),
    feSpaceRhs(feRhs)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace_T::CurFE_T::QR_T,
          typename FESpaceRhs_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    // if constexpr (!std::is_same<FESpace_T,FESpaceRhs_T>::value)
    {
      feSpaceRhs.curFE.reinit(elem);
    }
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFERhs_T = typename FESpaceRhs_T::CurFE_T;
    uint d = 0;
    for (auto const c: this->comp)
    {
      FVec<CurFERhs_T::numDOFs> localRhs;
      for(uint n=0; n<CurFERhs_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = feSpaceRhs.dof.getId(feSpaceRhs.curFE.elem->id, n, c);
        localRhs[n] = rhs[dofId];
      }
      auto Fec = Fe.template block<CurFE_T::numDOFs,1>(d*CurFE_T::numDOFs, 0);
      for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
      {
        Fec += coef * this->feSpace.curFE.JxW[q] *
               this->feSpace.curFE.phi[q] *
               (feSpaceRhs.curFE.phi[q].dot(localRhs));
      }
      d++;
    }
  }

  double coef;
  Vec const & rhs;
  FESpaceRhs_T & feSpaceRhs;
};

template <typename FESpace, typename FESpaceRhs>
struct AssemblyS2VProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyS2VProjection(double const c,
                                 Vec const & r,
                                 FESpace_T & fe,
                                 FESpaceRhs_T & feRhs):
    AssemblyVector<FESpace_T>(fe, {0}),
    coef(c),
    rhs(r),
    feSpaceRhs(feRhs)
  {
    static_assert(FESpace_T::RefFE_T::dim == FESpaceRhs_T::dim,
                  "the two fespaces are not of the same dimension");
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace_T::CurFE_T::QR_T,
          typename FESpaceRhs_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    feSpaceRhs.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFERhs_T = typename FESpaceRhs_T::CurFE_T;
    for (uint d=0; d<CurFE_T::RefFE_T::dim; ++d)
    {
      FVec<CurFERhs_T::numDOFs> localRhs;
      for(uint n=0; n<CurFERhs_T::numDOFs; ++n)
      {
        id_T const dofId = feSpaceRhs.dof.getId(feSpaceRhs.curFE.elem->id, n, d);
        localRhs[n] = rhs[dofId];
      }
      for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
      {
        Fe += coef * this->feSpace.curFE.JxW[q] *
              this->feSpace.curFE.phiVect[q].col(d) *
              (feSpaceRhs.curFE.phi[q].dot(localRhs)); // u*[q] = sum_k u*_k phi_k[q]
      }
    }
  }

  double const coef;
  Vec const & rhs;
  FESpaceRhs_T & feSpaceRhs;
};

template <typename FESpace, typename FESpaceRhs>
struct AssemblyV2SProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyV2SProjection(double const c,
                                 Vec const & r,
                                 FESpace_T & fe,
                                 FESpaceRhs_T & feRhs):
    AssemblyVector<FESpace_T>(fe, {0}),
    coef(c),
    rhs(r),
    feSpaceRhs(feRhs)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace_T::CurFE_T::QR_T,
          typename FESpaceRhs_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    feSpaceRhs.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFERhs_T = typename FESpaceRhs_T::CurFE_T;
    FVec<CurFERhs_T::numDOFs> localRhs;
    for(uint n=0; n<CurFERhs_T::RefFE_T::numFuns; ++n)
    {
      id_T const dofId = feSpaceRhs.dof.getId(feSpaceRhs.curFE.elem->id, n);
      localRhs[n] = rhs[dofId];
    }
    for (uint d=0; d<CurFE_T::RefFE_T::dim; ++d)
    {
      for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
      {
        Fe.template block<CurFE_T::numDOFs,1>(d*CurFE_T::numDOFs, 0) +=
            coef * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            (feSpaceRhs.curFE.phiVect[q].col(d).dot(localRhs)); // u*[q] = sum_k u*_k vec{psi}_k[q]
      }
    }
  }

  double const coef;
  Vec const & rhs;
  FESpaceRhs_T & feSpaceRhs;
};

template <typename FESpace, typename FESpaceRhs = FESpace>
struct AssemblyProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyProjection(double const c,
                              Vec const & r,
                              FESpace_T & fe,
                              AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
    AssemblyVector<FESpace_T>(fe, cmp)
  {
    if constexpr (FEDim<typename FESpace_T::RefFE_T>::value == FEDimType::SCALAR)
    {
      assembly = std::make_unique<AssemblyS2SProjection<FESpace_T>>(
                AssemblyS2SProjection<FESpace_T>{c, r, fe, cmp});
    }
    else
    {
      // no support for vector -> vector projection
      static_assert (dependent_false<FESpace>::value, "vector -> vector projection not yet implemented");
    }
  }

  explicit AssemblyProjection(double const c,
                              Vec const & r,
                              FESpace_T & fe,
                              FESpaceRhs_T & feRhs):
    AssemblyVector<FESpace_T>(fe, allComp<FESpace_T>())
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace_T::CurFE_T::QR_T,
          typename FESpaceRhs_T::CurFE_T::QR_T>::value,
          "the two quad rules are not the same");

    if constexpr (FEDim<typename FESpace_T::RefFE_T>::value == FEDimType::SCALAR &&
                  FEDim<typename FESpaceRhs_T::RefFE_T>::value == FEDimType::SCALAR)
    {
      assembly = std::make_unique<AssemblyS2SProjection<FESpace_T, FESpaceRhs_T>>(
              AssemblyS2SProjection<FESpace_T, FESpaceRhs_T>{c, r, fe, feRhs, allComp<FESpace_T>()});
    }
    else if constexpr (FEDim<typename FESpace_T::RefFE_T>::value == FEDimType::VECTOR &&
                       FEDim<typename FESpaceRhs_T::RefFE_T>::value == FEDimType::SCALAR)
    {
      assembly = std::make_unique<AssemblyS2VProjection<FESpace_T, FESpaceRhs_T>>(
              AssemblyS2VProjection<FESpace_T, FESpaceRhs_T>{c, r, fe, feRhs});
    }
    else if constexpr (FEDim<typename FESpace_T::RefFE_T>::value == FEDimType::SCALAR &&
                       FEDim<typename FESpaceRhs_T::RefFE_T>::value == FEDimType::VECTOR)
    {
      // static_assert (dependent_false<FESpace>::value, "not yet implemented");
      assembly = std::make_unique<AssemblyV2SProjection<FESpace_T, FESpaceRhs_T>>(
              AssemblyV2SProjection<FESpace_T, FESpaceRhs_T>{c, r, fe, feRhs});
    }
    else if constexpr (FEDim<typename FESpace_T::RefFE_T>::value == FEDimType::VECTOR &&
                       FEDim<typename FESpaceRhs_T::RefFE_T>::value == FEDimType::VECTOR)
    {
      // no support for vector -> vector projection
      static_assert (dependent_false<FESpace>::value, "vector -> vector projection not yet implemented");
    }
  }

  void reinit(GeoElem const & elem) const override
  {
    assembly->reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    assembly->build(Fe);
  }

  // LVec_T build() const
  // {
  //   return LVec_T{};
  // }

  // TODO: check own to avoid copies of the class
  // that would be implicitly deleted by the use of a unique_ptr
  std::shared_ptr<AssemblyVector<FESpace_T>> assembly;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyDivRhs: public AssemblyVector<FESpace1>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = AssemblyVector<FESpace1_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyDivRhs(double const c,
                          Vec const & vec,
                          FESpace1_T & fe1,
                          FESpace2_T & fe2,
                          AssemblyBase::CompList const & cmp = allComp<FESpace2_T>()):
    AssemblyVector<FESpace1_T>(fe1, cmp),
    coeff(c),
    feSpace2(fe2),
    data(vec)
  {}

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    feSpace2.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    FVec<CurFE2_T::numDOFs> localData;

    for (auto const c: this->comp)
    {
      for(uint n=0; n<CurFE2_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = feSpace2.dof.getId(feSpace2.curFE.elem->id, n, c);
        localData[n] = data[dofId];
      }

      for (uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
      {
        Fe += coeff * this->feSpace.curFE.JxW[q] *
              this->feSpace.curFE.phi[q] *
              (feSpace2.curFE.dphi[q].col(c).dot(localData));
      }
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coeff;
  FESpace2 & feSpace2;
  Vec const & data;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyGradRhs: public AssemblyVector<FESpace1>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = AssemblyVector<FESpace1_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyGradRhs(double const c,
                          Vec const & vec,
                          FESpace1_T & fe1,
                          FESpace2_T const & fe2,
                          AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
    AssemblyVector<FESpace1_T>(fe1, cmp),
    coeff(c),
    feSpace2(fe2),
    data(vec)
  {}

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    feSpace2.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    FVec<CurFE2_T::numDOFs * FESpace2_T::dim> localData;

    for (uint d2=0; d2<FESpace2_T::dim; ++d2)
    {
      for (uint n=0; n<CurFE2_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = feSpace2.dof.getId(feSpace2.curFE.elem->id, n, d2);
        localData[n + d2*CurFE2_T::numDOFs] = data[dofId];
      }
    }

    uint counter = 0;
    for (auto const c: this->comp)
    {
      uint const d1 = c % FESpace1_T::Mesh_T::Elem_T::dim;
      uint const d2 = c / FESpace1_T::Mesh_T::Elem_T::dim;
      for (uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
      {
        // d u_d2 / d x_d1
        double const gradDataQ =
            feSpace2.curFE.dphi[q].col(d1).dot(
              localData.template block<CurFE2_T::numDOFs, 1>(d2 * CurFE2_T::numDOFs, 0));
        Fe.template block<CurFE1_T::numDOFs, 1>(counter * CurFE1_T::numDOFs, 0) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            gradDataQ;
        // // u_d2
        // double const dataQ =
        //     feSpace2.curFE.phi[q].dot(
        //       localData.template block<CurFE2_T::numDOFs, 1>(d2 * CurFE2_T::numDOFs, 0));
        // Fe.template block<CurFE1_T::numDOFs, 1>(counter * CurFE1_T::numDOFs, 0) +=
        //     coeff * this->feSpace.curFE.JxW[q] *
        //     this->feSpace.curFE.dphi[q].col(d1) *
        //     dataQ;
      }
      counter++;
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coeff;
  FESpace2_T const & feSpace2;
  Vec const & data;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyGradRhs2: public AssemblyVector<FESpace1>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = AssemblyVector<FESpace1_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyGradRhs2(double const c,
                            Vec const & vec,
                            FESpace1_T & fe1,
                            FESpace2_T & fe2,
                            AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
    AssemblyVector<FESpace1_T>(fe1, cmp),
    coeff(c),
    feSpace2(fe2),
    data(vec)
  {}

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    feSpace2.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    FVec<CurFE2_T::numDOFs * FESpace2_T::dim> localData;

    for (uint d2=0; d2<FESpace2_T::dim; ++d2)
    {
      for (uint n=0; n<CurFE2_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = feSpace2.dof.getId(feSpace2.curFE.elem->id, n, d2);
        localData[n + d2*CurFE2_T::numDOFs] = data[dofId];
      }
    }

    uint counter = 0;
    for (auto const c: this->comp)
    {
      uint const d1 = c % FESpace1_T::Mesh_T::Elem_T::dim;
      uint const d2 = c / FESpace1_T::Mesh_T::Elem_T::dim;
      for (uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
      {
        // // d u_d2 / d x_d1
        // double const gradDataQ =
        //     feSpace2.curFE.dphi[q].col(d1).dot(
        //       localData.template block<CurFE2_T::numDOFs, 1>(d2 * CurFE2_T::numDOFs, 0));
        // Fe.template block<CurFE1_T::numDOFs, 1>(counter * CurFE1_T::numDOFs, 0) +=
        //     coeff * this->feSpace.curFE.JxW[q] *
        //     this->feSpace.curFE.phi[q] *
        //     gradDataQ;
        // u_d2
        double const dataQ =
            feSpace2.curFE.phi[q].dot(
              localData.template block<CurFE2_T::numDOFs, 1>(d2 * CurFE2_T::numDOFs, 0));
        Fe.template block<CurFE1_T::numDOFs, 1>(counter * CurFE1_T::numDOFs, 0) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.dphi[q].col(d1) *
            dataQ;
      }
      counter++;
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coeff;
  FESpace2_T & feSpace2;
  Vec const & data;
};

template <typename FESpace1, typename FESpace2 = FESpace1>
struct AssemblyStiffnessRhs: public AssemblyVector<FESpace1>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = AssemblyVector<FESpace1_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyStiffnessRhs(double const c,
                                Vec const & vec,
                                FESpace1_T & fe1,
                                FESpace2_T & fe2,
                                AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
    AssemblyVector<FESpace1_T>(fe1, cmp),
    coeff(c),
    feSpace2(fe2),
    data(vec)
  {}

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    feSpace2.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    FVec<CurFE2_T::numDOFs> localData;

    for(uint n=0; n<CurFE2_T::RefFE_T::numFuns; ++n)
    {
      id_T const dofId = feSpace2.dof.getId(feSpace2.curFE.elem->id, n);
      localData[n] = data[dofId];
    }
    uint d = 0;
    for (auto const c: this->comp)
    {
      for (uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
      {
        Fe.template block<CurFE1_T::numDOFs,1>(d*CurFE1_T::numDOFs, 0) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.dphi[q].col(c) *
            (feSpace2.curFE.dphi[q].col(c).dot(localData));
      }
      d++;
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coeff;
  FESpace2_T & feSpace2;
  Vec const & data;
};

template <typename FESpace1, typename FESpace2 = FESpace1>
struct AssemblyAdvectionRhs: public AssemblyVector<FESpace1>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = AssemblyVector<FESpace1_T>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyAdvectionRhs(double const c,
                                Vec const & u,
                                Vec const & vec,
                                FESpace1_T & fe1,
                                FESpace2_T & fe2,
                                AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
    AssemblyVector<FESpace1_T>(fe1, cmp),
    coeff(c),
    vel(u),
    feSpace2(fe2),
    data(vec)
  {}

  explicit AssemblyAdvectionRhs(double const c,
                                Vec const & u,
                                Vec const & vec,
                                FESpace1_T & fe1,
                                AssemblyBase::CompList const & cmp = allComp<FESpace1_T>()):
    AssemblyVector<FESpace1_T>(fe1, cmp),
    coeff(c),
    vel(u),
    feSpace2(fe1),
    data(vec)
  {}

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    feSpace2.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    FVec<CurFE2_T::numDOFs> localData;

    for(uint n=0; n<CurFE2_T::RefFE_T::numFuns; ++n)
    {
      id_T const dofId = feSpace2.dof.getId(feSpace2.curFE.elem->id, n);
      localData[n] = data[dofId];
    }
    uint d = 0;
    for (auto const c: this->comp)
    {
      for (uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
      {
        FVec<3> localVel = FVec<3>::Zero();
        for(uint n=0; n<CurFE1_T::RefFE_T::numFuns; ++n)
        {
          for (uint d=0; d<FESpace1_T::RefFE_T::dim; ++d)
          {
            id_T const dofId = this->feSpace.dof.getId(this->feSpace.curFE.elem->id, n, d);
            localVel[d] += vel[dofId] * this->feSpace.curFE.phi[q](n);
          }
        }

        Fe.template block<CurFE1_T::numDOFs,1>(d*CurFE1_T::numDOFs, 0) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            (feSpace2.curFE.dphi[q].col(c).dot(localData)) * localVel[d];
      }
      d++;
    }
  }

  LVec_T build(/*LVec_T & Fe*/) const
  {
    LVec_T Fe;
    this->build(Fe);
    return Fe;
  }

  double const coeff;
  Vec const & vel;
  FESpace2_T & feSpace2;
  Vec const & data;
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

  AssemblyBCNatural(Fun<FESpace_T::dim, 3> const r,
                    marker_T const m,
                    FESpace_T & fe,
                    AssemblyBase::CompList const cmp = allComp<FESpace_T>()):
    AssemblyVector<FESpace_T>(fe, std::move(cmp)),
    rhs(std::move(r)),
    marker(m)
  {}

  AssemblyBCNatural(scalarFun_T const r,
                    marker_T const m,
                    FESpace_T & fe,
                    AssemblyBase::CompList const cmp = allComp<FESpace_T>()):
    AssemblyBCNatural<FESpace_T>([r](Vec3 const & p) { return Vec1::Constant(r(p)); }, m, fe, cmp)
  {
    static_assert(FESpace_T::dim == 1);
  }

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    auto const & mesh = this->feSpace.mesh;
    auto const & e = *this->feSpace.curFE.elem;
    uint facetCounter = 0;
    for(auto const facetId: mesh.elemToFacet[e.id])
    {
      if(facetId != dofIdNotSet &&
         mesh.facetList[facetId].marker == marker)
      {
        auto const & facet = mesh.facetList[facetId];
        facetCurFE.reinit(facet);
        for(uint q=0; q<QR_T::numPts; ++q)
        {
          auto const value = rhs(facetCurFE.qpoint[q]);
          for (uint const d: this->comp)
          {
            for(uint i=0; i<FacetFE_T::numFuns; ++i)
            {
              auto const id = CurFE_T::RefFE_T::dofOnFacet[facetCounter][i] + d*CurFE_T::numDOFs;
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

  AssemblyBCNormal(scalarFun_T const r,
                   marker_T const m,
                   FESpace & fe,
                   AssemblyBase::CompList const cmp = allComp<FESpace>()):
    AssemblyVector<FESpace>(fe, std::move(cmp)),
    rhs(std::move(r)),
    marker(m)
  {}

  void build(LVec_T & Fe) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    auto const & mesh = this->feSpace.mesh;
    auto const & e = *this->feSpace.curFE.elem;
    uint facetCounter = 0;
    for(auto const facetId: mesh.elemToFacet[e.id])
    {
      if(facetId != dofIdNotSet &&
         mesh.facetList[facetId].marker == marker)
      {
        auto const & facet = mesh.facetList[facetId];
        auto const normal = facet.normal();
        facetCurFE.reinit(facet);
        for(uint q=0; q<QR_T::numPts; ++q)
        {
          for (auto const d: this->comp)
          {
            auto const value = rhs(facetCurFE.qpoint[q]);
            for(uint i=0; i<FacetFE_T::numFuns; ++i)
            {
              auto const id = CurFE_T::RefFE_T::dofOnFacet[facetCounter][i] + d*CurFE_T::numDOFs;
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

  AssemblyBCMixed(scalarFun_T const c,
                  marker_T const m,
                  FESpace & fe,
                  AssemblyBase::CompList const cmp = allComp<FESpace>()):
    Diagonal<FESpace>(fe, std::move(cmp)),
    coef(std::move(c)),
    marker(m)
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    auto const & mesh = this->feSpace.mesh;
    uint facetCounter = 0;
    for(auto const facetId: mesh.elemToFacet[this->feSpace.curFE.elem->id])
    {
      if(facetId != dofIdNotSet &&
         mesh.facetList[facetId].marker == marker)
      {
        auto const & facet = mesh.facetList[facetId];
        facetCurFE.reinit(facet);
        for(uint q=0; q<QR_T::numPts; ++q)
        {
          auto const localCoef = coef(facetCurFE.qpoint[q]);
          for (uint const d: this->comp)
          {
            for(uint i=0; i<FacetFE_T::numFuns; ++i)
            {
              auto const idI = CurFE_T::RefFE_T::dofOnFacet[facetCounter][i] + d*CurFE_T::numDOFs;
              for(uint j=0; j<FacetFE_T::numFuns; ++j)
              {
                auto const idJ = CurFE_T::RefFE_T::dofOnFacet[facetCounter][j] + d*CurFE_T::numDOFs;
                Ke(idI, idJ) +=
                    facetCurFE.JxW[q] * localCoef *
                    facetCurFE.phi[q](i) *
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
