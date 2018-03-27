#pragma once

#include "def.hpp"
#include "blockmatrix.hpp"
#include "fespace.hpp"

#include <type_traits>

struct AssemblyBase
{
  using CompList = std::vector<uint>;

  uint const offset_row;
  uint const offset_clm;
  CompList const comp;
};

template <typename FESpace>
struct Diagonal: public AssemblyBase
{
  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace_T::CurFE_T;
  using LMat_T = FMat<FESpace_T::dim*CurFE_T::size,FESpace_T::dim*CurFE_T::size>;
  using LVec_T = FVec<FESpace_T::dim*CurFE_T::size>;

  explicit Diagonal(FESpace_T & fe, uint offset_row, uint offset_clm, CompList const & comp):
    AssemblyBase{offset_row, offset_clm, comp},
    feSpace(fe)
  {}

  virtual void build(LMat_T & Ke) const = 0;
  virtual LMat_T build(/*LMat_T & Ke*/) const = 0;

  virtual void reinit(GeoElem const & elem) const final
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
  using LMat_T = FMat<FESpace1_T::dim*CurFE1_T::size,FESpace2_T::dim*CurFE2_T::size>;
  using LVec_T = FVec<FESpace1_T::dim*CurFE1_T::size>;

  explicit Coupling(FESpace1_T & fe1,
                    FESpace2 & fe2,
                    uint offset_row,
                    uint offset_clm,
                    CompList const & comp):
    AssemblyBase{offset_row, offset_clm, comp},
    feSpace1(fe1),
    feSpace2(fe2)
  {}

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
  using LMat_T = FMat<FESpace_T::dim*CurFE_T::size, FESpace_T::dim*CurFE_T::size>;
  using LVec_T = FVec<FESpace_T::dim*CurFE_T::size>;

  explicit AssemblyVector(FESpace_T & fe, uint offset_row, CompList const & comp):
    AssemblyBase{offset_row, 0, comp},
    feSpace(fe)
  {}

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
                    AssemblyBase::CompList const & comp = allComp<FESpace>(),
                    uint offset_row = 0,
                    uint offset_clm = 0):
    Diagonal<FESpace_T>(fe, offset_row, offset_clm, comp),
    coeff(c)
  {}

  void build(LMat_T & Ke) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::size,CurFE_T::size>(d*CurFE_T::size, d*CurFE_T::size) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.dphi[q] *
            this->feSpace.curFE.dphi[q].transpose();
      }
    }
  }

  LMat_T build(/*LMat_T & Ke*/) const
  {
    LMat_T Ke;
    this->build(Ke);
    return Ke;
  }

  double coeff;
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
                    AssemblyBase::CompList const & comp = allComp<FESpace>(),
                    uint offset_row = 0,
                    uint offset_clm = 0):
    Diagonal<FESpace_T>(fe, offset_row, offset_clm, comp),
    coeff(c)
  {}

  void build(LMat_T & Ke) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint di=0; di<FESpace_T::dim; ++di)
      {
        // d u_i / d x_j * d \phi_i / d x_j
        Ke.template block<CurFE_T::size,CurFE_T::size>(di*CurFE_T::size, di*CurFE_T::size) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.dphi[q] *
            this->feSpace.curFE.dphi[q].transpose();
        for (uint dj=0; dj<FESpace_T::dim; ++dj)
        {
          // d u_j / d x_i * d \phi_i / d x_j
          Ke.template block<CurFE_T::size,CurFE_T::size>(di*CurFE_T::size, dj*CurFE_T::size) +=
              coeff * this->feSpace.curFE.JxW[q] *
              this->feSpace.curFE.dphi[q].col(di) *
              this->feSpace.curFE.dphi[q].col(dj).transpose();
//          for (uint i=0; i<CurFE_T::size; ++i)
//            for (uint j=0; j<CurFE_T::size; ++j)
//            {
//              // d u_j / d x_i * d \phi_i / d x_j
//              Ke(i+di*CurFE_T::size, j+dj*CurFE_T::size) +=
//                  coeff * this->feSpace.curFE.JxW[q] *
//                  this->feSpace.curFE.dphi[q].row(j)(di) *
//                  this->feSpace.curFE.dphi[q].row(i)(dj);
//            }
        }
      }
    }
  }

  LMat_T build(/*LMat_T & Ke*/) const
  {
    LMat_T Ke;
    this->build(Ke);
    return Ke;
  }

  double coeff;
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
                        AssemblyBase::CompList const & comp = allComp<FESpace>(),
                        uint offset_row = 0,
                        uint offset_clm = 0):
     Diagonal<FESpace>(fe, offset_row, offset_clm, comp),
     coeff(c)
  {}

  void build(LMat_T & Ke) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::size,CurFE_T::size>(d*CurFE_T::size, d*CurFE_T::size) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            this->feSpace.curFE.phi[q].transpose();
      }
    }
  }

  LMat_T build(/*LMat_T & Ke*/) const
  {
    LMat_T Ke;
    this->build(Ke);
    return Ke;
  }

  double coeff;
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
                              AssemblyBase::CompList const & comp = allComp<FESpace>(),
                              uint offset_row = 0,
                              uint offset_clm = 0):
     Diagonal<FESpace>(fe, offset_row, offset_clm, comp),
     coeff(c)
  {}

  void build(LMat_T & Ke) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::size,CurFE_T::size>(d*CurFE_T::size, d*CurFE_T::size) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phiVect[q] *
            this->feSpace.curFE.phiVect[q].transpose();
      }
    }
  }

  LMat_T build(/*LMat_T & Ke*/) const
  {
    LMat_T Ke;
    this->build(Ke);
    return Ke;
  }

  double coeff;
};

template <typename FESpace>
struct AssemblyAnalyticRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyAnalyticRhs(Fun<FESpace_T::dim,3> const & r,
                      FESpace & fe,
                      AssemblyBase::CompList const & comp = allComp<FESpace>(),
                      uint offset_row = 0):
    AssemblyVector<FESpace>(fe, offset_row, comp),
    rhs(r)
  {}

  // this constructor is available only when FESpace::dim == 1
  template <typename FESpace1 = FESpace,
            std::enable_if_t<FESpace1::dim == 1, bool> = true>
  AssemblyAnalyticRhs(scalarFun_T const & r,
                      FESpace & fe,
                      AssemblyBase::CompList const & comp = allComp<FESpace>(),
                      uint offset_row = 0):
    AssemblyAnalyticRhs<FESpace>([r](Vec3 const & p) {return Vec1(r(p));}, fe, comp, offset_row)
  {}

  // otherwise, we fail
  template <typename FESpace1 = FESpace,
            std::enable_if_t<FESpace1::dim != 1, bool> = true>
  AssemblyAnalyticRhs(scalarFun_T const &,
                      FESpace & fe,
                      AssemblyBase::CompList const & comp = allComp<FESpace>(),
                      uint offset_row = 0):
    AssemblyVector<FESpace>(fe, offset_row, comp)
  {
    std::abort();
  }

  void build(LVec_T & Fe) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Fe.template block<CurFE_T::size,1>(d*CurFE_T::size, 0) +=
            this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            rhs(this->feSpace.curFE.qpoint[q])(d);
      }
    }
  }

  Fun<FESpace_T::dim,3> const rhs;
};

template <typename FESpace>
struct AssemblyVecRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyVecRhs(double const c,
                          Vec const & r,
                          FESpace_T & fe,
                          AssemblyBase::CompList const & comp = allComp<FESpace>(),
                          uint offset_row = 0):
    AssemblyVector<FESpace_T>(fe, offset_row, comp),
    coeff(c),
    rhs(r)
  {}

  void build(LVec_T & Fe) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    FMat<CurFE_T::size,FESpace::dim> localRhs;
    for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
    {
      id_T const dofId = this->feSpace.dof.elemMap[this->feSpace.curFE.e->id][n];
      uint d = 0;
      for (auto const c: this->comp)
      {
        localRhs(n, d) = rhs[dofId + c*this->feSpace.dof.totalNum];
        d++;
      }
    }
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<this->comp.size(); ++d)
      {
        Fe.template block<CurFE_T::size,1>(d*CurFE_T::size, 0) +=
            this->feSpace.curFE.JxW[q] * coeff *
            this->feSpace.curFE.phi[q] *
            (localRhs.col(d).dot(this->feSpace.curFE.phi[q])); // rhs[q] = sum_k u_k * phi_k[q]
      }
    }
  }

  double const coeff;
  Vec const & rhs;
};

template <typename FESpace>
struct AssemblyAdvection: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using VelVec_T = Eigen::Matrix<double, Eigen::Dynamic, FESpace_T::RefFE_T::dim>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyAdvection(double const c,
                             Vec const & u,
                             FESpace_T & fe,
                             AssemblyBase::CompList const & comp = allComp<FESpace>(),
                             uint offset_row = 0,
                             uint offset_clm = 0):
    Diagonal<FESpace>(fe, offset_row, offset_clm, comp),
    coeff(c),
    vel(u)
  {}

  void build(LMat_T & Ke) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      FVec<3> localVel = Vec3::Zero();
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = this->feSpace.dof.elemMap[this->feSpace.curFE.e->id][n];
        for (uint d=0; d<FESpace_T::RefFE_T::dim; ++d)
          localVel[d] += vel[dofId + d*this->feSpace.dof.totalNum] * this->feSpace.curFE.phi[q](n);
      }
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::size,CurFE_T::size>(d*CurFE_T::size, d*CurFE_T::size) +=
            this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q] *
            (this->feSpace.curFE.dphi[q] * localVel).transpose();
      }
    }
  }

  LMat_T build(/*LMat_T & Ke*/) const
  {
    LMat_T Ke;
    this->build(Ke);
    return Ke;
  }

  double const coeff;
  Vec const & vel;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyGrad: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyGrad(FESpace1_T & fe1,
                        FESpace2_T & fe2,
                        std::vector<uint> comp = allComp<FESpace1_T>(),
                        uint offset_row = 0,
                        uint offset_clm = 0):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, offset_row, offset_clm, comp),
    component(comp)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace1_T::CurFE_T::QR_T,
          typename FESpace2_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    uint d=0;
    for (auto const c: component)
    {
      auto Kec = Ke.template block<CurFE1_T::size,CurFE2_T::size>(d*CurFE1_T::size, 0);
      for(uint q=0; q<FESpace1_T::CurFE_T::QR_T::numPts; ++q)
      {
        Kec += this->feSpace1.curFE.JxW[q] *
               this->feSpace1.curFE.dphi[q].col(c)*
               this->feSpace2.curFE.phi[q].transpose();
      }
      d++;
    }
  }

  std::vector<uint> const component;
};

template <typename FESpace1, typename FESpace2>
AssemblyGrad<FESpace1, FESpace2> make_assemblyGrad(
    FESpace1 & fe1,
    FESpace2 & fe2,
    std::vector<uint> comp = allComp<FESpace1>(),
    uint offset_row = 0,
    uint offset_clm = 0)
{
  return AssemblyGrad<FESpace1, FESpace2>(fe1, fe2, comp, offset_row, offset_clm);
}

template <typename FESpace1, typename FESpace2>
struct AssemblyDiv: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyDiv(FESpace1_T & fe1,
                       FESpace2_T & fe2,
                       std::vector<uint> comp = allComp<FESpace2_T>(),
                       uint offset_row = 0,
                       uint offset_clm = 0):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, offset_row, offset_clm, comp),
    component(comp)
  {
    // this works only if the same quad rule is defined on both CurFE
    static_assert(
          std::is_same<
          typename FESpace1_T::CurFE_T::QR_T,
          typename FESpace2_T::CurFE_T::QR_T>::value,
          "the two quad rule are not the same");
  }

  void build(LMat_T & Ke) const
  {
    using CurFE1_T = typename FESpace1_T::CurFE_T;
    using CurFE2_T = typename FESpace2_T::CurFE_T;
    uint d = 0;
    for (auto const c: component)
    {
      auto Kec = Ke.template block<CurFE1_T::size,CurFE2_T::size>(0, d*CurFE2_T::size);
      for (uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
      {
        Kec += this->feSpace1.curFE.JxW[q] *
               this->feSpace1.curFE.phi[q] *
               this->feSpace2.curFE.dphi[q].col(c).transpose();
      }
      d++;
    }
  }

  std::vector<uint> const component;
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
                              std::vector<uint> comp = allComp<FESpace_T>(),
                              uint offset_row = 0):
    AssemblyVector<FESpace_T>(fe, offset_row, comp),
    coef(c),
    rhs(r),
    feSpaceRhs(fe)
  {}

  explicit AssemblyProjection(double const c,
                              Vec const & r,
                              FESpace_T & fe,
                              FESpaceRhs_T & feRhs,
                              std::vector<uint> comp = allComp<FESpace_T>(),
                              uint offset_row = 0):
    AssemblyVector<FESpace_T>(fe, offset_row, comp),
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

  void reinit(GeoElem const & elem) const
  {
    this->feSpace.curFE.reinit(elem);
    feSpaceRhs.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFERhs_T = typename FESpaceRhs_T::CurFE_T;
    uint d = 0;
    for (auto const c: this->comp)
    {
      FVec<CurFERhs_T::size> localRhs;
      for(uint n=0; n<CurFERhs_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = feSpaceRhs.dof.elemMap[feSpaceRhs.curFE.e->id][n] + c*feSpaceRhs.dof.totalNum;
        localRhs[n] = rhs[dofId];
      }
      auto Fec = Fe.template block<CurFE_T::size,1>(d*CurFE_T::size, 0);
      for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
      {
        Fec += coef * this->feSpace.curFE.JxW[q] *
               this->feSpace.curFE.phi[q] *
               (feSpaceRhs.curFE.phi[q].dot(localRhs));
      }
      d++;
    }
  }

  double const coef;
  Vec const & rhs;
  FESpaceRhs_T & feSpaceRhs;
};

template <typename FESpace, typename FESpaceRhs>
struct AssemblyVectorProjection: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceRhs_T = FESpaceRhs;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyVectorProjection(double const c,
                              Vec const & r,
                              FESpace_T & fe,
                              FESpaceRhs_T & feRhs,
                              uint offset_row = 0):
    AssemblyVector<FESpace_T>(fe, offset_row, {0}),
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

  void reinit(GeoElem const & elem) const
  {
    this->feSpace.curFE.reinit(elem);
    feSpaceRhs.curFE.reinit(elem);
  }

  void build(LVec_T & Fe) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    using CurFERhs_T = typename FESpaceRhs_T::CurFE_T;
    for (uint d=0; d<CurFE_T::RefFE_T::dim; ++d)
    {
      FVec<CurFERhs_T::size> localRhs;
      for(uint n=0; n<CurFERhs_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = feSpaceRhs.dof.elemMap[feSpaceRhs.curFE.e->id][n] + d*feSpaceRhs.dof.totalNum;
        localRhs[n] = rhs[dofId];
      }
      for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
      {
        Fe += coef * this->feSpace.curFE.JxW[q] *
              this->feSpace.curFE.phiVect[q].col(d) *
              (feSpaceRhs.curFE.phi[q].dot(localRhs));
      }
    }
  }

  double const coef;
  Vec const & rhs;
  FESpaceRhs_T & feSpaceRhs;
};
