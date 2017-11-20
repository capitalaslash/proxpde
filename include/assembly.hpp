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

  virtual void reinit(GeoElem const & elem) final
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

  virtual void reinit(GeoElem const & elem)
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
  using CurFE_T = typename FESpace_T::CurFE_T;
  using LMat_T = FMat<FESpace::dim*CurFE_T::size,FESpace::dim*CurFE_T::size>;
  using LVec_T = FVec<FESpace::dim*CurFE_T::size>;

  explicit AssemblyVector(FESpace_T & fe, uint offset_row, CompList const & comp):
    AssemblyBase{offset_row, 0, comp},
    feSpace(fe)
  {}

  virtual void build(LVec_T & Fe) const = 0;

  virtual void reinit(GeoElem const & elem)
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
                    AssemblyBase::CompList const & comp = {0},
                    uint offset_row = 0,
                    uint offset_clm = 0):
    Diagonal<FESpace_T>(fe, offset_row, offset_clm, comp),
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
            this->feSpace.curFE.dphi[q] *
            this->feSpace.curFE.dphi[q].transpose();
      }
    }
  }

  LMat_T build(/*LMat_T & Ke*/) const
  {
    LMat_T Ke;
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<CurFE_T::size,CurFE_T::size>(d*CurFE_T::size, d*CurFE_T::size) +=
            coeff * this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.dphi[q] *
            this->feSpace.curFE.dphi[q].transpose();
      }
    }
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
                        AssemblyBase::CompList const & comp = {0},
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
      Ke += coeff * this->feSpace.curFE.JxW[q] *
        this->feSpace.curFE.phi[q] *
        this->feSpace.curFE.phi[q].transpose();
    }
  }

  LMat_T build(/*LMat_T & Ke*/) const
  {
    LMat_T Ke;
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Ke += coeff * this->feSpace.curFE.JxW[q] *
          this->feSpace.curFE.phi[q] *
          this->feSpace.curFE.phi[q].transpose();
    }
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
                      AssemblyBase::CompList const & comp = {0},
                      uint offset_row = 0):
    AssemblyVector<FESpace>(fe, offset_row, comp),
    rhs(r)
  {}

  // this constructor is available only when FESpace::dim == 1
  template <typename FESpace1 = FESpace,
            std::enable_if_t<FESpace1::dim == 1, bool> = true>
  AssemblyAnalyticRhs(scalarFun_T const & r,
                      FESpace & fe,
                      AssemblyBase::CompList const & comp = {0},
                      uint offset_row = 0):
    AssemblyAnalyticRhs<FESpace>([r](Vec3 const & p) {return Vec1(r(p));}, fe, comp, offset_row)
  {}

  // otherwise, we fail
  template <typename FESpace1 = FESpace,
            std::enable_if_t<FESpace1::dim != 1, bool> = true>
  AssemblyAnalyticRhs(scalarFun_T const & r,
                      FESpace & fe,
                      AssemblyBase::CompList const & comp = {0},
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
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyVecRhs(Vec const & r,
                          FESpace_T & fe,
                          AssemblyBase::CompList const & comp = {0},
                          uint offset_row = 0):
    AssemblyVector<FESpace_T>(fe, offset_row, comp),
    rhs(r)
  {}

  void build(LVec_T & Fe) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double local_rhs = 0.0;
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = this->feSpace.dof.elemMap[this->feSpace.curFE.e->id][n];
        local_rhs += rhs(dofId) * this->feSpace.curFE.phi[q](n);
      }
      Fe += this->feSpace.curFE.JxW[q] *
          this->feSpace.curFE.phi[q] *
          local_rhs;
    }
  }

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

  explicit AssemblyAdvection(Field3 const u,
                             FESpace_T & fe,
                             AssemblyBase::CompList const & comp = {0},
                             uint offset_row = 0,
                             uint offset_clm = 0):
    Diagonal<FESpace>(fe, offset_row, offset_clm, comp),
    vel(u)
  {}

  void build(LMat_T & Ke) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Vec3 local_vel = Vec3::Zero();
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = this->feSpace.dof.elemMap[this->feSpace.curFE.e->id][n];
        local_vel += vel.row(dofId) * this->feSpace.curFE.phi[q](n);
      }
      Ke += this->feSpace.curFE.JxW[q] *
          this->feSpace.curFE.phi[q] *
          (this->feSpace.curFE.dphi[q] * local_vel).transpose();
      // for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      // {
      //   for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
      //   {
      //     Ke(i,j) += this->feSpace.curFE.JxW[q] *
      //         (local_vel.dot(this->feSpace.curFE.dphi[q].row(i).transpose())) *
      //         this->feSpace.curFE.phi[q](j);
      //   }
      // }
    }
  }

  LMat_T build(/*LMat_T & Ke*/) const
  {
    LMat_T Ke;
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Vec3 local_vel = Vec3::Zero();
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = this->feSpace.dof.elemMap[this->feSpace.curFE.e->id][n];
        local_vel += vel.row(dofId) * this->feSpace.curFE.phi[q](n);
      }
      Ke += this->feSpace.curFE.JxW[q] *
          this->feSpace.curFE.phi[q] *
          (this->feSpace.curFE.dphi[q] * local_vel).transpose();
      // for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      // {
      //   for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
      //   {
      //     Ke(i,j) += this->feSpace.curFE.JxW[q] *
      //         (local_vel.dot(this->feSpace.curFE.dphi[q].row(i).transpose())) *
      //         this->feSpace.curFE.phi[q](j);
      //   }
      // }
    }
    return Ke;
  }

  Field3 vel;
};

template <typename FESpace>
std::vector<uint> allComp()
{
  std::vector<uint> comp(FESpace::dim);
  std::iota(comp.begin(), comp.end(), 0);
  return comp;
}

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
