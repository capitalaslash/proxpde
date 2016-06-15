#pragma once

#include "def.hpp"
#include "fespace.hpp"

#include <type_traits>

struct AssemblyBase
{
  uint offset_row;
  uint offset_clm;
};

template <typename FESpace>
struct Diagonal: public AssemblyBase
{
  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace_T::CurFE_T;
  using LMat_T = FMat<CurFE_T::size(),CurFE_T::size()>;
  using LVec_T = FVec<CurFE_T::size()>;

  explicit Diagonal(FESpace_T & fe, uint offset_row, uint offset_clm):
    AssemblyBase{offset_row, offset_clm},
    feSpace(fe)
  {}

  virtual void build(LMat_T & Ke) const = 0;

  FESpace_T & feSpace;
};

template <typename FESpace1, typename FESpace2>
struct Coupling: public AssemblyBase
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using CurFE1_T = typename FESpace1_T::CurFE_T;
  using CurFE2_T = typename FESpace2_T::CurFE_T;
  using LMat_T = FMat<CurFE1_T::size(),CurFE2_T::size()>;
  using LVec_T = FVec<CurFE1_T::size()>;

  explicit Coupling(FESpace1_T & fe1,
                    FESpace2 & fe2,
                    uint offset_row,
                    uint offset_clm):
    AssemblyBase{offset_row, offset_clm},
    feSpace1(fe1),
    feSpace2(fe2)
  {}

  virtual void build(LMat_T & Ke) const = 0;

  FESpace1_T & feSpace1;
  FESpace2_T & feSpace2;
};

template <typename FESpace>
struct AssemblyVector: public AssemblyBase
{
  using FESpace_T = FESpace;
  using CurFE_T = typename FESpace_T::CurFE_T;
  using LMat_T = FMat<CurFE_T::size(),CurFE_T::size()>;
  using LVec_T = FVec<CurFE_T::size()>;

  explicit AssemblyVector(FESpace_T & fe, uint offset_row):
    AssemblyBase{offset_row, 0},
    feSpace(fe)
  {}

  virtual void build(LVec_T & Fe) const = 0;

  FESpace_T & feSpace;
};

template <typename FESpace>
struct AssemblyStiffness: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyStiffness(double const c,
                             FESpace_T & fe,
                             uint offset_row = 0,
                             uint offset_clm = 0):
    Diagonal<FESpace_T>(fe, offset_row, offset_clm),
    coeff(c)
  {}

  void build(LMat_T & Ke) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Ke += coeff * this->feSpace.curFE.JxW[q] *
          this->feSpace.curFE.dphi[q] *
          this->feSpace.curFE.dphi[q].transpose();
    }
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
                        uint offset_row = 0,
                        uint offset_clm = 0):
     Diagonal<FESpace>(fe, offset_row, offset_clm),
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

  double coeff;
};

template <typename FESpace>
struct AssemblyAnalyticRhs: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyAnalyticRhs(scalarFun_T const & r,
                               FESpace & fe,
                               uint offset_row = 0):
    AssemblyVector<FESpace>(fe, offset_row),
    rhs(r)
  {}

  void build(LVec_T & Fe) const
  {
    using CurFE_T = typename FESpace_T::CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Fe += this->feSpace.curFE.JxW[q] *
          this->feSpace.curFE.phi[q] *
          rhs(this->feSpace.curFE.qpoint[q]);
    }
  }

  scalarFun_T const & rhs;
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
                          uint offset_row = 0):
    AssemblyVector<FESpace_T>(fe, offset_row),
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

  explicit AssemblyAdvection(Vec3d const u,
                             FESpace_T & fe,
                             uint offset_row = 0,
                             uint offset_clm = 0):
    Diagonal<FESpace>(fe, offset_row, offset_clm),
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

  Vec3d vel;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyGrad: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyGrad(uint comp,
                        FESpace1_T & fe1,
                        FESpace2_T & fe2,
                        uint offset_row = 0,
                        uint offset_clm = 0):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, offset_row, offset_clm),
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
    for(uint q=0; q<FESpace1_T::CurFE_T::QR_T::numPts; ++q)
    {
      Ke -= this->feSpace1.curFE.JxW[q] *
          this->feSpace1.curFE.dphi[q].col(component)*
          this->feSpace2.curFE.phi[q].transpose();
    }
  }
  uint const component;
};

template <typename FESpace1, typename FESpace2>
AssemblyGrad<FESpace1, FESpace2> make_assemblyGrad(uint comp, FESpace1 & fe1, FESpace2 & fe2, uint offset_row = 0, uint offset_clm = 0)
{
  return AssemblyGrad<FESpace1, FESpace2>(comp, fe1, fe2, offset_row, offset_clm);
}

template <typename FESpace1, typename FESpace2>
struct AssemblyDiv: public Coupling<FESpace1, FESpace2>
{
  using FESpace1_T = FESpace1;
  using FESpace2_T = FESpace2;
  using Super_T = Coupling<FESpace1, FESpace2>;
  using LMat_T = typename Super_T::LMat_T;
  using LVec_T = typename Super_T::LVec_T;

  explicit AssemblyDiv(uint comp,
                       FESpace1_T & fe1,
                       FESpace2_T & fe2,
                       uint offset_row = 0,
                       uint offset_clm = 0):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, offset_row, offset_clm),
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
    for(uint q=0; q<FESpace1_T::CurFE_T::QR_T::numPts; ++q)
    {
      Ke -= this->feSpace1.curFE.JxW[q] *
          this->feSpace1.curFE.phi[q] *
          this->feSpace2.curFE.dphi[q].col(component).transpose();
    }
  }

  uint const component;
};
