#pragma once

#include "def.hpp"
#include "fespace.hpp"

struct AssemblyBase
{
  uint offset_row;
  uint offset_clm;
};

template <typename FESpace>
struct Diagonal: public AssemblyBase
{
  typedef FESpace FESpace_T;
  typedef typename FESpace_T::CurFE_T CurFE_T;
  typedef FMat<CurFE_T::size(),CurFE_T::size()> LMat_T;
  typedef FVec<CurFE_T::size()> LVec_T;

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
  typedef FESpace1 FESpace1_T;
  typedef FESpace2 FESpace2_T;
  typedef typename FESpace1_T::CurFE_T CurFE1_T;
  typedef typename FESpace2_T::CurFE_T CurFE2_T;
  typedef FMat<CurFE1_T::size(),CurFE2_T::size()> LMat_T;
  typedef FVec<CurFE1_T::size()> LVec_T;

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
  typedef FESpace FESpace_T;
  typedef typename FESpace_T::CurFE_T CurFE_T;
  typedef FMat<CurFE_T::size(),CurFE_T::size()> LMat_T;
  typedef FVec<CurFE_T::size()> LVec_T;

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
  typedef FESpace FESpace_T;
  typedef Diagonal<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyStiffness(FESpace_T & fe, uint offset_row = 0, uint offset_clm = 0):
    Diagonal<FESpace_T>(fe, offset_row, offset_clm)
  {}

  void build(LMat_T & Ke) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Ke +=
        this->feSpace.curFE.JxW[q] *
        this->feSpace.curFE.dphi[q] *
        this->feSpace.curFE.dphi[q].transpose();
    }
  }
};

template <typename FESpace>
struct AssemblyMass: public Diagonal<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Diagonal<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyMass(double const & c, FESpace_T & fe, uint offset_row = 0, uint offset_clm = 0):
     Diagonal<FESpace>(fe, offset_row, offset_clm),
     coeff(c)
  {}

  void build(LMat_T & Ke) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Ke +=
        coeff * this->feSpace.curFE.JxW[q] *
        this->feSpace.curFE.phi[q] *
        this->feSpace.curFE.phi[q].transpose();
    }
  }

  double coeff;
};

template <typename FESpace>
struct AssemblyAnalyticRhs: public AssemblyVector<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Diagonal<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyAnalyticRhs(scalarFun_T const & r,
                                 FESpace & fe,
                                 uint offset_row = 0):
    AssemblyVector<FESpace>(fe, offset_row),
    rhs(r)
  {}

  void build(LVec_T & Fe) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double const f = rhs(this->feSpace.curFE.qpoint[q]);
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        Fe(i) += this->feSpace.curFE.JxW[q] * this->feSpace.curFE.phi[q](i) * f;
      }
    }
  }

  scalarFun_T const & rhs;
};

template <typename FESpace>
struct AssemblyVecRhs: public AssemblyVector<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Diagonal<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyVecRhs(Vec const & r, FESpace_T & fe):
    Diagonal<FESpace_T>(fe),
    rhs(r)
  {}

  void build(LVec_T & Fe) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double local_rhs = 0.0;
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = this->feSpace.dof.elemMap[this->feSpace.curFE.e->id][n];
        local_rhs += rhs(dofId) * this->feSpace.curFE.phi(n, q);
      }
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        Fe(i) += this->feSpace.curFE.JxW[q] *
            this->feSpace.curFE.phi[q](i) *
            local_rhs;
      }
    }
  }

  Vec const & rhs;
};

template <typename FESpace>
struct AssemblyAdvection: public Diagonal<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Diagonal<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyAdvection(Vec const u,
                             FESpace_T & fe):
    Diagonal<FESpace>(fe),
    vel(u)
  {}

  void build(LMat_T & Ke) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Vec3 local_vel = Vec3::Zero();
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = this->feSpace.dof.elemMap[this->feSpace.curFE.e->id][n];
        local_vel += vel(dofId) * this->feSpace.curFE.phi[q](n);
      }
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace.curFE.JxW[q] *
              local_vel *
              this->feSpace.curFE.dphi[q].row(i) *
              this->feSpace.curFE.phi[q](j);
        }
      }
    }
  }

  Vec const & vel;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyGrad: public Coupling<FESpace1, FESpace2>
{
  typedef FESpace1 FESpace1_T;
  typedef FESpace2 FESpace2_T;
  typedef Coupling<FESpace1, FESpace2> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyGrad(uint comp,
                        FESpace1_T & fe1,
                        FESpace2_T & fe2,
                        uint offset_row = 0,
                        uint offset_clm = 0):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, offset_row, offset_clm),
    component(comp)
  {}

  void build(LMat_T & Ke) const
  {
    typedef typename FESpace1_T::CurFE_T CurFE1_T;
    typedef typename FESpace2_T::CurFE_T CurFE2_T;
    for(uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) -= this->feSpace1.curFE.JxW[q] *
              this->feSpace1.curFE.dphi[q](i, component) *
              this->feSpace2.curFE.phi[q](j);
        }
      }
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
  typedef FESpace1 FESpace1_T;
  typedef FESpace2 FESpace2_T;
  typedef Coupling<FESpace1, FESpace2> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyDiv(uint comp,
                       FESpace1_T & fe1,
                       FESpace2_T & fe2,
                       uint offset_row = 0,
                       uint offset_clm = 0):
    Coupling<FESpace1_T,FESpace2_T>(fe1, fe2, offset_row, offset_clm),
    component(comp)
  {}

  void build(LMat_T & Ke) const
  {
    typedef typename FESpace1_T::CurFE_T CurFE1_T;
    typedef typename FESpace2_T::CurFE_T CurFE2_T;
    for(uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) -= this->feSpace1.curFE.JxW[q] *
              this->feSpace1.curFE.phi[q](i) *
              this->feSpace2.curFE.dphi[q](j, component);
        }
      }
    }
  }

  uint const component;
};
