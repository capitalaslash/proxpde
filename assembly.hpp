#pragma once

#include "def.hpp"
#include "fespace.hpp"

template <typename FESpace1, typename FESpace2=FESpace1>
struct Assembly
{
  typedef FESpace1 FESpace1_T;
  typedef FESpace2 FESpace2_T;
  typedef typename FESpace1_T::CurFE_T CurFE1_T;
  typedef typename FESpace2_T::CurFE_T CurFE2_T;
  typedef Eigen::Matrix<double,CurFE1_T::size(),CurFE2_T::size()> LMat_T;
  typedef Eigen::Matrix<double,CurFE1_T::size(),1> LVec_T;

  explicit Assembly(FESpace1_T & fe1):
    feSpace1(fe1),
    feSpace2(fe1)
  {}

  Assembly(FESpace1_T & fe1, FESpace2_T & fe2):
    feSpace1(fe1),
    feSpace2(fe2)
  {}

  virtual void build(LMat_T & Ke, LVec_T & Fe) const = 0;

  FESpace1_T & feSpace1;
  FESpace2_T & feSpace2;
};

template <typename FESpace>
struct AssemblyStiffness: public Assembly<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Assembly<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyStiffness(FESpace_T & fe):
    Assembly<FESpace_T>(fe)
  {}

  void build(LMat_T & Ke, LVec_T & Fe) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace1.curFE.JxW[q] *
              (this->feSpace1.curFE.dphi(j, q).dot(this->feSpace2.curFE.dphi(i, q)));
        }
      }
    }
  }
};

template <typename FESpace>
struct AssemblyMass: public Assembly<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Assembly<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyMass(FESpace_T & fe):
     Assembly<FESpace>(fe)
  {}

  void build(LMat_T & Ke, LVec_T & Fe) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace1.curFE.JxW[q] *
              this->feSpace1.curFE.phi(j, q) * this->feSpace1.curFE.phi(i, q);
        }
      }
    }
  }
};

template <typename FESpace>
struct AssemblyAnalyticalRhs: public Assembly<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Assembly<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyAnalyticalRhs(scalarFun_T const & r, FESpace & fe):
    Assembly<FESpace>(fe),
    rhs(r)
  {}

  void build(LMat_T &, LVec_T & Fe) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double const f = rhs(this->feSpace1.curFE.qpoint[q]);
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        Fe(i) += this->feSpace1.curFE.JxW[q] * this->feSpace1.curFE.phi(i, q) * f;
      }
    }
  }

  scalarFun_T const & rhs;
};

template <typename FESpace>
struct AssemblyVecRhs: public Assembly<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Assembly<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyVecRhs(Vec const & r, FESpace_T & fe):
    Assembly<FESpace_T>(fe),
    rhs(r)
  {}

  void build(LMat_T &, LVec_T & Fe) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double local_rhs = 0.0;
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = this->feSpace1.dof.elemMap[this->feSpace1.curFE.e->id][n];
        local_rhs += rhs(dofId) * this->feSpace1.curFE.phi(n, q);
      }
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        Fe(i) += this->feSpace1.curFE.JxW[q] *
            this->feSpace1.curFE.phi(i, q) *
            local_rhs;
      }
    }
  }

  Vec const & rhs;
};

template <typename FESpace>
struct AssemblyAdvection: public Assembly<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Assembly<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyAdvection(Vec const u,
                             FESpace_T & fe):
    Assembly<FESpace>(fe),
    vel(u)
  {}

  void build(LMat_T & Ke, LVec_T & Fe) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      Vec3 local_vel = Vec3::Zero();
      for(uint n=0; n<CurFE_T::RefFE_T::numFuns; ++n)
      {
        id_T const dofId = this->feSpace1.dof.elemMap[this->feSpace1.curFE.e->id][n];
        local_vel += vel(dofId) * this->feSpace1.curFE.phi(n, q);
      }
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace1.curFE.JxW[q] *
              local_vel *
              this->feSpace1.curFE.dphi(i, q) *
              this->feSpace2.curFE.phi(j, q);
        }
      }
    }
  }

  Vec const & vel;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyGrad: public Assembly<FESpace1, FESpace2>
{
  typedef FESpace1 FESpace1_T;
  typedef FESpace2 FESpace2_T;
  typedef Assembly<FESpace1, FESpace2> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyGrad(uint comp, FESpace1_T & fe1, FESpace2_T & fe2):
    Assembly<FESpace1_T,FESpace2_T>(fe1, fe2),
    component(comp)
  {}

  void build(LMat_T & Ke, LVec_T &) const
  {
    typedef typename FESpace1_T::CurFE_T CurFE1_T;
    typedef typename FESpace2_T::CurFE_T CurFE2_T;
    for(uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace1.curFE.JxW[q] *
              this->feSpace1.curFE.dphi(i, q)(component) *
              this->feSpace2.curFE.phi(j, q);
        }
      }
    }
  }
  uint const component;
};

template <typename FESpace1, typename FESpace2>
struct AssemblyDiv: public Assembly<FESpace1, FESpace2>
{
  typedef FESpace1 FESpace1_T;
  typedef FESpace2 FESpace2_T;
  typedef Assembly<FESpace1, FESpace2> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyDiv(uint comp, FESpace1_T & fe1, FESpace2_T & fe2):
    Assembly<FESpace1_T,FESpace2_T>(fe1, fe2),
    component(comp)
  {}

  void build(LMat_T & Ke, LVec_T &) const
  {
    typedef typename FESpace1_T::CurFE_T CurFE1_T;
    typedef typename FESpace2_T::CurFE_T CurFE2_T;
    for(uint q=0; q<CurFE1_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace1.curFE.JxW[q] *
              this->feSpace1.curFE.phi(i, q) *
              this->feSpace2.curFE.dphi(j, q)(component);
        }
      }
    }
  }

  uint const component;
};

template <typename FESpace>
struct AssemblyPoisson: public Assembly<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Assembly<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyPoisson(scalarFun_T const & r,
                           FESpace_T & fe):
    Assembly<FESpace>(fe),
    assemblyLhs(fe),
    assemblyRhs(r, fe)
  {}

  void build(LMat_T & Ke, LVec_T & Fe) const
  {
    assemblyLhs.build(Ke, Fe);
    assemblyRhs.build(Ke, Fe);
  }

  AssemblyStiffness<FESpace> assemblyLhs;
  AssemblyAnalyticalRhs<FESpace> assemblyRhs;
};

struct Builder
{
  Builder(Mat & mat, Vec & vec):
    A(mat),
    b(vec)
  {}

  template <typename FESpace>
  void buildProblem(FESpace & feSpace,
                    Assembly<FESpace> & assembly,
                    bc_list<FESpace> const& bcs,
                    uint row_offset = 0,
                    uint clm_offset = 0)
  {
    buildProblem(feSpace, feSpace, assembly, bcs, row_offset, clm_offset);
  }

  template <typename FESpace1, typename FESpace2>
  void buildProblem(FESpace1 & feSpace1,
                    FESpace2 & feSpace2,
                    Assembly<FESpace1, FESpace2> & assembly,
                    bc_list<FESpace1> const& bcs,
                    uint row_offset = 0,
                    uint clm_offset = 0)
  {
    typedef typename FESpace1::CurFE_T CurFE1_T;
    typedef typename FESpace2::CurFE_T CurFE2_T;
    typedef typename Assembly<FESpace1, FESpace2>::LMat_T LMat_T;
    typedef typename Assembly<FESpace1, FESpace2>::LVec_T LVec_T;

    // FIXME - compute a proper sparsity pattern
    _triplets.reserve((2*CurFE1_T::RefFE_T::dim+1) * feSpace1.dof.totalNum);

    auto & curFE = feSpace1.curFE;

    for(auto &e: feSpace1.meshPtr->elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      curFE.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke, Fe);

      // --- apply bc ---
      // A_constrained = Cl A Cr
      // b_constrained = C^T (b-Ah)
      // Cl clears constrained rows
      // Cr clears constrained clms
      // h is the vector of local constraint values

      auto Ke_uncostrained = Ke;
      for(auto& bc: bcs)
      {
        typename CurFE1_T::LocalMat_T Cl = CurFE1_T::LocalMat_T::Identity();
        typename CurFE2_T::LocalMat_T Cr = CurFE2_T::LocalMat_T::Identity();
        typename CurFE2_T::LocalVec_T h = CurFE2_T::LocalVec_T::Zero();
        for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = feSpace1.dof.elemMap[e.id][i];
          if(bc.is_constrained(id))
          {
            Cl(i,i) = 0.;
          }
        }
        for(uint i=0; i<CurFE2_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = feSpace2.dof.elemMap[e.id][i];
          if(bc.is_constrained(id))
          {
            h(i) = bc.value(feSpace2.curFE.dofPts[i]);
            Cr(i,i) = 0.;
          }
        }
        Ke = Cl * Ke * Cr;
        Fe = Cl * (Fe - Ke_uncostrained * h);

        for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = feSpace1.dof.elemMap[e.id][i];
          if(bc.is_constrained(id))
          {
            Ke(i,i) = 1.0;
            Fe(i) = h[i];
          }
        }
      }

      // std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      // std::cout << "Ke:\n" << Ke << std::endl;
      // std::cout << "Fe:\n" << Fe << std::endl;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = feSpace1.dof.elemMap[e.id][i];
        b(row_offset+id_i) += Fe(i);
        for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
        {
          DOFid_T const id_j = feSpace2.dof.elemMap[e.id][j];
          _triplets.push_back(Tri(row_offset+id_i, clm_offset+id_j, Ke(i,j)));
        }
      }
    }
  }

  void closeMatrix()
  {
    A.setFromTriplets(_triplets.begin(), _triplets.end());
  }

  Mat & A;
  Vec & b;
  std::vector<Tri> _triplets;
};
