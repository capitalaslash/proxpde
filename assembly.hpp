#pragma once

#include "def.hpp"

template <typename CurFE>
struct Assembly
{
  typedef CurFE CurFE_T;
  typedef typename CurFE_T::LocalMat_T LMat_T;
  typedef typename CurFE_T::LocalVec_T LVec_T;

  explicit Assembly(CurFE_T & cfe):
    curFE(cfe)
  {}

  virtual void build(LMat_T & Ke, LVec_T & Fe) const = 0;

  CurFE_T & curFE;
};

template <typename CurFE>
struct AssemblyStiffness: public Assembly<CurFE>
{
  typedef CurFE CurFE_T;
  typedef Assembly<CurFE> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyStiffness(CurFE_T & cfe):
    Assembly<CurFE>(cfe)
  {}

  void build(LMat_T & Ke, LVec_T & Fe) const
  {
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->curFE.JxW[q] * (this->curFE.dphi(j, q).dot(this->curFE.dphi(i, q)));
        }
      }
    }
  }
};

template <typename CurFE>
struct AssemblyMass: public Assembly<CurFE>
{
  typedef CurFE CurFE_T;
  typedef Assembly<CurFE> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyMass(CurFE_T & cfe):
    Assembly<CurFE>(cfe)
  {}

  void build(LMat_T & Ke, LVec_T & Fe) const
  {
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->curFE.JxW[q] * this->curFE.phi(j, q) * this->curFE.phi(i, q);
        }
      }
    }
  }
};

template <typename CurFE>
struct AssemblyAnalyticalRhs: public Assembly<CurFE>
{
  typedef CurFE CurFE_T;
  typedef Assembly<CurFE> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyAnalyticalRhs(scalarFun_T const & r, CurFE_T & cfe):
    Assembly<CurFE>(cfe),
    rhs(r)
  {}

  void build(LMat_T &, LVec_T & Fe) const
  {
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double const f = this->rhs(this->curFE.qpoint[q]);
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        Fe(i) += this->curFE.JxW[q] * this->curFE.phi(i, q) * f;
      }
    }
  }

  scalarFun_T const & rhs;
};

template <typename CurFE>
struct AssemblyAdvection: public Assembly<CurFE>
{
  typedef CurFE CurFE_T;
  typedef Assembly<CurFE> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyAdvection(Vec3 const u, CurFE_T & cfe):
    Assembly<CurFE>(cfe),
    vel(u)
  {}

  void build(LMat_T & Ke, LVec_T & Fe) const
  {
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->curFE.JxW[q] * vel * this->curFE.dphi(i, q);
        }
      }
    }
  }

  Vec3 const vel;
};

template <typename CurFE>
struct AssemblyPoisson: public Assembly<CurFE>
{
  typedef CurFE CurFE_T;
  typedef Assembly<CurFE> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyPoisson(scalarFun_T const & r, CurFE_T & cfe):
    Assembly<CurFE>(cfe),
    assemblyLhs(cfe),
    assemblyRhs(r, cfe)
  {}

  void build(LMat_T & Ke, LVec_T & Fe) const
  {
    assemblyLhs.build(Ke, Fe);
    assemblyRhs.build(Ke, Fe);
  }

  AssemblyStiffness<CurFE> assemblyLhs;
  AssemblyAnalyticalRhs<CurFE> assemblyRhs;
};

template <typename FESpace>
void buildProblem(FESpace & feSpace,
                  Assembly<typename FESpace::CurFE_T> & assembly,
                  scalarFun_T const& rhs,
                  bc_list<FESpace> const& bcs,
                  Mat& A,
                  Vec& b)
{
  typedef typename FESpace::CurFE_T CurFE_T;
  typedef typename CurFE_T::LocalMat_T LMat_T;
  typedef typename CurFE_T::LocalVec_T LVec_T;

  auto & curFE = feSpace.curFE;

  std::vector<Tri> coefficients;
  // sparsity pattern
  coefficients.reserve(feSpace.meshPtr->pointList.size() * 3); // 3 = 2*dim+1

  for(auto &e: feSpace.meshPtr->elementList)
  {
    LMat_T Ke = LMat_T::Zero();
    LVec_T Fe = LVec_T::Zero();

    // --- set current fe ---
    curFE.reinit(e);

    // --- build local matrix and rhs ---
    assembly.build(Ke, Fe);

    // --- apply bc ---
    // A_constrained = C^T A C
    // b_constrained = C^T (b-Ah)
    // C clear constrained rows/cols
    // h is the vector of local constraint values

    for(auto& bc: bcs)
    {
      LMat_T C = LMat_T::Identity();
      LVec_T h = LVec_T::Zero();
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id = feSpace.dof._map[e.id][i];
        if(bc.is_constrained(id))
        {
          h(i) = bc.value(curFE.dofPts[i]);
          C(i,i) = 0.;
        }
      }
      Ke = C * Ke * C;
      Fe = C * (Fe - curFE.stiffMat * h);

      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id = feSpace.dof._map[e.id][i];
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
    for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
    {
      const DOFid_T id_i = feSpace.dof._map[e.id][i];
      b(id_i) += Fe(i);
      for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
      {
        const DOFid_T id_j = feSpace.dof._map[e.id][j];
        coefficients.push_back(Tri(id_i, id_j, Ke(i,j)));
      }
    }
  }
  A.setFromTriplets(coefficients.begin(), coefficients.end());
}
