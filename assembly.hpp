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
  typedef Eigen::Matrix<double,CurFE_T::size(),CurFE_T::size()> LMat_T;
  typedef Eigen::Matrix<double,CurFE_T::size(),1> LVec_T;

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
  typedef Eigen::Matrix<double,CurFE1_T::size(),CurFE2_T::size()> LMat_T;
  typedef Eigen::Matrix<double,CurFE1_T::size(),1> LVec_T;

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
  typedef Eigen::Matrix<double,CurFE_T::size(),CurFE_T::size()> LMat_T;
  typedef Eigen::Matrix<double,CurFE_T::size(),1> LVec_T;

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
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace.curFE.JxW[q] *
              (this->feSpace.curFE.dphi(j, q).dot(this->feSpace.curFE.dphi(i, q)));
        }
      }
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

  explicit AssemblyMass(FESpace_T & fe, uint offset_row = 0, uint offset_clm = 0):
     Diagonal<FESpace>(fe, offset_row, offset_clm)
  {}

  void build(LMat_T & Ke) const
  {
    typedef typename FESpace_T::CurFE_T CurFE_T;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace.curFE.JxW[q] *
              this->feSpace.curFE.phi(j, q) * this->feSpace.curFE.phi(i, q);
        }
      }
    }
  }
};

template <typename FESpace>
struct AssemblyAnalyticalRhs: public AssemblyVector<FESpace>
{
  typedef FESpace FESpace_T;
  typedef Diagonal<FESpace> Super_T;
  typedef typename Super_T::LMat_T LMat_T;
  typedef typename Super_T::LVec_T LVec_T;

  explicit AssemblyAnalyticalRhs(scalarFun_T const & r,
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
        Fe(i) += this->feSpace.curFE.JxW[q] * this->feSpace.curFE.phi(i, q) * f;
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
            this->feSpace.curFE.phi(i, q) *
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
        local_vel += vel(dofId) * this->feSpace.curFE.phi(n, q);
      }
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          Ke(i,j) += this->feSpace.curFE.JxW[q] *
              local_vel *
              this->feSpace.curFE.dphi(i, q) *
              this->feSpace.curFE.phi(j, q);
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
          Ke(i,j) += this->feSpace1.curFE.JxW[q] *
              this->feSpace1.curFE.phi(i, q) *
              this->feSpace2.curFE.dphi(j, q)(component);
        }
      }
    }
  }

  uint const component;
};

struct Builder
{
  Builder(Mat & mat, Vec & vec):
    A(mat),
    b(vec)
  {}

  template <typename FESpace>
  void buildProblem(Diagonal<FESpace> const & assembly,
                    bc_list<FESpace> const & bcs,
                    uint row_offset = 0,
                    uint clm_offset = 0)
  {
    typedef typename FESpace::CurFE_T CurFE_T;
    typedef typename Diagonal<FESpace>::LMat_T LMat_T;
    typedef typename Diagonal<FESpace>::LVec_T LVec_T;

    // FIXME - compute a proper sparsity pattern
    _triplets.reserve((2*CurFE_T::RefFE_T::dim+1) * assembly.feSpace.dof.totalNum);

    auto & curFE = assembly.feSpace.curFE;

    for(auto &e: assembly.feSpace.meshPtr->elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      curFE.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke);

      // --- apply bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      typename CurFE_T::LocalMat_T C = CurFE_T::LocalMat_T::Identity();
      typename CurFE_T::LocalVec_T h = CurFE_T::LocalVec_T::Zero();
      for(auto& bc: bcs)
      {
        for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i];
          if(bc.is_constrained(id))
          {
            C(i,i) = 0.;
          }
        }
        for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i];
          if(bc.is_constrained(id))
          {
            h(i) = bc.value(assembly.feSpace.curFE.dofPts[i]);
          }
        }
      }
      Fe = C * (Fe - Ke * h);
      Ke = C * Ke * C;

      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i];
        for(auto& bc: bcs)
        {
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
        DOFid_T const id_i = assembly.feSpace.dof.elemMap[e.id][i];
        b(assembly.offset_row+id_i) += Fe(i);
        for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
        {
          DOFid_T const id_j = assembly.feSpace.dof.elemMap[e.id][j];
          _triplets.push_back(Tri(assembly.offset_row+id_i, assembly.offset_clm+id_j, Ke(i,j)));
        }
      }
    }
  }

  template <typename FESpace1, typename FESpace2>
  void buildProblem(Coupling<FESpace1, FESpace2> const & assembly,
                    bc_list<FESpace1> const & bcs1,
                    bc_list<FESpace2> const & bcs2,
                    uint row_offset = 0,
                    uint clm_offset = 0)
  {
    typedef typename FESpace1::CurFE_T CurFE1_T;
    typedef typename FESpace2::CurFE_T CurFE2_T;
    typedef typename Coupling<FESpace1, FESpace2>::LMat_T LMat_T;
    typedef typename Coupling<FESpace1, FESpace2>::LVec_T LVec_T;

    // FIXME - compute a proper sparsity pattern
    _triplets.reserve((2*CurFE1_T::RefFE_T::dim+1) * assembly.feSpace1.dof.totalNum);

    auto & curFE1 = assembly.feSpace1.curFE;
    auto & curFE2 = assembly.feSpace2.curFE;

    for(auto &e: assembly.feSpace1.meshPtr->elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      curFE1.reinit(e);
      curFE2.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke);

      std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      std::cout << "Ke:\n" << Ke << std::endl;

      // --- apply bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      typename CurFE1_T::LocalMat_T Crow = CurFE1_T::LocalMat_T::Identity();
      for(auto& bc: bcs1)
      {
        for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace1.dof.elemMap[e.id][i];
          if(bc.is_constrained(id))
          {
            Crow(i,i) = 0.;
          }
        }
      }
      typename CurFE2_T::LocalMat_T Cclm = CurFE2_T::LocalMat_T::Identity();
      typename CurFE2_T::LocalVec_T h = CurFE2_T::LocalVec_T::Zero();
      for(auto& bc: bcs2)
      {
        for(uint i=0; i<CurFE2_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace2.dof.elemMap[e.id][i];
          if(bc.is_constrained(id))
          {
            Cclm(i,i) = 0.;
            h(i) = bc.value(assembly.feSpace2.curFE.dofPts[i]);
          }
        }
      }
      Fe = - Crow * Ke * h;
      Ke = Crow * Ke * Cclm;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = assembly.feSpace1.dof.elemMap[e.id][i];
        b(assembly.offset_row+id_i) += Fe(i);
        for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
        {
          DOFid_T const id_j = assembly.feSpace2.dof.elemMap[e.id][j];
          _triplets.push_back(Tri(assembly.offset_row+id_i, assembly.offset_clm+id_j, Ke(i,j)));
        }
      }
    }
  }

  template <typename FESpace>
  void buildProblem(AssemblyVector<FESpace> const & assembly,
                    bc_list<FESpace> const& bcs,
                    uint row_offset = 0,
                    uint clm_offset = 0)
  {
    typedef typename FESpace::CurFE_T CurFE_T;
    typedef typename AssemblyVector<FESpace>::LMat_T LMat_T;
    typedef typename AssemblyVector<FESpace>::LVec_T LVec_T;

    auto & curFE = assembly.feSpace.curFE;

    for(auto &e: assembly.feSpace.meshPtr->elementList)
    {
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      curFE.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Fe);

      // --- apply bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      LMat_T C = LMat_T::Identity();
      for(auto& bc: bcs)
      {
        for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i];
          if(bc.is_constrained(id))
          {
            C(i,i) = 0.;
          }
        }
      }
      Fe = C * Fe;

      std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      std::cout << "Fe:\n" << Fe << std::endl;

      assert(assembly.offset_row == row_offset);
      assert(assembly.offset_clm == clm_offset);
      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = assembly.feSpace.dof.elemMap[e.id][i];
        b(assembly.offset_row+id_i) += Fe(i);
      }
    }
  }

  void closeMatrix()
  {
    A.setFromTriplets(_triplets.begin(), _triplets.end());
  }

  Mat & A;
  Vec & b;
  std::array<std::vector<AssemblyBase*>,3> assemblies;
  std::vector<Tri> _triplets;
};
