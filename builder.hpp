#pragma once

#include "def.hpp"
#include "fespace.hpp"

struct Builder
{
  Builder(Mat & mat, Vec & vec):
    A(mat),
    b(vec)
  {}

  template <typename FESpace>
  void buildProblem(Diagonal<FESpace> const & assembly,
                    BCList<FESpace> const & bcs)
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
          _triplets.emplace_back(assembly.offset_row+id_i, assembly.offset_clm+id_j, Ke(i,j));
        }
      }
    }
  }

  template <typename FESpace1, typename FESpace2>
  void buildProblem(Coupling<FESpace1, FESpace2> const & assembly,
                    BCList<FESpace1> const & bcs1,
                    BCList<FESpace2> const & bcs2)
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

      // std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      // std::cout << "Ke:\n" << Ke << std::endl;

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
          _triplets.emplace_back(assembly.offset_row+id_i, assembly.offset_clm+id_j, Ke(i,j));
        }
      }
    }
  }

  template <typename FESpace>
  void buildProblem(AssemblyVector<FESpace> const & assembly,
                    BCList<FESpace> const& bcs)
  {
    typedef typename FESpace::CurFE_T CurFE_T;
    typedef typename AssemblyVector<FESpace>::LMat_T LMat_T;
    typedef typename AssemblyVector<FESpace>::LVec_T LVec_T;

    auto & curFE = assembly.feSpace.curFE;
    auto const & mesh = *assembly.feSpace.meshPtr;

    for(auto &e: mesh.elementList)
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

      for(auto & bc: bcs.bcNat_list)
      {
        uint facetCounter = 0;
        for(auto const facetId: mesh.elemToFacet[e.id])
        {
          if(facetId != DOFidNotSet &&
             mesh.facetList[facetId].marker == bc.marker)
          {
            auto const & facet = mesh.facetList[facetId];
            bc.curFE.reinit(facet);
            for(uint q=0; q<BCNat<FESpace>::CurFE_T::QR_T::numPts; ++q)
            {
              auto const value = bc.value(bc.curFE.qpoint[q]);
              typedef typename BCNat<FESpace>::CurFE_T::RefFE_T RefFESide_T;
              for(uint i=0; i<RefFESide_T::numFuns; ++i)
              {
                Fe(RefFESide_T::dofOnFacet[facetCounter][i]) += bc.curFE.JxW[q] *
                      bc.curFE.phi(i, q) *
                      value;
              }
            }
          }
          facetCounter++;
        }
      }

      std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      std::cout << "Fe:\n" << Fe << std::endl;

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
  std::vector<Triplet> _triplets;
};
