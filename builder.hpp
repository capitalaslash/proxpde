#pragma once

#include "def.hpp"
#include "array.hpp"
#include "blockmatrix.hpp"
#include "fespace.hpp"
#include "assembly.hpp"
#include "bc.hpp"

struct Builder
{
  Builder(Mat & mat, Vec & vec):
    A(mat),
    b(vec)
  {}

  template <typename FESpace>
  void buildProblem(Diagonal<FESpace> const & assembly,
                    BCList<FESpace> & bcs)
  {
    using CurFE_T = typename FESpace::CurFE_T;
    using LMat_T = typename Diagonal<FESpace>::LMat_T;
    using LVec_T = typename Diagonal<FESpace>::LVec_T;
    auto const gSize = assembly.feSpace.dof.totalNum;

    // FIXME: compute a proper sparsity pattern
    uint const approxEntryNum = (2*CurFE_T::RefFE_T::dim+1) * gSize;
    _triplets.reserve(approxEntryNum);

    auto & curFE = assembly.feSpace.curFE;
    auto const & mesh = *assembly.feSpace.meshPtr;

    uint entryNum = 0;
    for(auto &e: mesh.elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      curFE.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke);

      // --- apply Neumann bcs ---
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
            for(uint q=0; q<BCNat<FESpace>::QR_T::numPts; ++q)
            {
              for (uint d=0; d<assembly.feSpace.dim; ++d)
              {
                if (bc.hasComp(d))
                {
                  auto const value = bc.value(bc.curFE.qpoint[q])(d);
                  for(uint i=0; i<BCNat<FESpace>::RefFE_T::numFuns; ++i)
                  {
                    auto const id =
                        CurFE_T::RefFE_T::dofOnFacet[facetCounter][i];
                    Fe(id) += bc.curFE.JxW[q] *
                              bc.curFE.phi[q](i) *
                        value;
                  }
                }
              }
            }
          }
          facetCounter++;
        }
      }
      // std::cout << "Fe:\n" << Fe << std::endl;

      // --- apply Dirichlet bcs ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      LMat_T C = LMat_T ::Identity();
      LVec_T h = LVec_T::Zero();
      for(auto& bc: bcs.bcEss_list)
      {
        for (uint d=0; d<FESpace::dim; ++d)
        {
          for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
          {
            auto const pos = i+d*FESpace::RefFE_T::numFuns;
            DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][pos];
            if(bc.isConstrained(id))
            {
              C(pos, pos) = 0.;
              h(pos) = bc.value(assembly.feSpace.curFE.dofPts[i])(d);
            }
          }
        }
      }
      Fe = C * (Fe - Ke * h);
      Ke = C * Ke * C;

      for (uint d=0; d<FESpace::dim; ++d)
      {
        for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
        {
          auto const pos = i+d*FESpace::RefFE_T::numFuns;
          DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][pos];

          for(auto& bc: bcs.bcEss_list)
          {
            // dofs can be fixed with essential conditions only once!
            // all other bcs will be silently discarded
            if (bc.isConstrained(id))
            {
              Ke(pos, pos) = bc.diag;
              Fe(pos) = h[pos];
              // bc.fixedDofs.insert(id);
            }
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
        for (uint d=0; d<FESpace::dim; ++d)
        {
          b(assembly.offset_row + d*gSize + id_i) += Fe(i + d*FESpace::CurFE_T::size);
          for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
          {
            DOFid_T const id_j = assembly.feSpace.dof.elemMap[e.id][j];
            auto val = Ke(i+d*FESpace::CurFE_T::size,j+d*FESpace::CurFE_T::size);
            if(std::fabs(val) > 1.e-16)
            {
              _triplets.emplace_back(
                    assembly.offset_row + d*gSize + id_i,
                    assembly.offset_clm + d*gSize + id_j,
                    val);
              entryNum++;
            }
          }
        }
      }
    }
    std::cout << "entries: " << entryNum << " - expected " << approxEntryNum << std::endl;
  }

  template <typename FESpace1, typename FESpace2>
  void buildProblem(Coupling<FESpace1, FESpace2> const & assembly,
                    BCList<FESpace1> const & bcs1,
                    BCList<FESpace2> const & bcs2)
  {
    using CurFE1_T = typename FESpace1::CurFE_T;
    using CurFE2_T = typename FESpace2::CurFE_T;
    using LMat_T = typename Coupling<FESpace1, FESpace2>::LMat_T;
    using LVec_T = typename Coupling<FESpace1, FESpace2>::LVec_T;

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
      for(auto& bc: bcs1.bcEss_list)
      {
        for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace1.dof.elemMap[e.id][i];
          if(bc.isConstrained(id))
          {
            Crow(i,i) = 0.;
          }
        }
      }
      typename CurFE2_T::LocalMat_T Cclm = CurFE2_T::LocalMat_T::Identity();
      typename CurFE2_T::LocalVec_T h = CurFE2_T::LocalVec_T::Zero();
      for(auto& bc: bcs2.bcEss_list)
      {
        for(uint i=0; i<CurFE2_T::RefFE_T::numFuns; ++i)
        {
          DOFid_T const id = assembly.feSpace2.dof.elemMap[e.id][i];
          if(bc.isConstrained(id))
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
                    BCList<FESpace> const & bcs)
  {
    using CurFE_T = typename FESpace::CurFE_T;
    using LMat_T = typename AssemblyVector<FESpace>::LMat_T;
    using LVec_T = typename AssemblyVector<FESpace>::LVec_T;

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
      for(auto& bc: bcs.bcEss_list)
      {
        for (uint d=0; d<FESpace::dim; ++d)
        {
          for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
          {
            DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i+d*FESpace::RefFE_T::numFuns];
            if(bc.isConstrained(id))
            {
              C(i+d*FESpace::RefFE_T::numFuns,i+d*FESpace::RefFE_T::numFuns) = 0.;
            }
          }
        }
      }
      Fe = C * Fe;

      // std::cout << "\nelement" << e.id << "\n---------------" << std::endl;
      // std::cout << "Fe:\n" << Fe << std::endl;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = assembly.feSpace.dof.elemMap[e.id][i];
        for (uint d=0; d<FESpace::dim; ++d)
        {
          auto val = Fe(i+d*FESpace::CurFE_T::size);
          b(assembly.offset_row + d*assembly.feSpace.dof.totalNum + id_i) += val;
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
  std::vector<Triplet> _triplets;
  array<std::vector<AssemblyBase*>,3> assemblies;
};
