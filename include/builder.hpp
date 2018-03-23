#pragma once

#include "def.hpp"
#include "blockmatrix.hpp"
#include "fespace.hpp"
#include "assembly.hpp"
#include "bc.hpp"

template<class T>
inline constexpr T pow(const T base, unsigned const exponent)
{
    return exponent == 0 ? 1 : base * pow(base, exponent-1);
}

struct Builder
{
  Builder(uint const size):
    A{size, size},
    b{Vec::Zero(size)}
  {}

  template <typename FESpace>
  void buildProblem(Diagonal<FESpace> const & assembly,
                    BCList<FESpace> & bcs)
  {
    using CurFE_T = typename FESpace::CurFE_T;
    using LMat_T = typename Diagonal<FESpace>::LMat_T;
    using LVec_T = typename Diagonal<FESpace>::LVec_T;
    auto const gSize = assembly.feSpace.dof.totalNum;

    auto const & mesh = *assembly.feSpace.meshPtr;

    // FIXME: compute a proper sparsity pattern
    // approxEntryNum = n. localMat entries * n. elements
    uint const approxEntryNum =
        pow(CurFE_T::numDOFs*FESpace::dim,2) * mesh.elementList.size();
    _triplets.reserve(approxEntryNum);

    for(auto const & e: mesh.elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      assembly.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke);

      // --- apply Neumann bcs ---
      for(auto & bc: bcs.bcNatList)
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
      // filelog << "Ke:\n" << Ke << std::endl;
      // filelog << "Fe:\n" << Fe << std::endl;


      // --- apply Mixed bcs ---
      for(auto & bc: bcs.bcMixedList)
      {
        uint facetCounter = 0;
        for(auto const facetId: mesh.elemToFacet[e.id])
        {
          // facet exists && has the marker set in the bc
          if(facetId != DOFidNotSet &&
             mesh.facetList[facetId].marker == bc.marker)
          {
            auto const & facet = mesh.facetList[facetId];
            bc.curFE.reinit(facet);
            for(uint q=0; q<BCMixed<FESpace>::QR_T::numPts; ++q)
            {
              for (uint d=0; d<assembly.feSpace.dim; ++d)
              {
                if (bc.hasComp(d))
                {
                  auto const coeff = bc.coeff(bc.curFE.qpoint[q])(d);
                  for(uint i=0; i<BCMixed<FESpace>::RefFE_T::numFuns; ++i)
                  {
                    auto const id =
                        CurFE_T::RefFE_T::dofOnFacet[facetCounter][i];
                    for (uint j=0; j<BCMixed<FESpace>::RefFE_T::numFuns; ++j)
                    {
                      auto const jd =
                          CurFE_T::RefFE_T::dofOnFacet[facetCounter][j];
                      Ke(id, jd) += coeff * bc.curFE.JxW[q] *
                                    bc.curFE.phi[q](j) *
                                    bc.curFE.phi[q](i);
                    }
                  }
                }
              }
            }
          }
          facetCounter++;
        }
      }
      // filelog << "Ke:\n" << Ke << std::endl;
      // filelog << "Fe:\n" << Fe << std::endl;

      // --- apply Dirichlet bcs ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      LMat_T C = LMat_T ::Identity();
      LVec_T h = LVec_T::Zero();
      for(auto& bc: bcs.bcEssList)
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

          for(auto& bc: bcs.bcEssList)
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

      // filelog << "\nelement" << e.id << "\n---------------" << std::endl;
      // filelog << "Ke:\n" << Ke << std::endl;
      // filelog << "Fe:\n" << Fe << std::endl;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = assembly.feSpace.dof.elemMap[e.id][i];
        for (uint d1=0; d1<FESpace::dim; ++d1)
        {
          b(assembly.offset_row + d1*gSize + id_i) += Fe(i + d1*FESpace::CurFE_T::size);

          for(uint j=0; j<CurFE_T::RefFE_T::numFuns; ++j)
          {
            DOFid_T const id_j = assembly.feSpace.dof.elemMap[e.id][j];
            for (uint d2=0; d2<FESpace::dim; ++d2)
            {
              auto val = Ke(i+d1*FESpace::CurFE_T::size,j+d2*FESpace::CurFE_T::size);
              if(std::fabs(val) > 1.e-16)
              {
                _triplets.emplace_back(
                      assembly.offset_row + d1*gSize + id_i,
                      assembly.offset_clm + d2*gSize + id_j,
                      val);
              }
            }
          }
        }
      }
    }
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
    using SquareMat1_T = FMat<FESpace1::dim*CurFE1_T::size,FESpace1::dim*CurFE1_T::size>;
    using SquareMat2_T = FMat<FESpace2::dim*CurFE2_T::size,FESpace2::dim*CurFE2_T::size>;
    using Vec2_T = FVec<FESpace2::dim*CurFE2_T::size>;
    auto const gSize1 = assembly.feSpace1.dof.totalNum;
    auto const gSize2 = assembly.feSpace2.dof.totalNum;

    // FIXME: compute a proper sparsity pattern
    // approxEntryNum = n. localMat entries * n. elements
    uint const approxEntryNum =
        CurFE1_T::numDOFs*FESpace1::dim *
        CurFE2_T::numDOFs*FESpace2::dim *
        assembly.feSpace1.meshPtr->elementList.size();
    _triplets.reserve(approxEntryNum);

    for(auto &e: assembly.feSpace1.meshPtr->elementList)
    {
      LMat_T Ke = LMat_T::Zero();
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      assembly.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Ke);

      // filelog << "\nelement" << e.id << "\n---------------" << std::endl;
      // filelog << "Ke:\n" << Ke << std::endl;

      // --- apply bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      SquareMat1_T Crow = SquareMat1_T::Identity();
      for(auto& bc: bcs1.bcEssList)
      {
        for (uint d=0; d<FESpace1::dim; ++d)
        {
          for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
          {
            auto const pos = i+d*FESpace1::RefFE_T::numFuns;
            DOFid_T const id = assembly.feSpace1.dof.elemMap[e.id][pos];
            if(bc.isConstrained(id))
            {
              Crow(pos,pos) = 0.;
            }
          }
        }
      }
      SquareMat2_T Cclm = SquareMat2_T::Identity();
      Vec2_T h = Vec2_T::Zero();
      for(auto& bc: bcs2.bcEssList)
      {
        for (uint d=0; d<FESpace2::dim; ++d)
        {
          for(uint i=0; i<CurFE2_T::RefFE_T::numFuns; ++i)
          {
            auto const pos = i+d*FESpace2::RefFE_T::numFuns;
            DOFid_T const id = assembly.feSpace2.dof.elemMap[e.id][pos];
            if(bc.isConstrained(id))
            {
              Cclm(pos,pos) = 0.;
              h(pos) = bc.value(assembly.feSpace2.curFE.dofPts[i])[d];
            }
          }
        }
      }
      Fe = - Crow * Ke * h;
      Ke = Crow * Ke * Cclm;

      // --- store local values in global matrix and rhs ---
      for(uint i=0; i<CurFE1_T::RefFE_T::numFuns; ++i)
      {
        DOFid_T const id_i = assembly.feSpace1.dof.elemMap[e.id][i];
        for (uint d1=0; d1<FESpace1::dim; ++d1)
        {
          b(assembly.offset_row + d1*gSize1 + id_i) += Fe(i + d1*FESpace1::CurFE_T::size);

          for(uint j=0; j<CurFE2_T::RefFE_T::numFuns; ++j)
          {
            DOFid_T const id_j = assembly.feSpace2.dof.elemMap[e.id][j];
            for (uint d2=0; d2<FESpace2::dim; ++d2)
            {
              auto val = Ke(i+d1*FESpace1::CurFE_T::size,j+d2*FESpace2::CurFE_T::size);
              if(std::fabs(val) > 1.e-16)
              {
                _triplets.emplace_back(
                      assembly.offset_row + d1*gSize1 + id_i,
                      assembly.offset_clm + d2*gSize2 + id_j,
                      val);
              }
            }
          }
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

    auto const & mesh = *assembly.feSpace.meshPtr;

    for(auto &e: mesh.elementList)
    {
      LVec_T Fe = LVec_T::Zero();

      // --- set current fe ---
      assembly.reinit(e);

      // --- build local matrix and rhs ---
      assembly.build(Fe);

      // --- apply bc ---
      // A_constrained = C^T A C
      // b_constrained = C^T (b-Ah)
      // C^T clears constrained rows
      // C clears constrained clms
      // h is the vector of local constraint values

      LMat_T C = LMat_T::Identity();
      for(auto& bc: bcs.bcEssList)
      {
        for (uint d=0; d<FESpace::dim; ++d)
        {
          for(uint i=0; i<CurFE_T::RefFE_T::numFuns; ++i)
          {
            DOFid_T const id = assembly.feSpace.dof.elemMap[e.id][i+d*FESpace::RefFE_T::numFuns];
            if(bc.isConstrained(id))
            {
              C(i+d*FESpace::RefFE_T::numFuns, i+d*FESpace::RefFE_T::numFuns) = 0.;
            }
          }
        }
      }
      Fe = C * Fe;

      // filelog << "\nelement" << e.id << "\n---------------" << std::endl;
      // filelog << "Fe:\n" << Fe << std::endl;

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

  void clear()
  {
    // there is no need to clear A as setFromTriplets() discards any content of the matrix
    b = Vec::Zero(b.size());
    _triplets.clear();
  }

  Mat A;
  Vec b;
  std::vector<Triplet> _triplets;
  array<std::vector<AssemblyBase*>,3> assemblies;
};
