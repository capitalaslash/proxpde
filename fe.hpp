#pragma once

#include "def.hpp"
#include "qr.hpp"
#include "curfe.hpp"

template <typename Elem, uint order>
struct FEType {};

template <>
struct FEType<Line,1>
{
  typedef RefLineP1 RefFE_T;
};

template <>
struct FEType<Line,2>
{
  typedef RefLineP2 RefFE_T;
};

template <typename Mesh,
          typename RefFE,
          typename QR>
struct FESpace
{
  typedef Mesh Mesh_T;
  typedef RefFE RefFE_T;
  typedef QR QR_T;
  typedef CurFE<RefFE,QR> CurFE_T;

  explicit FESpace(std::shared_ptr<Mesh> const mesh):
    meshPtr(mesh)
  {}

  std::shared_ptr<Mesh> const meshPtr;
  CurFE_T curFE;
};

template <typename FESpace>
void buildProblem(FESpace feSpace,
                  scalarFun_T const& rhs,
                  bc_list<typename FESpace::Mesh_T> const& bcs,
                  Mat& A,
                  Vec& b)
{
  typedef typename FESpace::CurFE_T CurFE_T;
  typedef typename CurFE_T::LocalMat_T LMat_T;
  typedef typename CurFE_T::LocalVec_T LVec_T;

  CurFE_T & curFE = feSpace.curFE;

  std::vector<Tri> coefficients;
  // sparsity pattern
  coefficients.reserve(feSpace.meshPtr->pointList.size() * 3); // 3 = 2*dim+1

  for(auto &e: feSpace.meshPtr->elementList)
  {
    // --- set current fe ---
    curFE.reinit(e);
    LMat_T Ke = LMat_T::Zero();
    LVec_T Fe = LVec_T::Zero();

    // --- build local matrix and rhs ---
    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      double const f = rhs(curFE.qpoint[q]);
      for(uint i=0; i<CurFE_T::RefFE_T::numPts; ++i)
      {
        Fe(i) += curFE.JxW[q] * curFE.phi(i, q) * f;
        for(uint j=0; j<CurFE_T::RefFE_T::numPts; ++j)
        {
          Ke(i,j) += curFE.JxW[q] * (curFE.dphi(j, q).dot(curFE.dphi(i, q)));
        }
      }
    }

    // --- apply bc ---
    // A_constrained = C^T A C
    // b_constrained = C^T (b-Ah)
    // C clear constrained rows/cols
    // h is the vector of local constraint values

    for(auto& bc: bcs)
    {
      LMat_T C = LMat_T::Identity();
      LVec_T h = LVec_T::Zero();
      for(uint i=0; i<CurFE_T::RefFE_T::numPts; ++i)
      {
        Point const& p = *e.pointList[i];
        if(bc.is_constrained(p))
        {
          h(i) = bc.value(p.coord);
          C(i,i) = 0.;
        }
      }
      Ke = C * Ke * C;
      Fe = C * (Fe - curFE.stiffMat * h);

      for(uint i=0; i<CurFE_T::RefFE_T::numPts; ++i)
      {
        if(bc.is_constrained(*e.pointList[i]))
        {
          Ke(i,i) = curFE.Jm1;
          Fe(i) = h[i] * curFE.Jm1;
        }
      }
    }

    // --- store local values in global matrix and rhs ---
    for(uint i=0; i<CurFE_T::RefFE_T::numPts; ++i)
    {
      const id_T id_i = e.pointList[i]->id;
      b(id_i) += Fe(i);
      for(uint j=0; j<CurFE_T::RefFE_T::numPts; ++j)
      {
        const id_T id_j = e.pointList[j]->id;
        coefficients.push_back(Tri(id_i, id_j, Ke(i,j)));
      }
    }
  }
  A.setFromTriplets(coefficients.begin(), coefficients.end());
}
