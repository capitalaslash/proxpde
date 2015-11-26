#pragma once

#include "def.hpp"
#include "geo.hpp"
#include "qr.hpp"
#include "curfe.hpp"
void buildProblem(std::shared_ptr<Mesh<Line>> const meshPtr,
                  scalarFun_T const& rhs,
                  bc_list const& bcs,
                  std::vector<Tri>& coefficients,
                  Vec& b)
{
  typedef CurFE<RefLineP1,GaussQR<3>> curFE_T;
  curFE_T curFE;
  for(auto &e: meshPtr->elementList)
  {
    // --- set current fe ---
    curFE.reinit(e);
    curFE_T::LocalMat_T Ke = curFE_T::LocalMat_T::Zero();
    curFE_T::LocalVec_T Fe = curFE_T::LocalVec_T::Zero();

    // --- build local matrix and rhs ---
    for(uint q=0; q<curFE_T::QR_T::numPts; ++q)
    {
      double const f = rhs(curFE.qpoint[q]);
      for(uint i=0; i<curFE_T::RefFE_T::numPts; ++i)
      {
        Fe(i) += curFE.JxW[q] * curFE.phi(i, q) * f;
        for(uint j=0; j<curFE_T::RefFE_T::numPts; ++j)
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
      curFE_T::LocalMat_T C = curFE_T::LocalMat_T::Identity();
      curFE_T::LocalVec_T h = curFE_T::LocalVec_T::Zero();
      for(uint i=0; i<curFE_T::RefFE_T::numPts; ++i)
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

      for(uint i=0; i<curFE_T::RefFE_T::numPts; ++i)
      {
        if(bc.is_constrained(*e.pointList[i]))
        {
          Ke(i,i) = curFE.Jm1;
          Fe(i) = h[i] * curFE.Jm1;
        }
      }
    }

    // --- store local values in global matrix and rhs ---
    for(uint i=0; i<curFE_T::RefFE_T::numPts; ++i)
    {
      const id_T id_i = e.pointList[i]->id;
      b(id_i) += Fe(i);
      for(uint j=0; j<curFE_T::RefFE_T::numPts; ++j)
      {
        const id_T id_j = e.pointList[j]->id;
        coefficients.push_back(Tri(id_i, id_j, Ke(i,j)));
      }

      // // neumann bc
      // if(id_i == meshPtr->pointList.size()-1)
      // {
      //   b(id_i) += 1.;
      // }
    }
  }
}
