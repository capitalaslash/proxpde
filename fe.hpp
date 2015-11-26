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
  typedef CurFE<RefLineP1,GaussQR<3>> FE_T;
  FE_T curFE;
  for(auto &e: meshPtr->elementList)
  {
    // --- set current fe ---
    curFE.reinit(e);

    // --- build local matrix and rhs ---
    FE_T::LocalMat_T elemMat_c = curFE.stiffMat;
    FE_T::LocalVec_T elemRhs_c = FE_T::LocalVec_T::Zero();

    for(uint q=0; q<FE_T::QR_T::numPts; ++q)
    {
      double const f = rhs(curFE.qpoint[q]);
      for(uint i=0; i<FE_T::RefFE_T::numPts; ++i)
      {
        const id_T id_i = e.pointList[i]->id;
        b(id_i) += (1-bcs.vec[id_i])*curFE.JxW[q]*curFE.phi(i, q) * f;
      }
    }

    // --- apply bc ---
    // A_constrained = C^T A C
    // b_constrained = C^T (b-Ah)
    // C clear constrained rows/cols
    // h is the vector of local constraint values

    for(auto& bc: bcs)
    {
      FE_T::LocalMat_T C = FE_T::LocalMat_T::Identity();
      FE_T::LocalVec_T h = FE_T::LocalVec_T::Zero();
      for(uint i=0; i<FE_T::RefFE_T::numPts; ++i)
      {
        Point const& p = *e.pointList[i];
        if(bc.is_constrained(p))
        {
          h(i) = bc.value(p.coord);
          C(i,i) = 0.;
        }
      }
      elemMat_c = C * elemMat_c * C;
      elemRhs_c = C * (elemRhs_c - curFE.stiffMat * h);

      for(uint i=0; i<FE_T::RefFE_T::numPts; ++i)
      {
        if(bc.is_constrained(*e.pointList[i]))
        {
          elemMat_c(i,i) = curFE.Jm1;
          elemRhs_c(i) = h[i] * curFE.Jm1;
        }
      }
    }

    // --- store local values in global matrix and rhs ---
    for(uint i=0; i<FE_T::RefFE_T::numPts; ++i)
    {
      const id_T id_i = e.pointList[i]->id;
      b(id_i) += elemRhs_c(i);
      for(uint j=0; j<FE_T::RefFE_T::numPts; ++j)
      {
        const id_T id_j = e.pointList[j]->id;
        coefficients.push_back(Tri(id_i, id_j, elemMat_c(i,j)));
      }

      // // neumann bc
      // if(id_i == meshPtr->pointList.size()-1)
      // {
      //   b(id_i) += 1.;
      // }
    }
  }
}
