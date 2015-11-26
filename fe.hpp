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
    curFE.reinit(e);

    // double const J = curFE.J;

    Line::localVec_T elemRhs = Line::localVec_T::Zero(Line::numPts, 1);

    // A_constrained = C^T A C
    // b_constrained = C^T (b-Ah)
    // C clear constrained rows/cols
    // h is the vector of local constraint values

    Line::localMat_T elemMat_c = curFE.stiffMat;
    Line::localVec_T elemRhs_c = elemRhs;

    for(auto& bc: bcs)
    {
      Line::localMat_T C = Line::localMat_T::Identity(Line::numPts, Line::numPts);
      Line::localVec_T h = Line::localVec_T::Zero(Line::numPts, 1);
      for(uint i=0; i<Line::numPts; ++i)
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

      for(uint i=0; i<Line::numPts; ++i)
      {
        if(bc.is_constrained(*e.pointList[i]))
        {
          elemMat_c(i,i) = 1./curFE.J;
          elemRhs_c(i) = h[i]/curFE.J;
        }
      }
    }

    for(uint i=0; i<Line::numPts; ++i)
    {
      const id_T id_i = e.pointList[i]->id;
      b(id_i) += elemRhs_c(i);
      if(!bcs.is_constrained(*e.pointList[i]))
      {
        // b(id_i) += J * rectangleInt(e, i, rhs);
        // b(id_i) += J * gauss3Int(e, curFE.phi[i], rhs);
        for(uint q=0; q<FE_T::QR_T::numPts; ++q)
        {
          double const f = rhs(curFE.qpoint[q]);
          b(id_i) += curFE.JxW[q]*curFE.phi(i, q) * f;
        }
      }
    }
    for(uint i=0; i<Line::numPts; ++i)
    {
      const id_T id_i = e.pointList[i]->id;
      for(uint j=0; j<Line::numPts; ++j)
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
