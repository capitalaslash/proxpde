#pragma once

#include "def.hpp"
#include "geo.hpp"

void buildProblem(std::shared_ptr<Mesh1D> const meshPtr,
                  std::vector<bc_ess> const& bcs,
                  std::vector<Tri>& coefficients,
                  Vec& b)
{
  for(auto &e: meshPtr->elementList)
  {
    double J = Line::_refVolume / e.volume();

    Line::localVec_T elemRhs = Line::localVec_T::Zero(Line::numPts, 1);

    // A_constrained = C^T A C
    // b_constrained = C^T (b-Ah)
    // C clear constrained rows/cols
    // h is the vector of local constraint values

    Line::localMat_T elemMat_c = J*e.gradMat;
    Line::localVec_T elemRhs_c = elemRhs;

    for(auto& bc: bcs)
    {
      Line::localMat_T C = Line::localMat_T::Identity(Line::numPts, Line::numPts);
      Line::localVec_T h = Line::localVec_T::Zero(Line::numPts, 1);
      for(uint i=0; i<Line::numPts; ++i)
      {
        Point const& p = *e.pointList[i];
        id_T const gid = p.id;
        if(bc.vec[gid])
        {
          h(i) = bc.value(p);
          C(i,i) = 0.;
        }
      }
      elemMat_c = C * elemMat_c * C;
      elemRhs_c = C * (elemRhs_c - J*e.gradMat * h);

      for(uint i=0; i<Line::numPts; ++i)
      {
        Point const& p = *e.pointList[i];
        id_T const gid = p.id;
        if(bc.vec[gid])
        {
          elemMat_c(i,i) = J;
          elemRhs_c(i) = J*h[i];
        }
      }
    }

    for(uint i=0; i<Line::numPts; ++i)
    {
      const id_T id_i = e.pointList[i]->id;
      b(id_i) += elemRhs_c(i);
      for(uint j=0; j<Line::numPts; ++j)
      {
        const id_T id_j = e.pointList[j]->id;
        coefficients.push_back(Tri(id_i, id_j, elemMat_c(i,j)));
      }

      // // neumann bc
      // if(right.vec[id_i])
      // {
      //   b(id_i) += right.value(*e.pointList[i]);
      // }
    }
  }
}
