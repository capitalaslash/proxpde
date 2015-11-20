#pragma once

#include "def.hpp"
#include "geo.hpp"

double rectangleInt(Line const& e, id_T i, scalarFun_T rhs)
{
  // volume * shape fun midpoint * rhs fun midpoint
  return 2.*Line::shapeFuns[i](Vec3::Zero(3))*rhs(e.midpoint());
}

static Eigen::Vector2d const gauss2w = Eigen::Vector2d::Constant(1.L);
static double constexpr sqrt13rd = std::sqrt(1.L/3);
static std::array<Vec3,2> const gauss2n
{
  Vec3(-sqrt13rd, 0., 0.),
  Vec3( sqrt13rd, 0., 0.)
};

double gauss2Int(Line const& e, id_T i, scalarFun_T rhs)
{
  double sum = 0.;
  for(uint g=0; g<2; ++g)
  {
    Vec3 const node = e.midpoint() + 0.5*gauss2n[g]*e.volume();
    sum += gauss2w(g) * Line::shapeFuns[i](gauss2n[g]) * rhs(node);
  }
  return sum;
}

static Eigen::Vector3d const gauss3w = {5.L/9, 8.L/9, 5.L/9};
static double constexpr sqrt35th = std::sqrt(3.L/5);
static std::array<Vec3,3> const gauss3n
{
  Vec3(-sqrt35th, 0., 0.),
  Vec3(       0., 0., 0.),
  Vec3( sqrt35th, 0., 0.)
};

double gauss3Int(Line const& e, id_T i, scalarFun_T rhs)
{
  double sum = 0.;
  for(uint g=0; g<3; ++g)
  {
    Vec3 const node = e.midpoint() + 0.5*gauss3n[g]*e.volume();
    sum += gauss3w(g) * Line::shapeFuns[i](gauss3n[g]) * rhs(node);
  }
  return sum;
}

void buildProblem(std::shared_ptr<Mesh1D> const meshPtr,
                  scalarFun_T const& rhs,
                  bc_list const& bcs,
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
        if(bc.is_constrained(p))
        {
          h(i) = bc.value(p.coord);
          C(i,i) = 0.;
        }
      }
      elemMat_c = C * elemMat_c * C;
      elemRhs_c = C * (elemRhs_c - J*e.gradMat * h);

      for(uint i=0; i<Line::numPts; ++i)
      {
        if(bc.is_constrained(*e.pointList[i]))
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
      if(!bcs.is_constrained(*e.pointList[i]))
      {
        // b(id_i) += rectangleInt(e, i, rhs)/J;
        b(id_i) += gauss3Int(e, i, rhs)/J;
      }
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
