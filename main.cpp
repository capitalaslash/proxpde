#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "geo.hpp"
#include "bc.hpp"
// #include "poisson2d.hpp"

void buildMesh1D(std::shared_ptr<Mesh1D> meshPtr,
                   Point const& origin,
                   Point const& length,
                   uint const numPts)
{
  Point const h(length.coord / (numPts-1));
  meshPtr->pointList.reserve(numPts);
  for(uint p=0; p<numPts; ++p)
  {
    meshPtr->pointList.emplace_back(origin.coord + p * h.coord, p);
  }

  uint const numElems = numPts-1;
  meshPtr->elementList.reserve(numElems);
  for(uint e=0; e<numElems; ++e)
  {
    meshPtr->elementList.emplace_back(
      std::array<Point*,Line::numPts>{
        &meshPtr->pointList[e], &meshPtr->pointList[e+1]});
  }

  meshPtr->buildConnectivity();
}

int main()
{
  uint const numPts = 3;
  Point const origin(0., 0., 0.);
  Point const length(1., 0., 0.);

  std::shared_ptr<Mesh1D> meshPtr(new Mesh1D);

  buildMesh1D(meshPtr, origin, length, numPts);

  // bc setup
  bc_ess left, right;
  left.insert(0);
  left.init(numPts);
  left.value = [] (Point const&) {return 10.;};
  right.insert(numPts-1);
  right.init(numPts);
  right.value = [] (Point const&) {return 5.;};

  std::vector<bc_ess> bcs {left, right};

  // uint const n = 5;
  // uint const m = n*n;
  uint const m = numPts;

  std::vector<Tri> coefficients;
  // coefficients.reserve(lines.size()*Line::numPts*Line::numPts);
  Vec b = Vec::Zero(m);

  //buildProblem(coefficients, b, n);

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
          elemMat_c(i,i) = 1.;
          elemRhs_c(i) = h[i];
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
        coefficients.push_back(
          Tri(id_i, id_j, elemMat_c(i,j)) );
      }
    }
  }

  std::cout << "b:\n" << b << std::endl;

  Mat A(m,m);
  A.setFromTriplets(coefficients.begin(), coefficients.end());
  std::cout << "A:\n" << A << std::endl;

  // std::ofstream fout("output.m");
  // for( int k=0; k<A.outerSize(); k++)
  // {
  //   for (Mat::InnerIterator it(A,k); it; ++it)
  //   {
  //     std::cout << it.row() << " " << it.col() << " " << it.value() << " " << it.index() << std::endl;
  //     fout << it.row()+1 << " " << it.col()+1 << " " << it.value() << std::endl;
  //   }
  //   std::cout << "-----" << std::endl;
  // }
  // std::cout << "=====" << std::endl;
  // fout.close();

  // Eigen::SimplicialCholesky<Mat> solver(A);
  Eigen::BiCGSTAB<Mat> solver(A);
  Vec x = solver.solve(b);

  std::cout<< "x:\n" << x << std::endl;

  return 0;
}
