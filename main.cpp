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
  right.value = [] (Point const&) {return 0.;};

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

    // A_constrained = C^T A C
    // b_constrained = C^T (b-Ah)
    Eigen::Array<bool,Line::numPts,1> bcs_local =
      Eigen::Array<bool,Line::numPts,1>::Constant(Line::numPts, false);
    Line::localMat_T bc_mask = Eigen::Matrix2d::Identity(2,2);

    for(auto& bc: bcs)
    {
      Eigen::Array<bool,Line::numPts,1> bc_local =
        Eigen::Array<bool,Line::numPts,1>::Constant(Line::numPts, false);
      for(uint i=0; i<Line::numPts; ++i)
      {
        bc_local(i) = bc_local(i) || bc.vec(e.pointList[i]->id);
        bcs_local(i) = bcs_local(i) || bc.vec(e.pointList[i]->id);
      }

      for(uint i=0; i<Line::numPts; ++i)
      {
        if(bc_local[i])
        {
          bc_mask(i,i) = 0;
          for(uint j=0; j<Line::numPts; ++j)
          {
            if(j!=i)
              b[e.pointList[j]->id] -= J*e.gradMat(j,i)*bc.value(Point());
            else
              b[e.pointList[j]->id] = bc.value(Point());
          }
        }
      }
    }

    Line::localMat_T localMat = bc_mask * J*e.gradMat * bc_mask;

    for(uint i=0; i<Line::numPts; ++i)
    {
      if(bcs_local[i])
        localMat(i,i) = 1.0;
    }

    for(uint i=0; i<Line::numPts; ++i)
    {
      const id_T id_i = e.pointList[i]->id;
      for(uint j=0; j<Line::numPts; ++j)
      {
        const id_T id_j = e.pointList[j]->id;
        coefficients.push_back(
          Tri(id_i, id_j, localMat(i,j)) );
      }
    }
  }

  std::cout << "b:\n" << b << std::endl;

  Mat A(m,m);
  A.setFromTriplets(coefficients.begin(), coefficients.end());

  std::cout << "------\nBC" << std::endl;

  // const double penalty = 1.e6;

  // left
  // A.coeffRef(0, 0) += 1.*penalty;
  // b(0) += 10.*penalty;
  // A.coeffRef(0, 0) = 1.;
  // A.coeffRef(0, 1) = 0.;
  // b(0) = 10.;

  // right
  // A.coeffRef(m-1, m-1) += 1.*penalty;
  // b(m-1) += 0.*penalty;
  // A.coeffRef(m-1, m-1) = 1.;
  // A.coeffRef(m-1, m-2) = 0.;

  std::cout << "------\nBC" << std::endl;

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

  // saveAsBitmap(x, n, "output.jpg");

  return 0;
}
