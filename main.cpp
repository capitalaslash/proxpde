#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "geo.hpp"
#include "bc.hpp"
#include "fe.hpp"

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

  std::vector<Tri> coefficients;
  coefficients.reserve((numPts-1)*Line::numPts*Line::numPts);
  Vec b = Vec::Zero(numPts);

  buildProblem(meshPtr, bcs, coefficients, b);

  std::cout << "b:\n" << b << std::endl;

  Mat A(numPts,numPts);
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
