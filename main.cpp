#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "geo.hpp"
#include "bc.hpp"
#include "fe.hpp"

enum SolverType
{
  CHOLESKY,
  BICGSTAB,
  SPARSELU
};

int main()
{
  SolverType solver_type = SPARSELU;

  uint const numPts = 11;
  Vec3 const origin(0., 0., 0.);
  Vec3 const length(1., 0., 0.);

  std::shared_ptr<Mesh1D> meshPtr(new Mesh1D);

  buildMesh1D(meshPtr, origin, length, numPts);

  // bc setup
  bc_ess left, right;
  left.insert(0);
  left.init(numPts);
  left.value = [] (Vec3 const&) {return 0.;};
  right.insert(numPts-1);
  right.init(numPts);
  right.value = [] (Vec3 const&) {return 0.;};

  std::vector<bc_ess> bcs {left};

  std::vector<Tri> coefficients;
  coefficients.reserve((numPts-1)*Line::numPts*Line::numPts);
  Vec b = Vec::Zero(numPts);

  scalarFun_T rhs = [] (Vec3 const& p) {return M_PI*std::sin(M_PI*p(0));};

  buildProblem(meshPtr, rhs, bcs, coefficients, b);

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

  Vec sol;
  switch(solver_type)
  {
    case CHOLESKY:
    {
      Eigen::SimplicialCholesky<Mat> solver(A);
      sol = solver.solve(b);
      break;
    }
    case BICGSTAB:
    {
      Eigen::SimplicialCholesky<Mat> solver(A);
      sol = solver.solve(b);
    }
    case SPARSELU:
    {
      Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
      // Compute the ordering permutation vector from the structural pattern of A
      solver.analyzePattern(A);
      // Compute the numerical factorization
      solver.factorize(A);

      sol = solver.solve(b);
    }
  }
  std::cout<< "sol:\n" << sol << std::endl;

  Vec exact = Eigen::ArrayXd::LinSpaced(numPts, origin(0), origin(0)+length(0));
  for(uint i=0; i<numPts; ++i)
  {
    double x = exact(i);
    exact(i) = std::sin(M_PI*x)/M_PI + x;
  }

  Vec error = sol - exact;
  std::cout << "error: " << error.norm() << std::endl;

  return 0;
}
