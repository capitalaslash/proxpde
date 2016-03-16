#include "def.hpp"
#include "geo.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"

#include <iostream>

typedef Line Elem_T;
typedef Mesh<Elem_T> Mesh_T;
typedef FESpace<
          Mesh_T,
          FEType<Elem_T,1>::RefFE_T,
          GaussQR<Elem_T,3>> FESpace_T;

scalarFun_T rhs = [] (Vec3 const& p)
{
  return M_PI*std::sin(M_PI*p(0));
};
scalarFun_T exact_sol = [] (Vec3 const& p)
{
  return std::sin(M_PI*p(0))/M_PI + p(0);
};

int main(int argc, char* argv[])
{
  uint const numPts = (argc < 2)? 21 : std::stoi(argv[1]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {numPts, 0, 0});

  bc_ess<Mesh_T> left(*meshPtr, side::LEFT, [] (Vec3 const&) {return 0.;});
  bc_list<Mesh_T> bcs{left};
  bcs.init(numPts);

  Mat A(numPts, numPts);
  Vec b = Vec::Zero(numPts);

  FESpace_T feSpace(meshPtr);

  AssemblyPoisson<FESpace_T::CurFE_T> assembly(rhs, feSpace.curFE);

  buildProblem(feSpace, assembly, rhs, bcs, A, b);

  Vec sol;
  Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  sol = solver.solve(b);

  Vec exact = Vec::Zero(numPts);
  for(uint i=0; i<numPts; ++i)
  {
    exact(i) = exact_sol(meshPtr->pointList[i].coord);
  }
  double norm = (sol - exact).norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 2.61664e-11) > 1.e-10)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
