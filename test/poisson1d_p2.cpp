#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "iomanager.hpp"

#include <iostream>

typedef Line Elem_T;
typedef Mesh<Elem_T> Mesh_T;
typedef FESpace<
          Mesh_T,
          FEType<Elem_T,2>::RefFE_T,
          FEType<Elem_T,2>::RecommendedQR> FESpace_T;

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

  FESpace_T feSpace(meshPtr);

  bc_ess<FESpace_T> left(feSpace, side::LEFT, [] (Vec3 const&) {return 0.;});
  bc_list<FESpace_T> bcs{feSpace, {left}};
  bcs.init();

  Mat A(feSpace.dof.totalNum, feSpace.dof.totalNum);
  Vec b = Vec::Zero(feSpace.dof.totalNum);

  AssemblyPoisson<FESpace_T::CurFE_T> assembly(rhs, feSpace.curFE);

  buildProblem(feSpace, assembly, rhs, bcs, A, b);
//  std::cout << "A:\n" << A << std::endl;

  Vec sol;
  Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  sol = solver.solve(b);
//  std::cout << "sol:\n" << sol << std::endl;

  IOManager<FESpace_T> io(feSpace);
  io.print(sol);

  Vec exact = Vec::Zero(feSpace.dof.totalNum);
  interpolateAnalyticalFunction(exact_sol, feSpace, exact);
  double norm = (sol - exact).norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 3.18935e-07) > 1.e-9)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
