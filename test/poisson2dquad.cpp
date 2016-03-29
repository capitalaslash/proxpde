#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "iomanager.hpp"

#include <iostream>

typedef Quad Elem_T;
typedef Mesh<Elem_T> Mesh_T;
typedef FESpace<
          Mesh_T,
          FEType<Elem_T,1>::RefFE_T,
          FEType<Elem_T,1>::RecommendedQR> FESpace_T;

scalarFun_T rhs = [] (Vec3 const& p)
{
  return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};
scalarFun_T exact_sol = [] (Vec3 const& p)
{
  return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};

int main(int argc, char* argv[])
{
  uint const numPts_x = (argc < 3)? 11 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 11 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {numPts_x, numPts_y, 0});

  FESpace_T feSpace(meshPtr);

  bc_ess<FESpace_T> left(feSpace, side::LEFT, [] (Vec3 const&) {return 0.;});
  bc_ess<FESpace_T> bottom(feSpace, side::BOTTOM, [] (Vec3 const&) {return 0.;});
  bc_list<FESpace_T> bcs{feSpace, {left, bottom}};
  bcs.init();

  Mat A(feSpace.dof.totalNum, feSpace.dof.totalNum);
  Vec b = Vec::Zero(feSpace.dof.totalNum);

  AssemblyPoisson<FESpace_T::CurFE_T> assembly(rhs, feSpace.curFE);

  buildProblem(feSpace, assembly, rhs, bcs, A, b);

  Var sol{"u"};
  Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  sol.data = solver.solve(b);

  IOManager<FESpace_T> io{"sol_poisson2dquad.xmf", feSpace};
  io.print({sol});

  Vec exact = Vec::Zero(feSpace.dof.totalNum);
  interpolateAnalyticalFunction(exact_sol, feSpace, exact);
  double norm = (sol.data - exact).norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 0.020304) > 1.e-5)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
