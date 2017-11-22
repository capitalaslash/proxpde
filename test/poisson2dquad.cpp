#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

#include <iostream>

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          FEType<Elem_T,1>::RecommendedQR>;

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
  MilliTimer t;
  uint const numPts_x = (argc < 3)? 11 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 11 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  t.start();
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace(meshPtr);
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList<FESpace_T> bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [] (Vec3 const&) {return 0.;});
  bcs.addEssentialBC(side::BOTTOM, [] (Vec3 const&) {return 0.;});
  std::cout << "bcs: " << t << " ms" << std::endl;


  t.start();
  AssemblyStiffness<FESpace_T> stiffness(1.0, feSpace);
  Builder builder{feSpace.dof.totalNum};
  builder.buildProblem(stiffness, bcs);
  builder.buildProblem(AssemblyAnalyticRhs<FESpace_T>(rhs, feSpace), bcs);
  builder.closeMatrix();
  std::cout << "fe build: " << t << " ms" << std::endl;

  t.start();
  Var sol{"u"};
  Eigen::SparseLU<Mat, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  std::cout << "solve: " << t << " ms" << std::endl;

  Var exact{"exact", feSpace.dof.totalNum};
  interpolateAnalyticFunction(exact_sol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start();
  IOManager<FESpace_T> io{feSpace, "sol_poisson2dquad.xmf"};
  io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 0.020304) > 1.e-5)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
