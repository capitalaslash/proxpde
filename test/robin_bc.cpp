#include "def.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

#include <iostream>

using Elem_T = Line;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          FEType<Elem_T,1>::RecommendedQR>;

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
  MilliTimer t;
  uint const numPts = (argc < 2)? 21 : std::stoi(argv[1]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  t.start();
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts, 0, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  // rotation matrix
  double theta = 0. * M_PI / 3.;
  FMat<3,3> R;
  R << std::cos(theta), std::sin(theta), 0.0,
      -std::sin(theta), std::cos(theta), 0.0,
      0.0, 0.0, 1.0;
  auto Rt = R.transpose();

  // rotate mesh
  for (auto & p: meshPtr->pointList)
  {
    p.coord = R * p.coord;
  }

  t.start();
  FESpace_T feSpace(meshPtr);
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList<FESpace_T> bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [](Vec3 const &){return 0.;});
  bcs.addMixedBC(
        side::RIGHT,
        [](Vec3 const &){return 1.;},
        [](Vec3 const &){return 1.;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  AssemblyStiffness<FESpace_T> stiffness(1.0, feSpace);

  auto rotatedRhs = [&Rt] (Vec3 const& p) {return rhs(Rt * p);};
  AssemblyAnalyticRhs<FESpace_T> f(rotatedRhs, feSpace);

  t.start();
  Builder builder{feSpace.dof.size};
  builder.buildProblem(stiffness, bcs);
  builder.buildProblem(f, bcs);
  builder.closeMatrix();
  std::cout << "fe build: " << t << " ms" << std::endl;

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  t.start();
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  std::cout << "solve: " << t << " ms" << std::endl;

  Var exact{"exact"};
  auto rotatedESol = [&Rt] (Vec3 const& p) {return exact_sol(Rt * p);};
  interpolateAnalyticFunction(rotatedESol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start();
  IOManager<FESpace_T> io{feSpace, "sol_robinbc"};
  io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 2.61664e-11) > 1.e-10)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
