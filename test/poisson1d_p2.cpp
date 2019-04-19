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

using Elem_T = Line;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,2>::RefFE_T,
                          FEType<Elem_T,2>::RecommendedQR>;

static scalarFun_T rhs = [] (Vec3 const& p)
{
  return M_PI*std::sin(M_PI*p(0));
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return std::sin(M_PI*p(0))/M_PI + p(0);
};

int main(int argc, char* argv[])
{
  MilliTimer t;
  uint const numElems = (argc < 2)? 20 : std::stoi(argv[1]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start();
  buildHyperCube(*mesh, origin, length, {{numElems, 0, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList bcs{feSpace};
  bcs.addBC(BCEss{feSpace, side::LEFT, [] (Vec3 const&) {return 0.;}});
  std::cout << "bcs: " << t << " ms" << std::endl;

  AssemblyStiffness stiffness(1.0, feSpace);
  AssemblyAnalyticRhs f(rhs, feSpace);

  t.start();
  Builder builder{feSpace.dof.size};
  builder.buildProblem(stiffness, bcs);
  builder.buildProblem(f, bcs);
  builder.closeMatrix();
  std::cout << "fe build: " << t << " ms" << std::endl;

  t.start();
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  std::cout << "solve: " << t << " ms" << std::endl;

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start();
  IOManager io{feSpace, "output/sol_poisson1d_p2"};
  io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << norm << std::endl;
  if(std::fabs(norm - 3.18935e-07) > 1.e-9)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
