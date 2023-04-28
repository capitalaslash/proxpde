#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 2>::RefFE_T,
      LagrangeFE<Elem_T, 2>::RecommendedQR>;

  scalarFun_T const rhs = [](Vec3 const & p) { return M_PI * std::sin(M_PI * p(0)); };
  scalarFun_T const exactSol = [](Vec3 const & p)
  { return std::sin(M_PI * p(0)) / M_PI + p(0); };

  MilliTimer t;
  uint const numElems = (argc < 2) ? 20 : std::stoi(argv[1]);

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
  auto bc = BCEss{feSpace, side::LEFT};
  bc << [](Vec3 const &) { return 0.; };
  auto const bcs = std::tuple{bc};
  std::cout << "bcs: " << t << " ms" << std::endl;

  AssemblyStiffness stiffness(1.0, feSpace);
  AssemblyAnalyticRhs f(rhs, feSpace);

  t.start();
  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{stiffness}, bcs);
  builder.buildRhs(std::tuple{f}, bcs);
  builder.closeMatrix();
  std::cout << "fe build: " << t << " ms" << std::endl;

  t.start();
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  std::cout << "solve: " << t << " ms" << std::endl;

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start();
  IOManager io{feSpace, "output/sol_poisson1d_p2"};
  io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm
            << std::endl;
  return checkError({norm}, {3.18934615592152e-07});
}
