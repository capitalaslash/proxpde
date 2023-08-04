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

  using Elem_T = Tetrahedron;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  scalarFun_T const rhs = [](Vec3 const & p)
  {
    return 2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
    // return 0.;
  };

  scalarFun_T const exactSol = [](Vec3 const & p)
  {
    return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1));
    // return 1.;
  };

  MilliTimer t;

  auto const configFile = (argc > 1) ? argv[1] : "poisson3dtet.yaml";
  auto const config = YAML::LoadFile(configFile);

  t.start();
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, ParameterDict{config["mesh"]});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  // face refs with z-axis that exits from the plane, x-axis towards the right
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::vector{bcLeft, bcBottom};
  std::cout << "bcs: " << t << " ms" << std::endl;

  t.start();
  MilliTimer tBuild;
  tBuild.start();
  Builder builder{feSpace.dof.size};
  std::cout << tBuild << std::endl;
  tBuild.start();
  builder.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  std::cout << tBuild << std::endl;
  tBuild.start();
  // builder.buildRhs(AssemblyAnalyticRhs(rhs, feSpace), bcs);
  // using an interpolated rhs makes its quality independent of the chosen qr
  std::cout << tBuild << std::endl;
  tBuild.start();
  Vec rhsProj;
  interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  std::cout << tBuild << std::endl;
  tBuild.start();
  builder.buildRhs(std::tuple{AssemblyProjection(1.0, rhsProj, feSpace)}, bcs);
  std::cout << tBuild << std::endl;
  tBuild.start();
  builder.closeMatrix();
  std::cout << tBuild << std::endl;
  tBuild.start();
  std::cout << "fe build: " << t << " ms" << std::endl;

  t.start();
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  std::cout << "solve: " << t << " ms" << std::endl;

  t.start();
  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  IOManager io{feSpace, "output_poisson3dtet/sol"};
  io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(12) << norm
            << std::endl;
  return checkError({norm}, {config["expected_error"].as<double>()});
}
