#include "def.hpp"

// stl
#include <bitset>

// local
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

using namespace proxpde;

int test(YAML::Node const & config)
{
  using Elem_T = Hexahedron;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  scalarFun_T const rhs = [](Vec3 const & p)
  {
    return 2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
    // return 2.0;
    // return 0.;
  };

  scalarFun_T const exactSol = [](Vec3 const & p)
  {
    return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1));
    // return 2.*p(0) - p(0)*p(0);
    // return 1.;
  };

  MilliTimer t;

  t.start("mesh build");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(
      *mesh,
      {0., 0., 0.},
      {1., 1., 1.},
      {config["n"].as<uint>(), config["n"].as<uint>(), config["n"].as<uint>()});
  t.stop();

  t.start("fe space");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  // face refs with z-axis that exits from the plane, x-axis towards the right
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::vector{bcLeft, bcBottom};
  t.stop();

  t.start("fe build");
  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  builder.buildRhs(std::tuple{AssemblyRhsAnalytic(rhs, feSpace)}, bcs);
  // using an interpolated rhs makes its quality independent of the chosen qr
  // `Vec rhsProj;
  // `interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  // `builder.buildRhs(AssemblyProjection(1.0, rhsProj, feSpace), bcs);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  std::cout << "sol:\n" << sol.data << std::endl;

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output_poisson3dhex/sol"};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  double norm = error.data.norm();

  return checkError({norm}, {config["expected_error"].as<double>()});
}

int main()
{
  // {
  //   YAML::Node config;
  //   config["n"] = 2;
  //   config["expected_error"] = 0.0;
  //   test(config);
  // }
  std::bitset<3> tests;
  {
    YAML::Node config;
    config["n"] = 4;
    config["expected_error"] = 0.1352139631934435;
    tests[0] = test(config);
  }
  // {
  //   YAML::Node config;
  //   config["n"] = 8;
  //   config["expected_error"] = 0.0789237040676312;
  //   tests[1] = test(config);
  // }
  // {
  //   YAML::Node config;
  //   config["n"] = 16;
  //   config["expected_error"] = 0.05080511744764222;
  //   tests[2] = test(config);
  // }
  // {
  //   YAML::Node config;
  //   config["n"] = 32;
  //   config["expected_error"] = 0.0;
  //   test(config);
  // }
  return tests.any();
}
