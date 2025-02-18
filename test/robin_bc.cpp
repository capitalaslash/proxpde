#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int test(YAML::Node const & config)
{
  using namespace proxpde;

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  MilliTimer t;

  auto const hConv = config["hConv"].as<double>();
  auto const temp0 = config["temp0"].as<double>();
  auto const tempA = config["tempA"].as<double>();

  std::cout << "test setup:\n"
            << "  - hConv = " << hConv << "\n"
            << "  - temp0 = " << temp0 << "\n"
            << "  - tempA = " << tempA << std::endl;

  const scalarFun_T rhsFun = [](Vec3 const & p)
  { return M_PI * std::sin(M_PI * p(0)); };

  const scalarFun_T exactSol = [hConv, temp0, tempA](Vec3 const & p)
  {
    return std::sin(M_PI * p(0)) / M_PI +
           (1. - hConv * (temp0 - tempA)) * p(0) / (hConv + 1.) + temp0;
  };

  auto const numElems = config["n"].as<uint>();

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start();
  buildHyperCube(*mesh, {0., 0., 0.}, {1., 0., 0.}, {{numElems, 0, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto bc = BCEss{feSpace, side::LEFT};
  bc << [temp0](Vec3 const &) { return temp0; };
  auto const bcs = std::vector{bc};
  std::cout << "bcs: " << t << " ms" << std::endl;

  t.start();
  AssemblyStiffness stiffness{1.0, feSpace};
  AssemblyRhsAnalytic rhsAssembly{rhsFun, feSpace};
  Builder builder{feSpace.dof.size};

  // mixed bc: a u + \nabla u = b
  // - \lap u = f
  // (\nabla u, \nabla v) - <\nabla u, v> = (f, v)
  // (., .) == integration on volume
  // <., .> == integration on boundary
  // (\nabla u, \nabla v) + <a*u, v> = (f, v) + <b, v>
  // a >= \eps > 0 cohercitivity to guarantee solution
  // - \nabla u = hConv (u - tempA)
  // -> a = hConv, b = hConv * tempA
  // hConv -> 0: \nabla u = 0, Neumann homogeneous
  // hConv -> inf: u = b / a = tempA, Dirichlet
  auto const lhs = std::tuple{
      stiffness,
      AssemblyBCMixed{[hConv](Vec3 const &) { return hConv; }, side::RIGHT, feSpace}};
  builder.buildLhs(lhs, bcs);
  auto const rhs = std::tuple{
      rhsAssembly,
      AssemblyBCNaturalAnalytic{
          [hConv, tempA](Vec3 const &) { return hConv * tempA; },
          side::RIGHT,
          feSpace}};
  builder.buildRhs(rhs, bcs);
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

  // std::cout << "u:\n" << sol.data << std::endl;

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start();
  IOManager io{
      feSpace, std::filesystem::path{"output"} / config["filename"].as<std::string>()};
  io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  auto const errorNorm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << errorNorm
            << std::endl;
  return checkError({errorNorm}, {config["expected_error"].as<double>()});
}

int main()
{
  std::bitset<4> tests;
  {
    YAML::Node config;
    config["n"] = 20;
    config["hConv"] = 1.0;
    config["temp0"] = 2.0;
    config["tempA"] = 1.0;
    config["filename"] = "sol_robin_test0";
    config["expected_error"] = 3.106760735256447e-07;

    tests[0] = test(config);
  }

  {
    YAML::Node config;
    config["n"] = 20;
    config["hConv"] = 0.0;
    config["temp0"] = 2.0;
    config["tempA"] = 1.0;
    config["filename"] = "sol_robin_test1";
    config["expected_error"] = 2.877854011280693e-07;

    tests[1] = test(config);
  }

  {
    YAML::Node config;
    config["n"] = 20;
    config["hConv"] = 1e20;
    config["temp0"] = 2.0;
    config["tempA"] = 1.0;
    config["filename"] = "sol_robin_test2";
    config["expected_error"] = 4.261223260744055e-07;

    tests[2] = test(config);
  }

  {
    YAML::Node config;
    config["n"] = 20;
    config["hConv"] = 1.0;
    config["temp0"] = 1.0;
    config["tempA"] = 0.0;
    config["filename"] = "sol_robin_test3";
    config["expected_error"] = 3.106760453810386e-07;

    tests[3] = test(config);
  }

  return tests.any();
}
