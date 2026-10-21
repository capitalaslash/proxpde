#include <bitset>      // for bitset
#include <cmath>       // for sin, M_PI, sqrt
#include <cstdlib>     // for uint
#include <filesystem>  // for path, operator/
#include <memory>      // for unique_ptr, construct_at, destr...
#include <optional>    // for optional
#include <string>      // for char_traits, basic_string, string
#include <tuple>       // for tuple, get
#include <type_traits> // for is_base_of_v
#include <vector>      // for vector

#include <fmt/base.h>      // for println
#include <yaml-cpp/yaml.h> // for Node

#include <proxpde/assembly_bc.hpp>   // for AssemblyBCMixed, AssemblyBCNatu...
#include <proxpde/assembly_lhs.hpp>  // for AssemblyStiffness
#include <proxpde/assembly_rhs.hpp>  // for AssemblyRhsAnalytic
#include <proxpde/bc.hpp>            // for BCEss
#include <proxpde/builder.hpp>       // for Builder
#include <proxpde/def.hpp>           // for Vec3, checkError, scalarFun_T
#include <proxpde/fe.hpp>            // for LagrangeFE
#include <proxpde/fespace.hpp>       // for FESpace, interpolateAnalyticFun...
#include <proxpde/geo.hpp>           // for Line
#include <proxpde/iomanager.hpp>     // for IOManager
#include <proxpde/mesh.hpp>          // for side, Mesh, buildHyperCube
#include <proxpde/reffe.hpp>         // for RefLineP1
#include <proxpde/sparse_matrix.hpp> // for LUSolver
#include <proxpde/timer.hpp>         // for operator<<, MilliTimer
#include <proxpde/var.hpp>           // for Var

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

  fmt::println("test setup:\n");
  fmt::println("  - hConv = \n  - temp0 = \n  - tempA = ", hConv, temp0, tempA);

  const scalarFun_T rhsFun = [](Vec3 const & p)
  { return M_PI * std::sin(M_PI * p(0)); };

  const scalarFun_T exactSol = [hConv, temp0, tempA](Vec3 const & p)
  {
    return std::sin(M_PI * p(0)) / M_PI +
           (1. - hConv * (temp0 - tempA)) * p(0) / (hConv + 1.) + temp0;
  };

  auto const numElems = config["n"].as<uint>();

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start("mesh");
  buildHyperCube(*mesh, {0., 0., 0.}, {1., 0., 0.}, {{numElems, 0, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  auto const bcs = std::vector{
      BCEss{feSpace, side::LEFT, [temp0](Vec3 const &) { return temp0; }},
  };
  t.stop();

  t.start("fe build");
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
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  t.start("solve");
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "u:\n" << sol.data << std::endl;

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{
      feSpace, std::filesystem::path{"output"} / config["filename"].as<std::string>()};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  auto const errorNorm = error.data.norm();
  fmt::println("the norm of the error is {:.16e}", errorNorm);
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
