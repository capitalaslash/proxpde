#include "def.hpp"

// libfmt
#include <fmt/std.h>

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

using namespace proxpde;

// =====================================================================
struct TestQR
{
  using Real_T = double;
  using GeoElem_T = Quad;
  using Vec_T = Vec2;
  using Weights_T = FVec<4u>;
  static short_T constexpr numPts = 4u;
  static short_T const bestPt;

  static Weights_T const weight;
  static std::array<Vec_T, 4> const node;
};

TestQR::Weights_T const TestQR::weight = TestQR::Weights_T::Constant(1.L);
std::array<TestQR::Vec_T, 4u> const TestQR::node = {{
    Vec2{0.0, -2. / 3},
    Vec2{2. / 3, 0.0},
    Vec2{0.0, 2. / 3},
    Vec2{-2. / 3, 0.0},
}};

short_T const TestQR::bestPt = 0u;

// =====================================================================
using Elem_T = Quad;

template <typename QR>
int test(ParameterDict const & config)
{
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<Mesh_T, LagrangeFE<Elem_T, 1>::RefFE_T, QR>;

  scalarFun_T const rhs = [](Vec3 const & p)
  {
    return 2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
  };

  scalarFun_T const exactSol = [](Vec3 const & p)
  { return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1)); };

  MilliTimer t;

  t.start("mesh build");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, config["mesh"]);
  t.stop();

  t.start("fespace");
  // FESpace_T feSpace{*mesh};
  FESpace_T feSpace;
  feSpace.init(*mesh);
  t.stop();

  t.start("bcs");
  // auto bcLeft = BCEss{feSpace, side::LEFT};
  auto bcLeft = BCEss<FESpace_T>{};
  bcLeft.init(feSpace, side::LEFT);
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::vector{bcLeft, bcBottom};
  t.stop();

  t.start("fe space");
  AssemblyStiffness stiffness(1.0, feSpace);
  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{stiffness}, bcs);
  builder.buildRhs(std::tuple{AssemblyRhsAnalytic(rhs, feSpace)}, bcs);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output_poisson2dquad/sol"};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  double norm = error.data.norm();
  fmt::print("the norm of the error is {:.16e}\n", norm);
  return checkError({norm}, {config["expected_error"].as<double>()});
}

int main(int argc, char * argv[])
{
  if (argc > 1)
  {
    ParameterDict config;
    config = YAML::LoadFile(argv[1]);
    auto const result = test<LagrangeFE<Elem_T, 1u>::RecommendedQR>(config);
    return result;
  }
  else
  {
    std::bitset<6u> tests;

    {
      ParameterDict config;
      config["mesh"]["type"] = MeshType::STRUCTURED;
      config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
      config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
      config["mesh"]["n"] = std::array{10u, 10u, 0u};
      config["expected_error"] = 2.049777877938e-02;
      tests[0] = test<LagrangeFE<Elem_T, 1u>::RecommendedQR>(config);
    }

    {
      ParameterDict config;
      config["mesh"]["type"] = MeshType::STRUCTURED;
      config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
      config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
      config["mesh"]["n"] = std::array{20u, 20u, 0u};
      config["expected_error"] = 9.732187610621e-03;
      tests[1] = test<LagrangeFE<Elem_T, 1u>::RecommendedQR>(config);
    }

    {
      ParameterDict config;
      config["mesh"]["type"] = MeshType::STRUCTURED;
      config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
      config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
      config["mesh"]["n"] = std::array{10u, 10u, 0u};
      config["expected_error"] = 7.925676692262e-03;
      tests[2] = test<TestQR>(config);
    }

    {
      ParameterDict config;
      config["mesh"]["type"] = MeshType::STRUCTURED;
      config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
      config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
      config["mesh"]["n"] = std::array{20u, 20u, 0u};
      config["expected_error"] = 3.914800457522e-03;
      tests[3] = test<TestQR>(config);
    }

    {
      ParameterDict config;
      config["mesh"]["type"] = MeshType::GMSH;
      config["mesh"]["filename"] = "square_q.msh";
      config["expected_error"] = 4.138121542136e-02;
      tests[4] = test<LagrangeFE<Elem_T, 1u>::RecommendedQR>(config);
    }

    {
      ParameterDict config;
      config["mesh"]["type"] = MeshType::GMSH;
      config["mesh"]["filename"] = "square_q.msh";
      config["expected_error"] = 3.987220124118e-02;
      tests[5] = test<TestQR>(config);
    }

    fmt::print("tests: {}\n", tests);
    return tests.any();
  }
}
