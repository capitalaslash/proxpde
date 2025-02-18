#include "def.hpp"

#include "assembly.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

#include <fmt/ranges.h>

using namespace proxpde;

enum class TimeIntegrationMethod : uint8_t
{
  BDF1 = 1,
  BDF2 = 2,
};

std::string_view TIMtoString(TimeIntegrationMethod m)
{
  switch (m)
  {
  case TimeIntegrationMethod::BDF1:
    return "BDF1";
  case TimeIntegrationMethod::BDF2:
    return "BDF2";
  default:
    abort();
  }
}

namespace YAML
{
template <>
struct convert<TimeIntegrationMethod>
{
  static Node encode(TimeIntegrationMethod const & rhs)
  {
    Node node;
    node.push_back(static_cast<int>(rhs));
    return node;
  }

  static bool decode(Node const & node, TimeIntegrationMethod & rhs)
  {
    if (!node.IsSequence() || node.size() != 1)
    {
      return false;
    }
    rhs = static_cast<TimeIntegrationMethod>(node[0].as<int>());
    return true;
  }
};
} // namespace YAML

template <typename FESpace, TimeIntegrationMethod M>
struct TimeDer
{};

// BDF1: (u - uOld) / dt + A * u = 0 =>
// BDF1: u / dt + A * u = uOld / dt
template <typename FESpace>
struct TimeDer<FESpace, TimeIntegrationMethod::BDF1>
{
  static TimeIntegrationMethod constexpr method = TimeIntegrationMethod::BDF1;

  TimeDer(double idt, FESpace const & feSpace):
      lhs{AssemblyScalarMass{idt, feSpace}},
      rhs{AssemblyProjection{idt, uOld, feSpace},
          AssemblyProjection{0.0, uOld, feSpace}}
  {}

  void update(Vec const & data) { uOld = data; }

  auto & getLhs(uint /*iter*/) { return lhs; }

  auto & getRhs(uint /*iter*/) { return rhs; }

  Vec uOld;
  std::tuple<AssemblyScalarMass<FESpace>> lhs;
  std::tuple<AssemblyProjection<FESpace>, AssemblyProjection<FESpace>> rhs;
};

// BDF2: (3 u - 4 uOld + uOldOld) / 2 dt + A * u = 0 =>
// BDF2: 3 u / 2 dt + A * u = (2 uOld - 1/2 uOldOld ) / dt
template <typename FESpace>
struct TimeDer<FESpace, TimeIntegrationMethod::BDF2>
{
  static TimeIntegrationMethod constexpr method = TimeIntegrationMethod::BDF2;

  TimeDer(double idt, FESpace const & feSpace):
      lhs{AssemblyScalarMass{1.5 * idt, feSpace}},
      rhs{AssemblyProjection{2.0 * idt, uOld, feSpace},
          AssemblyProjection{-0.5 * idt, uOldOld, feSpace}},
      lhsBdf1{AssemblyScalarMass{idt, feSpace}},
      rhsBdf1{
          AssemblyProjection{idt, uOld, feSpace},
          AssemblyProjection{0.0, uOldOld, feSpace}}
  {}

  void update(Vec const & data)
  {
    uOldOld = uOld;
    uOld = data;
  }

  auto & getLhs(uint iter)
  {
    if (iter < 1U)
    {
      return lhsBdf1;
    }
    else
    {
      return lhs;
    }
  }

  auto & getRhs(uint iter)
  {
    if (iter < 1U)
    {
      return rhsBdf1;
    }
    else
    {
      return rhs;
    }
  }

  Vec uOld;
  Vec uOldOld;
  std::tuple<AssemblyScalarMass<FESpace>> lhs;
  std::tuple<AssemblyProjection<FESpace>, AssemblyProjection<FESpace>> rhs;
  std::tuple<AssemblyScalarMass<FESpace>> lhsBdf1;
  std::tuple<AssemblyProjection<FESpace>, AssemblyProjection<FESpace>> rhsBdf1;
};

template <TimeIntegrationMethod Method>
std::tuple<bool, double> test(YAML::Node const & config)
{
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 0>::RefFE_T,
      LagrangeFE<Elem_T, 0>::RecommendedQR>;

  MilliTimer t;

  // du / dt + A u = f
  // u = 1 - exp(-t) -> f = (1 - A) exp(-t) + A
  auto const A = config["A"].as<double>();
  auto const rhs = [A](double const t) { return (1. - A) * std::exp(-t) + A; };
  auto const exactSol = [](double const t) { return 1. - std::exp(-t); };

  t.start("mesh build");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, {0., 0., 0.}, {1., 1., 0.}, {{1, 1, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  // we don't need any boundary condition
  t.stop();

  t.start("ode setup");
  auto const dt = config["dt"].as<double>();
  TimeDer<FESpace_T, Method> timeDer{1.0 / dt, feSpace};
  AssemblyScalarMass reaction{A, feSpace};
  Builder builder{feSpace.dof.size};
  t.stop();

  Var u("u", feSpace.dof.size);

  Var exact{"exact"};
  interpolateAnalyticFunction(
      [&exactSol](Vec3 const &) { return exactSol(0.0); }, feSpace, exact.data);

  // initial conditions
  u.data = exact.data;
  timeDer.update(u.data);
  Var error{"e"};
  error.data = u.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output_timeder/sol"};
  io.print({u, exact, error});
  t.stop();

  LUSolver solver;
  auto const nSteps = config["n_steps"].as<uint>();
  Vec uTime{nSteps};
  Vec exactTime{nSteps};

  fmt::print("\ntime discretization method: {}, ", TIMtoString(Method));
  fmt::print("n steps: {}\n", nSteps);

  double time = 0.0;
  for (uint itime = 0; itime < nSteps; itime++)
  {
    time += dt;
    // fmt::print("solving timestep {}, time = {}\n", itime, time);

    timeDer.update(u.data);
    auto const fTime =
        AssemblyRhsAnalytic{[&rhs, time](Vec3 const &) { return rhs(time); }, feSpace};
    builder.buildLhs(std::tuple_cat(timeDer.getLhs(itime), std::tuple{reaction}));
    builder.buildRhs(std::tuple_cat(timeDer.getRhs(itime), std::tuple{fTime}));
    builder.closeMatrix();

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    u.data = solver.solve(builder.b);
    // Vec res = builder.A * u.data - builder.b;
    // std::cout << "residual: " << std::setprecision(16) << res.norm() << std::endl;
    // std::cout << "A: " << builder.A << std::endl;
    // std::cout << "b: " << builder.b << std::endl;
    // std::cout << "u: " << u.data << std::endl;

    interpolateAnalyticFunction(
        [&exactSol, time](Vec3 const &) { return exactSol(time); },
        feSpace,
        exact.data);
    error.data = u.data - exact.data;
    uTime[itime] = u.data[0];
    exactTime[itime] = exact.data[0];

    // print
    io.print({u, exact, error}, time);

    builder.clear();
  }

  // t.print();

  double norm = error.data.norm();
  fmt::print("the norm of the error is {:.16e}\n", norm);
  fmt::print("the norm of the time error is {:.16e}\n", (uTime - exactTime).norm());
  return {checkError({norm}, {config["expected error"].as<double>()}), norm};
}

int main()
{
  uint constexpr testSize = 4U;
  std::bitset<2 * testSize> tests;

  auto fout = std::fopen("timeder.py", "w");
  fmt::print(fout, "from matplotlib import pyplot\n\n");

  double const totalTime = 10.0;
  YAML::Node config;
  config["A"] = 0.5;

  std::vector<double> const bdf1ExpectedErrors = {
      7.0847897464965115e-04,
      3.4448211337523293e-04,
      1.6978330397354746e-04,
      8.4274977572129650e-05,
  };

  std::vector<double> const bdf2ExpectedErrors = {
      4.2604048825500840e-06,
      1.2720539949162557e-06,
      3.4391182679449628e-07,
      8.9218863963402839e-08,
  };

  std::vector<double> dts(testSize);
  std::vector<double> bdf1Errors(testSize);
  std::vector<double> bdf2Errors(testSize);

  {
    uint nSteps = 100;
    for (uint i = 0U; i < testSize; i++)
    {
      dts[i] = totalTime / nSteps;
      config["dt"] = totalTime / nSteps;
      config["n_steps"] = nSteps;

      config["expected error"] = bdf1ExpectedErrors[i];
      auto const [success1, error1] = test<TimeIntegrationMethod::BDF1>(config);
      tests[i] = success1;
      bdf1Errors[i] = error1;

      config["expected error"] = bdf2ExpectedErrors[i];
      auto const [success2, error2] = test<TimeIntegrationMethod::BDF2>(config);
      tests[i + testSize] = success2;
      bdf2Errors[i] = error2;

      nSteps *= 2;
    }
  }

  fmt::print(fout, "dt = {}\n", dts);
  fmt::print(fout, "bdf1 = {::.12e}\n\n", bdf1Errors);
  fmt::print(fout, "slope1 = [{:.12e}]\n", bdf1Errors[0] * 0.8);
  fmt::print(fout, "for i in range({}):\n", testSize - 1);
  fmt::print(fout, "    slope1.append(0.5 * slope1[i])\n\n");
  fmt::print(fout, "pyplot.loglog(dt, bdf1, 'o-', label='BDF1')\n");
  fmt::print(fout, "pyplot.loglog(dt, slope1, 'k-', label='slope1')\n\n");
  fmt::print(fout, "bdf2 = {::.12e}\n\n", bdf2Errors);
  fmt::print(fout, "slope2 = [{:.12e}]\n", bdf2Errors[0] * 0.8);
  fmt::print(fout, "for i in range({}):\n", testSize - 1);
  fmt::print(fout, "    slope2.append(0.25 * slope2[i])\n\n");
  fmt::print(fout, "pyplot.loglog(dt, bdf2, '^-', label='BDF2')\n");
  fmt::print(fout, "pyplot.loglog(dt, slope2, 'k--', label='slope2')\n\n");
  fmt::print(fout, "pyplot.legend()\n");
  fmt::print(fout, "pyplot.show()\n");

  if (tests.any())
    fmt::print(std::cerr, "tests: {}\n", tests.to_string());
  return tests.any();
}
