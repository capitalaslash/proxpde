#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T, 0>::RefFE_T,
                          FEType<Elem_T, 0>::RecommendedQR>;

using TimeFun_T = std::function<double(double const)>;

enum class TimeIntegrationMethod: char
{
  BDF1 = 1,
  BDF2 = 2,
};

namespace YAML {
template<>
struct convert<TimeIntegrationMethod> {
  static Node encode(TimeIntegrationMethod const & rhs)
  {
    Node node;
    node.push_back(static_cast<int>(rhs));
    return node;
  }

  static bool decode(Node const & node, TimeIntegrationMethod & rhs)
  {
    if(!node.IsSequence() || node.size() != 1)
    {
      return false;
    }
    rhs = static_cast<TimeIntegrationMethod>(node[0].as<int>());
    return true;
  }
};
}

int test(YAML::Node const & config)
{
  MilliTimer t;

  // du / dt + A u = f
  // u = 1 - exp(-t) -> f = (1 - A) exp(-t) + A

  auto const A = config["A"].as<double>();

  TimeFun_T rhs = [A] (double const t)
  {
    return (1. - A) * std::exp(-t) + A;
  };

  TimeFun_T exactSol = [] (double const t)
  {
    return 1. - std::exp(-t);
  };

  t.start("mesh build");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, {0., 0., 0.}, {1., 1., 0.}, {{1, 1, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  // we don't need any boundary condition
  auto const bcs = std::make_tuple();
  t.stop();

  t.start("fe build");
  Vec uOld;
  Vec uOldOld;
  auto const dt = config["dt"].as<double>();
  AssemblyMass bdf1Lhs{1.0 / dt, feSpace};
  AssemblyMass bdf2Lhs{0.5 / dt, feSpace};
  AssemblyMass mass{A, feSpace};
  AssemblyProjection bdf1Rhs{1.0 / dt, uOld, feSpace};
  AssemblyProjection bdf2Rhs1{1.0 / dt, uOld, feSpace};
  AssemblyProjection bdf2Rhs2{-0.5 / dt, uOldOld, feSpace};
  Builder builder{feSpace.dof.size};
  // builder.buildLhs(timeder, bcs);
  // builder.buildLhs(mass, bcs);
  // builder.closeMatrix();
  t.stop();

  Var u("u", feSpace.dof.size);

  Var exact{"exact"};
  Vec exactOld;
  Var error{"e"};

  interpolateAnalyticFunction([&exactSol] (Vec3 const &) { return exactSol(0.); }, feSpace, exact.data);
  exactOld = exact.data;

  // initial conditions
  u.data = exact.data;
  uOld = u.data;
  uOldOld = u.data;
  error.data = u.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output_timeder/sol"};
  io.print({u, exact, error});
  t.stop();

  LUSolver solver;
  auto const method = config["method"].as<TimeIntegrationMethod>();
  auto const ntime = config["ntime"].as<uint>();
  double time = 0.;
  Vec uTime{ntime};
  Vec exactTime{ntime};
  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    // std::cout << "solving timestep " << itime << ", time = " << time << std::endl;

    uOldOld = uOld;
    uOld = u.data;
    builder.buildLhs(std::tuple{bdf1Lhs, mass}, bcs);
    builder.buildRhs(std::tuple{bdf1Rhs}, bcs);
    if (method == TimeIntegrationMethod::BDF2 && itime > 1)
    {
      builder.buildLhs(std::tuple{bdf2Lhs}, bcs);
      builder.buildRhs(std::tuple{bdf2Rhs1, bdf2Rhs2}, bcs);
    }
    builder.buildRhs(std::tuple{AssemblyAnalyticRhs{[&rhs, time] (Vec3 const &) { return rhs(time); }, feSpace}}, bcs);
    builder.closeMatrix();

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    u.data = solver.solve(builder.b);
    Vec res = builder.A * u.data - builder.b;
    // std::cout << "residual: " << std::setprecision(16) << res.norm() << std::endl;
    // std::cout << "A: " << builder.A << std::endl;
    // std::cout << "b: " << builder.b << std::endl;
    //vstd::cout << "u: " << u.data << std::endl;

    exactOld = exact.data;
    interpolateAnalyticFunction([&exactSol, time] (Vec3 const &) { return exactSol(time); }, feSpace, exact.data);
    error.data = u.data - exact.data;
    uTime[itime] = u.data[0];
    exactTime[itime] = exact.data[0];

    // print
    io.print({u, exact, error}, time);

    builder.clear();
  }

  // t.print();

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  std::cout << "the norm of the time error is " << std::setprecision(16) << (uTime - exactTime).norm() << std::endl;
  if(std::fabs(norm - config["expected error"].as<double>()) > 1.e-12)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}

int main()
{
  std::bitset<6> tests;

  {
    YAML::Node config;
    config["method"] = TimeIntegrationMethod::BDF1;
    config["A"] = 0.5;

    config["dt"] = 0.1;
    config["ntime"] = 100;
    config["expected error"] = 0.0007084789746496511;
    tests[0] = test(config);

    config["dt"] = 0.05;
    config["ntime"] = 200;
    config["expected error"] = 0.0003444821133752329;
    tests[1] = test(config);

    config["dt"] = 0.025;
    config["ntime"] = 400;
    config["expected error"] = 0.0001697833039735475;
    tests[2] = test(config);
  }

  {
    YAML::Node config;
    config["method"] = TimeIntegrationMethod::BDF2;
    config["A"] = 0.5;

    config["dt"] = 0.1;
    config["ntime"] = 100;
    config["expected error"] = 3.555874653571323e-05;
    tests[3] = test(config);

    config["dt"] = 0.05;
    config["ntime"] = 200;
    config["expected error"] = 9.38725893884218e-06;
    tests[4] = test(config);

    config["dt"] = 0.025;
    config["ntime"] = 400;
    config["expected error"] = 2.410583358369855e-06;
    tests[5] = test(config);
  }

  return tests.any();
}
