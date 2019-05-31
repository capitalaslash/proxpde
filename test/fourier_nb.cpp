#include "def.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

#include <yaml-cpp/yaml.h>

using Elem_T = Line;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          FEType<Elem_T,1>::RecommendedQR>;
using RecFESpace_T = FESpace<Mesh_T,
                             FEType<Elem_T, 1>::RefFE_T,
                             FEType<Elem_T, 1>::ReconstructionQR>;

int test(YAML::Node const & config)
{
  MilliTimer t;

  auto const hConv = config["hConv"].as<double>();
  auto const tempA = config["tempA"].as<double>();

  std::cout << "test setup:\n"
            << "  - hConv = " << hConv << "\n"
            << "  - tempA = " << tempA << std::endl;

  const scalarFun_T rhs = [] (Vec3 const &)
  {
    return 2.;
  };

  const scalarFun_T exactSol = [] (Vec3 const & p)
  {
    return 4. - p(0) * p(0);
  };

  const scalarFun_T ic = [tempA] (Vec3 const & p)
  {
    return tempA;
  };

  uint const numElems = config["n"].as<uint>();

  t.start();
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, {0., 0., 0.}, {1., 0., 0.}, {{numElems, 0, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace{*mesh};
  RecFESpace_T feSpaceRec{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList bcs{feSpace};
  // bcs.addBC(BCEss{feSpace, side::LEFT, [temp0](Vec3 const &){return temp0;}});
  // mixed bc: a u + \nabla u = b
  // - \nabla u = hConv (u - tempA)
  // -> a = hConv, b = hConv * tempA
  // hConv -> 0: \nabla u = 0, Neumann homogeneous
  // hConv -> inf: u = b / a = tempA, Dirichlet
  // bcs.addMixedBC(
  //       side::RIGHT,
  //       [hConv] (Vec3 const &) { return hConv; },
  //       [hConv, tempA] (Vec3 const &) { return hConv * tempA; });
  std::cout << "bcs: " << t << " ms" << std::endl;

  uint const size = feSpace.dof.size;
  Var sol{"u"};
  interpolateAnalyticFunction(ic, feSpace, sol.data);
  Vec solOld;
  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  LUSolver solver;
  double const steps = 10;
  double const dt = 2.0;

  AssemblyMass timeDer{1. / dt, feSpace};
  AssemblyStiffness stiffness{1.0, feSpace};
  AssemblyAnalyticRhs f{rhs, feSpace};
  AssemblyProjection timeRhs{1. / dt, solOld, feSpace};
  AssemblyBCMixed mixBC{[hConv] (Vec3 const &) { return hConv; }, side::RIGHT, feSpace};
  AssemblyBCNatural natBC{[hConv, tempA] (Vec3 const &) { return hConv * tempA; }, side::RIGHT, feSpace};

  Var flux{"flux", size};

  Builder builder{size};
  double time = 0.;
  IOManager io{
    feSpace,
    fs::path{"output_fourier"} / config["filename"].as<std::string>()
  };
  io.print({sol, flux});
  for (uint itime = 0; itime < steps; ++itime)
  {
    time += dt;
    std::cout << "solving time " << time << std::endl;

    solOld = sol.data;

    builder.clear();
    builder.buildProblem(timeDer, bcs);
    builder.buildProblem(stiffness, bcs);
    builder.buildProblem(timeRhs, bcs);
    builder.buildProblem(f, bcs);
    builder.buildProblem(natBC, bcs);
    builder.buildProblem(mixBC, bcs);
    builder.closeMatrix();
    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    sol.data = solver.solve(builder.b);
    // std::cout << "u:\n" << sol.data << std::endl;

    reconstructGradient(sol.data, feSpaceRec, flux.data);
    std::cout << "u|1 = " << sol.data[size-1] << ", exact = " << exact.data[size-1] << std::endl;
    std::cout << "du / dx |1 = " << flux.data[size-1] << std::endl;
    std::cout << "h (u|1 - temp0)  = " << hConv * (sol.data[size-1] - tempA) << std::endl;

    io.print({sol, flux}, time);
  }

  t.start();
  std::cout << "fe build: " << t << " ms" << std::endl;

  t.start();
  std::cout << "solve: " << t << " ms" << std::endl;

  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start();
  std::cout << "output: " << t << " ms" << std::endl;

  auto const errorNorm = error.data.norm();
  std::cout << "the norm of the error is "<< std::setprecision(16) << errorNorm << std::endl;
  if(std::fabs(errorNorm - config["expected_error"].as<double>()) > 1.e-15)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }
  return 0;
}

int main()
{
  std::bitset<4> tests;
  {
    YAML::Node config;
    config["n"] = 20;
    config["hConv"] = 1.0;
    config["tempA"] = 1.0;
    config["filename"] = "sol_fourier_test0";
    config["expected_error"] = 0.001390624752302502;

    tests[0] = test(config);
  }

  // {
  //   YAML::Node config;
  //   config["n"] = 20;
  //   config["hConv"] = 0.0;
  //   config["temp0"] = 2.0;
  //   config["tempA"] = 1.0;
  //   config["filename"] = "sol_robin_test1";
  //   config["expected_error"] = 2.877854011280693e-07;
  //
  //   tests[1] = test(config);
  // }
  //
  // {
  //   YAML::Node config;
  //   config["n"] = 20;
  //   config["hConv"] = 1e20;
  //   config["temp0"] = 2.0;
  //   config["tempA"] = 1.0;
  //   config["filename"] = "sol_robin_test2";
  //   config["expected_error"] = 4.261223260744055e-07;
  //
  //   tests[2] = test(config);
  // }
  //
  // {
  //   YAML::Node config;
  //   config["n"] = 20;
  //   config["hConv"] = 1.0;
  //   config["temp0"] = 1.0;
  //   config["tempA"] = 0.0;
  //   config["filename"] = "sol_robin_test3";
  //   config["expected_error"] = 3.106760453810386e-07;
  //
  //   tests[2] = test(config);
  // }

  return tests.any();
}
