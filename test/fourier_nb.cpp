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

  const scalarFun_T ic = [tempA] (Vec3 const & /*p*/)
  {
    return tempA;
  };

  uint const numElems = config["n"].as<uint>();

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, {0., 0., 0.}, {1., 0., 0.}, {{numElems, 0, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  RecFESpace_T feSpaceRec{*mesh};
  t.stop();

  t.start("bcs");
  auto const bcs = std::make_tuple();
  // bcs.addBC(BCEss{feSpace, side::LEFT, [](Vec3 const &){return 4.;}});
  // bcs.addBC(BCEss{feSpace, side::RIGHT, [](Vec3 const &){return 3.;}});
  t.stop();

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
  // mixed bc: a u + \nabla u = b
  // - \nabla u = hConv (u - tempA)
  // -> a = hConv, b = hConv * tempA
  // hConv -> 0: \nabla u = 0, Neumann homogeneous
  // hConv -> inf: u = b / a = tempA, Dirichlet
  // the matrix block and the rhs block must be added separatly
  AssemblyBCMixed mixBC{
    [hConv] (Vec3 const &) { return hConv; },
    side::RIGHT,
    feSpace};
  AssemblyBCNatural natBC{
    [hConv, tempA] (Vec3 const &) { return hConv * tempA;},
    side::RIGHT,
    feSpace
  };

  Var flux{"flux", size};

  Builder builder{size};
  double time = 0.;
  IOManager io{
    feSpace,
    fs::path{"output_fourier"} / config["filename"].as<std::string>()
  };
  io.print({sol, flux, exact});
  for (uint itime = 0; itime < steps; ++itime)
  {
    time += dt;
    std::cout << "solving time " << time << std::endl;

    solOld = sol.data;

    t.start("build");
    builder.clear();
    builder.buildLhs(std::tuple{timeDer, stiffness, mixBC}, bcs);
    builder.buildRhs(timeRhs, bcs);
    builder.buildRhs(f, bcs);
    builder.buildRhs(natBC, bcs);
    builder.closeMatrix();
    t.stop();
    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;

    t.start("solve");
    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    sol.data = solver.solve(builder.b);
    t.stop();
    // std::cout << "u:\n" << sol.data << std::endl;

    t.start("gradient");
    reconstructGradient(flux.data, feSpaceRec, sol.data, feSpace);
    t.stop();
    std::cout << "u|1 = " << sol.data[size-1] << ", exact = " << exact.data[size-1] << std::endl;
    std::cout << "du / dx |1 = " << flux.data[size-1] << std::endl;
    std::cout << "h (u|1 - temp0)  = " << hConv * (sol.data[size-1] - tempA) << std::endl;

    t.start("print");
    io.print({sol, flux, exact}, time);
    t.stop();
  }

  t.print();

  Var error{"e"};
  error.data = sol.data - exact.data;

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
