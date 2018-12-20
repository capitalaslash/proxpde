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

double test(YAML::Node const & config)
{
  MilliTimer t;

  const double hConv = config["hConv"].as<double>();
  const double temp0 = config["temp0"].as<double>();
  const double tempA = config["tempA"].as<double>();

  std::cout << "test setup:\n"
            << "  - hConv = " << hConv << "\n"
            << "  - temp0 = " << temp0 << "\n"
            << "  - tempA = " << tempA << std::endl;

  const scalarFun_T rhs = [] (Vec3 const& p)
  {
    return M_PI*std::sin(M_PI*p(0));
  };

  const scalarFun_T exactSol = [hConv, temp0, tempA] (Vec3 const& p)
  {
    return std::sin(M_PI*p(0))/M_PI
            + (1. - hConv * (temp0 - tempA)) * p(0)/(hConv + 1.)
                    + temp0;
  };

  uint const numPts = config["n"].as<uint>()+1;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start();
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, origin, length, {{numPts, 0, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [temp0](Vec3 const &){return temp0;});
  // mixed bc: a u + \nabla u = b
  // - \nabla u = hConv (u - tempA)
  // -> a = hConv, b = hConv * tempA
  // hConv -> 0: \nabla u = 0, Neumann homogeneous
  // hConv -> inf: u = b / a = tempA, Dirichlet
  bcs.addMixedBC(
        side::RIGHT,
        [hConv](Vec3 const &){return hConv;},
        [hConv, tempA](Vec3 const &){return hConv * tempA;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  AssemblyStiffness stiffness{1.0, feSpace};
  AssemblyAnalyticRhs f{rhs, feSpace};

  t.start();
  Builder builder{feSpace.dof.size};
  builder.buildProblem(stiffness, bcs);
  builder.buildProblem(f, bcs);
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
  IOManager io{feSpace, fs::path{"output"} / config["filename"].as<std::string>()};
  io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  return error.data.norm();
}

int main()
{
  {
    YAML::Node config;
    config["n"] = 20;
    config["hConv"] = 1.0;
    config["temp0"] = 2.0;
    config["tempA"] = 1.0;
    config["filename"] = "sol_robin_test1";

    auto const error = test(config);
    std::cout << "test1: the norm of the error is "<< std::setprecision(16) << error << std::endl;
    if(std::fabs(error - 3.071207174712583e-11) > 1.e-15)
    {
      std::cerr << "the norm of the error is not the prescribed value" << std::endl;
      return 1;
    }
  }

  {
    YAML::Node config;
    config["n"] = 20;
    config["hConv"] = 0.0;
    config["temp0"] = 2.0;
    config["tempA"] = 1.0;
    config["filename"] = "sol_robin_test2";

    auto const error = test(config);
    std::cout << "test2: the norm of the error is " << error << std::endl;
    if(std::fabs(error - 2.607139656824781e-11) > 1.e-15)
    {
      std::cerr << "the norm of the error is not the prescribed value" << std::endl;
      return 2;
    }
  }

  {
    YAML::Node config;
    config["n"] = 20;
    config["hConv"] = 1e20;
    config["temp0"] = 2.0;
    config["tempA"] = 1.0;
    config["filename"] = "sol_robin_test3";

    auto const error = test(config);
    std::cout << "test3: the norm of the error is " << error << std::endl;
    if(std::fabs(error - 3.753004639986381e-11) > 1.e-15)
    {
      std::cerr << "the norm of the error is not the prescribed value" << std::endl;
      return 3;
    }
  }

  return 0;
}
