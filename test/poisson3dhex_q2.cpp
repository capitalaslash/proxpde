#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

#include <yaml-cpp/yaml.h>
#include <bitset>

using Elem_T = Hexahedron;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T, 2>::RefFE_T,
                          FEType<Elem_T, 2>::RecommendedQR>;

static scalarFun_T rhs = [] (Vec3 const& p)
{
  return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  // return 2.0;
  // return 0.;
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  // return 2.*p(0) - p(0)*p(0);
  // return 1.;
};

int test(YAML::Node const & config)
{
  MilliTimer t;

  t.start("mesh build");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(
    *mesh,
    {0., 0., 0.},
    {1., 1., 1.},
    {config["n"].as<uint>()+1, config["n"].as<uint>()+1, config["n"].as<uint>()+1});
  t.stop();

  t.start("fe space");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  BCList bcs{feSpace};
  // face refs with z-axis that exits from the plane, x-axis towards the right
  bcs.addEssentialBC(side::LEFT, [] (Vec3 const&) {return 0.;});
  bcs.addEssentialBC(side::BOTTOM, [] (Vec3 const&) {return 0.;});
  t.stop();

  t.start("fe build");
  Builder builder{feSpace.dof.size};
  builder.buildProblem(AssemblyStiffness(1.0, feSpace), bcs);
  builder.buildProblem(AssemblyAnalyticRhs(rhs, feSpace), bcs);
  // using an interpolated rhs makes its quality independent of the chosen qr
  // `Vec rhsProj;
  // `interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  // `builder.buildProblem(AssemblyProjection(1.0, rhsProj, feSpace), bcs);
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
  IOManager io{feSpace, "output_poisson3dhex_q2/sol"};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  if (std::fabs(norm - config["expected_error"].as<double>()) > 1.e-15)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
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
    config["n"] = 1;
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