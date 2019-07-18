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

using Elem_T = Tetrahedron;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T, 2>::RefFE_T,
                          FEType<Elem_T, 2>::RecommendedQR>;

static scalarFun_T rhs = [] (Vec3 const& p)
{
  return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  // return 0.;
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  // return 1.;
};

int test(YAML::Node const & config)
{
  MilliTimer t;

  t.start("mesh build");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, config);
  t.stop();

  t.start("fe space");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  // face refs with z-axis that exits from the plane, x-axis towards the right
  auto const bcs = std::make_tuple(
        BCEss{feSpace, side::LEFT,   [] (Vec3 const &) { return 0.; }},
        BCEss{feSpace, side::BOTTOM, [] (Vec3 const &) { return 0.; }});
  t.stop();

  // std::cout << bcs << std::endl;
  // for (auto const & bc: bcs.bcEssList)
  // {
  //   std::cout << "bc on " << bc.marker << std::endl;
  //   for (auto const & id: bc.constrainedDOFSet)
  //   {
  //     std::cout << id << ": " << feSpace.findCoords(id).transpose() << std::endl;
  //   }
  // }

  t.start("fe build");
  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  // builder.buildRhs(AssemblyAnalyticRhs(rhs, feSpace), bcs);
  // using an interpolated rhs makes its quality independent of the chosen qr
  Vec rhsProj;
  interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  builder.buildRhs(std::tuple{AssemblyProjection(1.0, rhsProj, feSpace)}, bcs);
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
  IOManager io{feSpace, "output_poisson3dtet_p2/sol"};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  return checkError({norm}, {config["expected_error"].as<double>()});
}

int main(int argc, char * argv[])
{
  std::bitset<2> tests;
  if (argc > 1)
  {
    auto const config = YAML::LoadFile(argv[1]);
    tests[0] = test(config);
  }
  else
  {
    {
      YAML::Node config;
      config["mesh_type"] = "structured";
      auto const n = 4;
      config["nx"] = n;
      config["ny"] = n;
      config["nz"] = n;
      config["expected_error"] = 0.08097818482375403;
      tests[0] = test(config);
    }
    {
      YAML::Node config;
      config["mesh_type"] = "structured";
      auto const n = 8;
      config["nx"] = n;
      config["ny"] = n;
      config["nz"] = n;
      config["expected_error"] = 0.02265721054353444;
      tests[0] = test(config);
    }
  }
  return tests.any();
}
