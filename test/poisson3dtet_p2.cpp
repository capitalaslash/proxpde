#include <array>       // for array
#include <bitset>      // for bitset
#include <cmath>       // for sin, M_PI, sqrt
#include <filesystem>  // for path
#include <memory>      // for unique_ptr, construct_at, destr...
#include <optional>    // for optional
#include <string>      // for char_traits, basic_string
#include <tuple>       // for tuple, get
#include <type_traits> // for is_base_of_v
#include <vector>      // for vector

#include <fmt/base.h>      // for println
#include <yaml-cpp/yaml.h> // for Node, LoadFile

#include <proxpde/assembly_lhs.hpp>  // for AssemblyStiffness
#include <proxpde/assembly_rhs.hpp>  // for AssemblyProjection
#include <proxpde/bc.hpp>            // for BCEss
#include <proxpde/builder.hpp>       // for Builder
#include <proxpde/def.hpp>           // for Vec3, ParameterDict, checkError
#include <proxpde/fe.hpp>            // for LagrangeFE
#include <proxpde/fespace.hpp>       // for interpolateAnalyticFunction
#include <proxpde/geo.hpp>           // for Tetrahedron
#include <proxpde/iomanager.hpp>     // for IOManager
#include <proxpde/mesh.hpp>          // for Mesh, MeshType, side, readMesh
#include <proxpde/reffe.hpp>         // for RefTetrahedronP2
#include <proxpde/sparse_matrix.hpp> // for LUSolver
#include <proxpde/timer.hpp>         // for MilliTimer
#include <proxpde/var.hpp>           // for Var

using namespace proxpde;

int test(YAML::Node const & config)
{
  using Elem_T = Tetrahedron;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 2>::RefFE_T,
      LagrangeFE<Elem_T, 2>::RecommendedQR>;

  scalarFun_T const rhs = [](Vec3 const & p)
  {
    return 2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
    // return 0.;
  };

  scalarFun_T const exactSol = [](Vec3 const & p)
  {
    return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1));
    // return 1.;
  };

  MilliTimer t;

  t.start("mesh build");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.start("fe space");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  // face refs with z-axis that exits from the plane, x-axis towards the right
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::vector{bcLeft, bcBottom};
  t.stop();

  t.start("fe build");
  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  // builder.buildRhs(AssemblyRhsAnalytic(rhs, feSpace), bcs);
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
  fmt::println("the norm of the error is {:.16e}", norm);
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
      config["mesh"]["type"] = MeshType::STRUCTURED;
      auto const n = 4U;
      config["mesh"]["n"] = std::array{n, n, n};
      config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
      config["mesh"]["length"] = Vec3{1.0, 1.0, 1.0};
      config["expected_error"] = 0.08097818482375403;
      tests[0] = test(config);
    }
    {
      YAML::Node config;
      config["mesh"]["type"] = MeshType::STRUCTURED;
      auto const n = 8U;
      config["mesh"]["n"] = std::array{n, n, n};
      config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
      config["mesh"]["length"] = Vec3{1.0, 1.0, 1.0};
      config["expected_error"] = 0.02265721054353444;
      tests[1] = test(config);
    }
  }
  return tests.any();
}
