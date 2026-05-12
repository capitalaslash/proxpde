#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Tetrahedron;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

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

  auto const configFile = (argc > 1) ? argv[1] : "poisson3dtet.yaml";
  auto const config = YAML::LoadFile(configFile);

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.start("fespace");
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

  MilliTimer tBuild;

  tBuild.start("builder ctor");
  Builder builder{feSpace.dof.size};
  tBuild.stop();

  tBuild.start("stiffness");
  builder.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  tBuild.stop();

  tBuild.start("rhs");
  // builder.buildRhs(AssemblyRhsAnalytic(rhs, feSpace), bcs);
  // using an interpolated rhs makes its quality independent of the chosen qr
  tBuild.stop();

  tBuild.start("interp");
  Vec rhsProj;
  interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  tBuild.stop();

  tBuild.start("proj");
  builder.buildRhs(std::tuple{AssemblyProjection(1.0, rhsProj, feSpace)}, bcs);
  tBuild.stop();

  tBuild.start("close");
  builder.closeMatrix();
  tBuild.stop();

  tBuild.print();

  t.stop();

  t.start("solve");
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  t.start("check");
  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;
  t.stop();

  t.start("output");
  IOManager io{feSpace, "output_poisson3dtet/sol"};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  double const norm = error.data.norm();
  fmt::println("the norm of the error is {:.16e}", norm);
  return checkError({norm}, {config["expected_error"].as<double>()});
}
