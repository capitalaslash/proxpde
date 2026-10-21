#include <cmath>       // for sqrt
#include <cstdlib>     // for uint
#include <filesystem>  // for path
#include <memory>      // for construct_at, destroy_at, uniqu...
#include <optional>    // for optional
#include <string>      // for basic_string, stoi
#include <tuple>       // for tuple, get
#include <type_traits> // for is_base_of_v
#include <vector>      // for vector, stderr

#include <fmt/base.h> // for println

#include <proxpde/assembly_lhs.hpp>  // for AssemblyMass, AssemblyStiffness
#include <proxpde/assembly_rhs.hpp>  // for AssemblyProjection, AssemblyRhs...
#include <proxpde/bc.hpp>            // for BCEss
#include <proxpde/builder.hpp>       // for Builder
#include <proxpde/def.hpp>           // for Vec3, Vec, checkError, scalarFun_T
#include <proxpde/fe.hpp>            // for LagrangeFE
#include <proxpde/fespace.hpp>       // for interpolateAnalyticFunction
#include <proxpde/geo.hpp>           // for Line
#include <proxpde/iomanager.hpp>     // for IOManager
#include <proxpde/mesh.hpp>          // for side, Mesh, buildHyperCube
#include <proxpde/reffe.hpp>         // for RefLineP1
#include <proxpde/sparse_matrix.hpp> // for LUSolver
#include <proxpde/timer.hpp>         // for MilliTimer
#include <proxpde/var.hpp>           // for Var

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1u>::RefFE_T,
      LagrangeFE<Elem_T, 1u>::RecommendedQR>;

  scalarFun_T const qiii = [](Vec3 const &) { return 10.0; };
  double const k = 1.0;
  // double const qiiLeft = 1.0;
  double const tLeft = 1.0;
  // double const qiiRight = 0.0; // equivalent to no bc
  double const tRight = 0.0;
  double const dt = 1e-01;
  uint const ntime = 20u;
  double const toll = 1e-12;

  scalarFun_T const exactSol = [](Vec3 const & p)
  { return -5.0 * p[0] * p[0] + 4.0 * p[0] + 1.0; };

  MilliTimer t;

  uint const numElems = (argc < 2) ? 20u : std::stoi(argv[1]);

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, {numElems, 0u, 0u});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  Var exact{"tempExact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);

  t.start("bcs");
  auto const bcs = std::vector{
      BCEss{feSpace, side::LEFT, [tLeft](Vec3 const &) { return tLeft; }},
      BCEss{feSpace, side::RIGHT, [tRight](Vec3 const &) { return tRight; }},
  };
  t.stop();

  t.start("assembly def");
  Var temp{"temp"};
  // tempOld can be removed using temp in its place
  Vec tempOld = Vec::Zero(feSpace.dof.size);
  AssemblyMass mass{1.0 / dt, feSpace};
  AssemblyStiffness diffusion{k, feSpace};
  // the lhs terms are the same for full and incremental versions
  auto const lhs = std::tuple{mass, diffusion};
  AssemblyRhsAnalytic heatGeneration{qiii, feSpace};
  AssemblyProjection massOld{1.0 / dt, tempOld, feSpace};
  auto const rhs = std::tuple{massOld, heatGeneration};

  Var tempInc{"tempInc"};
  Var dTemp{"dTemp"};
  // tempOldInc can be removed using temp in its place
  Vec tempOldInc = Vec::Zero(feSpace.dof.size);
  AssemblyProjection massOldInc{1.0 / dt, tempOldInc, feSpace};
  auto const rhsInc = std::tuple{massOldInc, heatGeneration};
  t.stop();

  t.start("ic");
  auto const ic = [](Vec3 const &) { return 1.0; };
  interpolateAnalyticFunction(ic, feSpace, temp.data);
  // the ic must be compatible with the original bcs for the incremental solution!
  interpolateAnalyticFunction(ic, feSpace, tempInc.data);
  dTemp.data = Vec::Zero(feSpace.dof.size);
  tempOldInc = tempInc.data;
  t.stop();

  Builder builder{feSpace.dof.size};
  Builder builderInc{feSpace.dof.size};

  LUSolver solver;
  LUSolver solverInc;

  t.start("print");
  double time = 0.0;
  IOManager io{feSpace, "output_heat1d/temp"};
  io.print({temp, tempInc, dTemp, exact}, time);
  t.stop();

  for (auto itime = 0u; itime < ntime; ++itime)
  {
    time += dt;
    fmt::println("{}solving timestep {}, time = {}", Utils::separator, itime + 1, time);

    t.start("update");
    tempOld = temp.data;
    fmt::println("tempOld norm: {:e}", tempOld.norm());
    t.stop();

    t.start("build");
    builder.clear();
    builder.buildLhs(lhs, bcs);
    builder.closeMatrix();
    builder.buildRhs(rhs, bcs);
    t.stop();

    t.start("solve");
    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    temp.data = solver.solve(builder.b);
    t.stop();

    t.start("update inc");
    tempOldInc += dTemp.data;
    fmt::println("tempOldInc norm: {:e}", tempOldInc.norm());
    t.stop();

    t.start("build inc");
    builderInc.clear();
    builderInc.buildLhs(lhs, bcs);
    builderInc.closeMatrix();
    builderInc.buildRhs(rhsInc, bcs);
    builderInc.b -= builderInc.A * tempOldInc;
    t.stop();

    t.start("solve inc");
    solverInc.analyzePattern(builderInc.A);
    solverInc.factorize(builderInc.A);
    dTemp.data = solver.solve(builderInc.b);
    tempInc.data += dTemp.data;
    t.stop();

    fmt::println("stationary check - increment norm: {:e}", dTemp.data.norm());
    t.start("output");
    auto const diffNorm = (temp.data - tempInc.data).norm();
    fmt::println("diffNorm: {:e}", diffNorm);
    if (diffNorm > toll)
    {
      fmt::println(stderr, "the norm of the difference is too big, aborting");
      return 2;
    }
    io.print({temp, tempInc, dTemp, exact}, time);
    t.stop();
  }

  // fmt::println("temp:\n{}", temp.data);

  Var error{"tempError"};
  error.data = temp.data - exact.data;

  t.print();

  auto const errorNorm = error.data.norm();
  fmt::println("the norm of the error is {:.16e}", errorNorm);
  return checkError({errorNorm}, {2.196219431925e-06});
}
