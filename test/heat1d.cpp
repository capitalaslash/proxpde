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

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  scalarFun_T const qiii = [](Vec3 const &) { return 10.; };
  const double k = 1.0;
  // double const qiiLeft = 1.;
  double const tLeft = 1.;
  // double const qiiRight = 0.; // equivalent to no bc
  double const tRight = 0.;
  double const dt = 0.1;
  uint const ntime = 20U;
  double const toll = 1.e-13;

  scalarFun_T const exactSol = [](Vec3 const & p)
  { return -5. * p[0] * p[0] + 4. * p[0] + 1; };

  MilliTimer t;

  uint const numElems = (argc < 2) ? 20U : std::stoi(argv[1]);

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{1., 0., 0.}, {{numElems, 0, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [tLeft](Vec3 const &) { return tLeft; };
  auto bcRight = BCEss{feSpace, side::RIGHT};
  bcRight << [tRight](Vec3 const &) { return tRight; };
  auto const bcs = std::vector{bcLeft, bcRight};
  t.stop();

  t.start("assembly def");
  Var temp{"temp"};
  // tempOld can be removed using temp in its place
  Vec tempOld = Vec::Zero(feSpace.dof.size);
  AssemblyMass mass{1. / dt, feSpace};
  AssemblyStiffness diffusion{k, feSpace};
  // the lhs terms are the same for full and incremental versions
  auto const lhs = std::tuple{mass, diffusion};
  AssemblyRhsAnalytic heatGeneration{qiii, feSpace};
  AssemblyProjection massOld{1. / dt, tempOld, feSpace};
  auto const rhs = std::tuple{massOld, heatGeneration};

  Var tempInc{"tempInc"};
  Var dTemp{"dTemp"};
  // tempOldInc can be removed using temp in its place
  Vec tempOldInc = Vec::Zero(feSpace.dof.size);
  AssemblyProjection massOldInc{1. / dt, tempOldInc, feSpace};
  auto const rhsInc = std::tuple{massOldInc, heatGeneration};
  t.stop();

  t.start("ic");
  auto const ic = [](Vec3 const &) { return 1.; };
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
  double time = 0.;
  IOManager io{feSpace, "output_heat1d/temp"};
  io.print({temp, tempInc, dTemp}, time);
  t.stop();

  for (uint itime = 0; itime < ntime; ++itime)
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
    io.print({temp, tempInc, dTemp}, time);
    t.stop();
  }

  // fmt::println("temp:\n{}", temp.data);

  Vec exact;
  interpolateAnalyticFunction(exactSol, feSpace, exact);
  Vec error;
  error = temp.data - exact;

  t.print();

  double const norm = error.norm();
  fmt::println("the norm of the error is {:.16e}", norm);
  return checkError({norm}, {2.196219431925e-06});
}
