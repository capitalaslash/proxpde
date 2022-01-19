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
  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  scalarFun_T const qiii = [](Vec3 const &) { return 0.; };
  const double k = 1.0;
  // double const qiiLeft = 1.;
  double const tLeft = 1.;
  // double const qiiRight = 0.; // equivalent to no bc
  double const tRight = 0.;
  double const dt = 0.1;
  uint const ntime = 10U;
  double const toll = 1.e-13;

  scalarFun_T const exactSol = [](Vec3 const & p) { return 1. - p[0]; };

  MilliTimer t;

  uint const numElems = (argc < 2) ? 20U : std::stoi(argv[1]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start("mesh build");
  buildHyperCube(*mesh, origin, length, {{numElems, 0, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [tLeft](Vec3 const &) { return tLeft; };
  auto bcRight = BCEss{feSpace, side::RIGHT};
  bcRight << [tRight](Vec3 const &) { return tRight; };
  auto const bcs = std::tuple{bcLeft, bcRight};
  // set increment BCs to zero for all Dirichlet boundaries
  auto const zero = [](Vec3 const &) { return 0.; };
  auto bcsInc = bcs;
  static_for(bcsInc, [&zero]([[maybe_unused]] auto const i, auto & bc) { bc << zero; });
  t.stop();

  t.start("fe build");
  Var temp{"temp"};
  Var dTemp{"dTemp"};
  Var tempInc{"tempInc"};
  // tempOldInc can be removed using temp in its place
  Vec tempOld = Vec::Zero(feSpace.dof.size);
  // tempOldInc can be removed using tempInc in its place
  Vec tempOldInc = Vec::Zero(feSpace.dof.size);
  AssemblyMass mass{1. / dt, feSpace};
  AssemblyStiffness diffusion{k, feSpace};
  AssemblyAnalyticRhs heatGeneration{qiii, feSpace};
  AssemblyProjection massOld{1. / dt, tempOld, feSpace};
  AssemblyProjection massOldInc{1. / dt, tempOldInc, feSpace};
  t.stop();

  t.start("ic");
  auto const ic = [](Vec3 const &) { return 1.; };
  interpolateAnalyticFunction(ic, feSpace, temp.data);
  // the ic must be compatible with the original bcs for the incremental solution!
  dTemp.data = Vec::Zero(feSpace.dof.size);
  tempInc.data = temp.data;
  applyBCs(tempInc.data, bcs);
  tempOldInc = tempInc.data;
  // this is not strictly required for the full solution, but it greatly improves the
  // similarity between the solutions
  temp.data = tempInc.data;
  t.stop();

  Builder builder{feSpace.dof.size};
  Builder builderInc{feSpace.dof.size};
  Builder builderTmp{feSpace.dof.size};

  LUSolver solver;

  t.start("print");
  double time = 0.;
  IOManager io{feSpace, "output_heat1d/temp"};
  io.print({temp, tempInc, dTemp}, time);
  t.stop();

  // the lhs terms are the same for full and incremental versions
  auto const lhs = std::tuple{mass, diffusion};
  for (uint itime = 0; itime < ntime; ++itime)
  {
    time += dt;
    std::cout << Utils::separator << "time = " << time << std::endl;

    t.start("build");
    tempOld = temp.data;
    builder.clear();
    builder.buildLhs(lhs, bcs);
    builder.closeMatrix();
    builder.buildRhs(std::tuple{heatGeneration, massOld}, bcs);
    t.stop();

    t.start("build inc");
    tempOldInc += dTemp.data;
    builderInc.clear();
    builderInc.buildLhs(lhs, bcsInc);
    builderInc.closeMatrix();

    // need a matrix where no Dirichlet bcs are set
    // TODO: add buildLHS() split in assembly + apply_bc globally?
    builderTmp.clear();
    builderTmp.buildLhs(lhs, std::tuple{});
    builderTmp.closeMatrix();
    // std::cout << "K:\n" << builderTmp.A << std::endl;
    Vec res = builderTmp.A * tempOldInc;
    // residual on Dirichlet bc must be set to zero
    applyBCs(res, bcsInc);

    builderInc.buildRhs(std::tuple{heatGeneration, massOldInc}, bcsInc);
    builderInc.b -= res;
    t.stop();

    t.start("solve");
    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    temp.data = solver.solve(builder.b);
    t.stop();

    t.start("solve inc");
    solver.analyzePattern(builderInc.A);
    solver.factorize(builderInc.A);
    dTemp.data = solver.solve(builderInc.b);
    // std::cout << "dTemp: " << dTemp.data.transpose() << std::endl;
    t.stop();

    std::cout << "stationary check - increment norm: " << dTemp.data.norm()
              << std::endl;
    t.start("output");
    tempInc.data += dTemp.data;
    auto const diffNorm = (temp.data - tempInc.data).norm();
    std::cout << "diffNorm: " << diffNorm << std::endl;
    assert(diffNorm < toll);
    io.print({temp, tempInc, dTemp}, time);
    t.stop();
  }

  // std::cout << "temp:\n" << temp.data << std::endl;

  Vec exact;
  interpolateAnalyticFunction(exactSol, feSpace, exact);
  Vec error;
  error = temp.data - exact;

  t.print();

  double norm = error.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm
            << std::endl;
  return checkError({norm}, {2.073273836549e-03});
  return 0;
}
