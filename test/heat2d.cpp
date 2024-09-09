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

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;
  using FESpaceVel_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR,
      Elem_T::dim>;

  // input parameters ---
  // mesh size
  double const l = 1.0;
  double const h = 10 * l;
  // thermal conductivity
  double const kappa = 1.0;
  // density and heat capacity
  double const rhoCp = 1.0;
  // thermal diffusivity
  double const alpha = kappa / rhoCp;
  // advection velocity
  double const vMean = 1.0;
  Fun<2, 3> const velocity = [vMean](Vec3 const & p)
  {
    return Vec2{0.0, vMean};
    // the parabolic profile in T holds for T_bulk when the velocity is not constant
    // return Vec2{0.0, 1.5 * vMean * (1.0 - p[0] * p[0])};
  };
  // volumetric heat source
  scalarFun_T const qiii = [rhoCp](Vec3 const &) { return 0.0 / rhoCp; };
  // surface heat source
  scalarFun_T const qiiRight = [rhoCp](Vec3 const &) { return -1.0 / rhoCp; };
  // outlet bc
  scalarFun_T const outlet = [qiiRight, vMean, l](Vec3 const & p)
  { return qiiRight(p) / (vMean * l); };
  // inlet temperature
  double const tMInlet = 0.0;
  scalarFun_T const tInlet = [tMInlet](Vec3 const & p)
  { return tMInlet + 1.0 / 6 - 0.5 * p[0] * p[0]; };
  // initial temperature
  double const tInit = 0.0;

  double const dt = 4.0;
  uint const nsteps = 20U;
  // ---

  double const toll = 1.e-6;

  scalarFun_T const exactSol = [tMInlet](Vec3 const & p)
  { return tMInlet + 1.0 / 6 - p[1] - 0.5 * p[0] * p[0]; };

  MilliTimer t;

  uint const numElems = (argc < 2) ? 20U : std::stoi(argv[1]);

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, Vec3{0., 0., 0.}, Vec3{l, h, 0.}, {{numElems, numElems, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  FESpaceVel_T feSpaceVel{*mesh};
  FEVar vel{"vel", feSpaceVel};
  vel << velocity;
  t.stop();

  t.start("bcs");
  auto const bcBottom = BCEss{feSpace, side::BOTTOM, tInlet};
  auto const bcs = std::vector{bcBottom};
  t.stop();

  t.start("assembly def");
  Var temp{"temp"};
  // tempOld can be removed using temp in its place
  Vec tempOld = Vec::Zero(feSpace.dof.size);
  AssemblyMass mass{1.0 / dt, feSpace};
  AssemblyStiffness diffusion{alpha, feSpace};
  AssemblyAdvectionFE advection{vel, feSpace};
  // the lhs terms are the same for full and incremental versions
  auto const lhs = std::tuple{mass, diffusion, advection};
  AssemblyProjection massOld{1.0 / dt, tempOld, feSpace};
  AssemblyAnalyticRhs heatGeneration{qiii, feSpace};
  AssemblyBCNatural heatSurface{qiiRight, side::RIGHT, feSpace};
  AssemblyBCNatural bcOutlet{outlet, side::TOP, feSpace};
  auto const rhs = std::tuple{massOld, heatGeneration, heatSurface, bcOutlet};

  Var tempInc{"tempInc"};
  Var dTemp{"dTemp"};
  // tempOldInc can be removed using temp in its place
  Vec tempOldInc = Vec::Zero(feSpace.dof.size);
  AssemblyProjection massOldInc{1. / dt, tempOldInc, feSpace};
  auto const rhsInc = std::tuple{massOldInc, heatGeneration, heatSurface, bcOutlet};
  t.stop();

  t.start("ic");
  auto const ic = [tInit](Vec3 const &) { return tInit; };
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
  IOManager io{feSpace, "output_heat2d/temp"};
  io.print({temp, tempInc, dTemp});
  IOManager ioVel{feSpaceVel, "output_heat2d/vel"};
  ioVel.print(std::tuple{vel});
  t.stop();

  // TODO: compute tBulk and check it against analytic solution
  FEVar tBulk{"tBulk", feSpace};
  // integrate up to coord x using fe var
  // \Sum 0.5 * (v[l]*temp[l] + v[r]*temp[r]) * dx / \Sum 0.5 * (v[l] + v[r]) * dx

  double time = 0.0;
  for (uint istep = 0; istep < nsteps; ++istep)
  {
    time += dt;
    fmt::print("{}solving timestep {}, time = {}\n", Utils::separator, istep + 1, time);

    t.start("update");
    tempOld = temp.data;
    fmt::print("tempOld norm: {}\n", tempOld.norm());
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
    fmt::print("tempOldInc norm: {}\n", tempOldInc.norm());
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

    fmt::print("stationary check - increment norm: {}\n", dTemp.data.norm());

    t.start("output");
    auto const diffNorm = (temp.data - tempInc.data).norm();
    fmt::print("diffNorm: {}\n", diffNorm);
    if (diffNorm > toll)
    {
      fmt::print(stderr, "the norm of the difference is too big, aborting.\n");
      return 2;
    }
    io.print({temp, tempInc, dTemp}, time);
    t.stop();
  }

  // std::cout << "temp:\n" << temp.data << std::endl;

  Vec exact;
  interpolateAnalyticFunction(exactSol, feSpace, exact);
  Vec error;
  error = temp.data - exact;

  t.print();

  double const norm = error.norm();
  fmt::print("the norm of the error is {:.16e}\n", norm);
  return checkError({norm}, {2.138630383415e-05});
}
