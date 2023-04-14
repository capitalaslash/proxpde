#include "def.hpp"

#include "assembler.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"
#include "var.hpp"

int main(int argc, char * argv[])
{
  static constexpr uint dim = 2;
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpaceVel_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 2>::RefFE_T,
      LagrangeFE<Elem_T, 2>::RecommendedQR,
      dim>;
  using FESpaceP_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 2>::RecommendedQR>;

  MilliTimer t;

  t.start("config file");

  ParameterDict config;

  // default config
  config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
  config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
  config["mesh"]["n"] = std::array{4U, 4U, 0U};
  config["mesh"]["flags"] = MeshFlags::BOUNDARY_FACETS;
  config["dt"] = 0.1;
  config["ntime"] = 10U;
  config["nu"] = 0.1;
  config["printStep"] = 1U;
  config["toll"] = 1.e-11;

  if (argc > 1)
  {
    config.override(argv[1]);
  }
  config.validate({"mesh", "dt", "ntime", "nu", "printStep", "toll"});
  t.stop();

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh, dim * feSpaceVel.dof.size};
  t.stop();

  t.start("bcs");
  auto zero = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  auto bcsVel = std::tuple{
      BCEss{feSpaceVel, side::BOTTOM},
      BCEss{feSpaceVel, side::RIGHT},
      BCEss{feSpaceVel, side::LEFT},
      BCEss{feSpaceVel, side::TOP},
  };
  std::get<0>(bcsVel) << zero;
  std::get<1>(bcsVel) << zero;
  std::get<2>(bcsVel) << zero;
  std::get<3>(bcsVel) << [](Vec3 const &) { return Vec2{1.0, 0.0}; };

  // select the point(s) on the bottom boundary in the middle
  auto const o = config["mesh"]["origin"].as<Vec3>();
  auto const l = config["mesh"]["length"].as<Vec3>();
  auto const hx = l[0] / config["mesh"]["n"].as<std::array<uint, 3>>()[0];
  auto const xm = o[0] + 0.5 * l[0];
  auto const ymin = o[1];
  DOFCoordSet pinSet{
      feSpaceP,
      [xm, ymin, hx](Vec3 const & p)
      { return std::fabs(p[0] - xm) < hx && std::fabs(p[1] - ymin) < 1e-12; },
  };
  auto bcPin = BCEss{feSpaceP, pinSet.ids};
  bcPin << [](Vec3 const &) { return 0.; };
  auto const bcsP = std::tuple{bcPin};
  t.stop();

  t.start("assembly def");
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU * dim + dofP;
  auto const nu = config["nu"].as<double>();
  auto const dt = config["dt"].as<double>();

  Var vel{"vel", dofU * dim};
  Vec velOld{dofU * dim};
  Var p{"p"};
  p.data = Vec::Zero(dofP);
  AssemblyScalarMass timeDer(1. / dt, feSpaceVel);
  AssemblyAdvection advection(1.0, velOld, feSpaceVel, feSpaceVel);
  AssemblyTensorStiffness diffusion(nu, feSpaceVel);
  // AssemblyStiffness diffusion(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);
  AssemblyProjection timeDerRhs(1. / dt, velOld, feSpaceVel);
  auto const lhs = std::tuple{timeDer, advection, diffusion};
  auto const rhs = std::tuple{timeDerRhs};

  // we need this in order to properly apply the pinning bc on the pressure
  AssemblyDummy dummy{feSpaceP};
  t.stop();

  t.start("ic");
  auto ic = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  interpolateAnalyticFunction(ic, feSpaceVel, vel.data);
  velOld = vel.data;
  t.stop();

  Builder builder{numDOFs};
  Builder builderFixed{numDOFs};

  t.start("assembly fixed");
  builderFixed.buildLhs(std::tuple{timeDer, diffusion}, bcsVel);
  builderFixed.buildCoupling(grad, bcsVel, bcsP);
  builderFixed.buildCoupling(div, bcsP, bcsVel);
  builderFixed.buildLhs(std::tuple{dummy}, bcsP);
  builderFixed.closeMatrix();
  auto const matFixed = builderFixed.A;
  auto const rhsFixed = builderFixed.b;
  t.stop();

  LUSolver solver;
  LUSolver solverFixed;

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_cavitytime/sol_v"};
  ioVel.print({vel});
  IOManager ioP{feSpaceP, "output_cavitytime/sol_p"};
  ioP.print({p});
  t.stop();

  auto const ntime = config["ntime"].as<uint>();
  uint const printStep = config["printStep"].as<uint>();
  double const toll = config["toll"].as<double>();
  Vec sol = Vec::Zero(numDOFs);
  Vec solFixed = Vec::Zero(numDOFs);
  double time = 0.0;
  MilliTimer timerStep;
  for (uint itime = 0; itime < ntime; itime++)
  {
    timerStep.start();
    time += dt;
    std::cout << Utils::separator << "solving timestep " << itime + 1
              << ", time = " << time << std::endl;

    t.start("update");
    velOld = vel.data;
    t.stop();

    t.start("build");
    builder.clear();
    builder.buildLhs(lhs, bcsVel);
    builder.buildCoupling(grad, bcsVel, bcsP);
    builder.buildCoupling(div, bcsP, bcsVel);
    builder.buildLhs(std::tuple{dummy}, bcsP);
    builder.closeMatrix();
    builder.buildRhs(rhs, bcsVel);
    t.stop();

    t.start("build fixed");
    builderFixed.clear();
    builderFixed.buildRhs(std::tuple{timeDerRhs}, bcsVel);
    builderFixed.buildLhs(std::tuple{advection}, bcsVel);
    builderFixed.closeMatrix();
    builderFixed.A += matFixed;
    builderFixed.b += rhsFixed;
    t.stop();

    // auto const diffMat = builder.A - builderFixed.A;
    // std::cout << "diffMat norm: " << diffMat.norm() << std::endl;

    t.start("solve");
    solver.compute(builder.A);
    sol = solver.solve(builder.b);
    vel.data = sol.head(dofU * dim);
    p.data = sol.tail(dofP);
    t.stop();

    t.start("res");
    Vec const res = builder.A * sol - builder.b;
    std::cout << "residual norm: " << res.norm() << std::endl;
    t.stop();

    t.start("solve fixed");
    solverFixed.compute(builderFixed.A);
    solFixed = solverFixed.solve(builderFixed.b);
    t.stop();

    t.start("res fixed");
    Vec const resFixed = builderFixed.A * sol - builderFixed.b;
    std::cout << "residual norm (fixed): " << resFixed.norm() << std::endl;
    t.stop();

    t.start("check");
    auto const solDiffNorm = (sol - solFixed).norm();
    t.stop();

    std::cout << "solution difference norm: " << solDiffNorm << std::endl;
    if (solDiffNorm > toll)
    {
      std::cerr << "the 2 solutions differ" << std::endl;
      return 2;
    }

    t.start("print");
    if (itime % printStep == 0)
    {
      ioVel.print({vel}, time);
      ioP.print({p}, time);
    }
    t.stop();

    std::cout << "time required: " << timerStep << " ms" << std::endl;
  }

  t.start("norm");
  auto const solNorm = sol.norm();
  std::cout << "solution norm: " << std::setprecision(16) << solNorm << std::endl;
  t.stop();

  t.print();

  return checkError({solNorm}, {4.409436658784684}, 1.e-12);
}
