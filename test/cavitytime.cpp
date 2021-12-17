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
  using QuadraticRefFE = LagrangeFE<Elem_T, dim>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, dim>::RecommendedQR;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, dim>;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;

  MilliTimer t;

  t.start("config file");

  ParameterDict config;

  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
    config["mesh"]["n"] = std::array{4U, 4U, 0U};
    config["mesh"]["flags"] = MeshFlags::BOUNDARY_FACETS;
    config["dt"] = 0.1;
    config["ntime"] = 10U;
    config["nu"] = 0.1;
    config["printStep"] = 1U;
    config["toll"] = 1.e-11;
  }
  config.validate({"mesh", "dt", "ntime", "nu", "printStep", "toll"});
  t.stop();

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh, feSpaceVel.dof.size * dim};
  t.stop();

  t.start("bcs");
  auto zero = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  auto bcsVel = std::tuple{
      BCEss{feSpaceVel, side::RIGHT},
      BCEss{feSpaceVel, side::LEFT},
      BCEss{feSpaceVel, side::BOTTOM},
      BCEss{feSpaceVel, side::TOP}};
  std::get<0>(bcsVel) << zero;
  std::get<1>(bcsVel) << zero;
  std::get<2>(bcsVel) << zero;
  std::get<3>(bcsVel) << [](Vec3 const &) { return Vec2{1.0, 0.0}; };
  // select the point on the bottom boundary in the middle
  DOFCoordSet pinSet{feSpaceP, [](Vec3 const & p) {
                       return std::fabs(p[0] - 0.5) < 1e-12 && std::fabs(p[1]) < 1e-12;
                     }};
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
  Vec velOld{dofU * dim};

  AssemblyTensorStiffness stiffness(nu, feSpaceVel);
  // AssemblyStiffness stiffness(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);
  AssemblyScalarMass timeder(1. / dt, feSpaceVel);
  AssemblyProjection timeder_rhs(1. / dt, velOld, feSpaceVel);
  AssemblyAdvection advection(1.0, velOld, feSpaceVel, feSpaceVel);
  // we need this in order to properly apply the pinning bc on the pressure
  AssemblyDummy dummy{feSpaceP};
  t.stop();

  t.start("ic");
  Var sol{"vel", numDOFs};
  Vec fixedSol;
  auto ic = [](Vec3 const &) { return Vec2(1., 0.); };
  interpolateAnalyticFunction(ic, feSpaceVel, sol.data);
  fixedSol = sol.data;
  t.stop();

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_cavitytime/sol_v"};
  ioVel.print({sol});
  IOManager ioP{feSpaceP, "output_cavitytime/sol_p"};
  Var p{"p", sol.data, dofU * dim, dofP};
  ioP.print({p});
  t.stop();

  Builder<StorageType::RowMajor> builder{numDOFs};

  t.start("assembly fixed");
  Builder<StorageType::RowMajor> fixedBuilder{numDOFs};
  fixedBuilder.buildLhs(std::tuple{timeder, stiffness}, bcsVel);
  fixedBuilder.buildCoupling(grad, bcsVel, bcsP);
  fixedBuilder.buildCoupling(div, bcsP, bcsVel);
  fixedBuilder.buildLhs(std::tuple{dummy}, bcsP);
  fixedBuilder.closeMatrix();
  auto const fixedMat = fixedBuilder.A;
  auto const fixedRhs = fixedBuilder.b;
  t.stop();

  IterSolver solver;
  IterSolver fixedSolver;
  auto const ntime = config["ntime"].as<uint>();
  uint const printStep = config["printStep"].as<uint>();
  double const toll = config["toll"].as<double>();
  double time = 0.0;
  MilliTimer timerStep;
  auto const lhs = std::tuple{advection, timeder, stiffness};
  auto const rhs = std::tuple{timeder_rhs};
  for (uint itime = 0; itime < ntime; itime++)
  {
    timerStep.start();
    time += dt;
    std::cout << "solving timestep " << itime << ", time = " << time << std::endl;

    t.start("update");
    velOld = sol.data;
    t.stop();

    t.start("assembly");
    builder.clear();
    builder.buildRhs(rhs, bcsVel);
    builder.buildLhs(lhs, bcsVel);
    builder.buildCoupling(grad, bcsVel, bcsP);
    builder.buildCoupling(div, bcsP, bcsVel);
    builder.buildLhs(std::tuple{dummy}, bcsP);
    builder.closeMatrix();
    t.stop();

    t.start("assembly fixed");
    fixedBuilder.clear();
    fixedBuilder.buildRhs(std::tuple{timeder_rhs}, bcsVel);
    fixedBuilder.buildLhs(std::tuple{advection}, bcsVel);
    fixedBuilder.closeMatrix();
    fixedBuilder.A += fixedMat;
    fixedBuilder.b += fixedRhs;
    t.stop();

    // auto const diffMat = builder.A - fixedBuilder.A;
    // std::cout << "diffMat norm: " << diffMat.norm() << std::endl;

    t.start("solve");
    solver.compute(builder.A);
    sol.data = solver.solve(builder.b);
    t.stop();

    t.start("res");
    auto res = fixedBuilder.A * sol.data - fixedBuilder.b;
    std::cout << "residual norm: " << res.norm() << std::endl;
    t.stop();

    t.start("solve fixed");
    fixedSolver.compute(fixedBuilder.A);
    fixedSol = fixedSolver.solve(fixedBuilder.b);
    t.stop();

    t.start("check");
    auto const solDiffNorm = (sol.data - fixedSol).norm();
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
      ioVel.print({sol}, time);
      p.data = sol.data.block(dofU * dim, 0, dofP, 1);
      ioP.print({p}, time);
    }
    t.stop();

    std::cout << "time required: " << timerStep << " ms" << std::endl;
  }

  t.start("norm");
  auto const solNorm = sol.data.norm();
  std::cout << "solution norm: " << std::setprecision(16) << solNorm << std::endl;
  t.stop();

  t.print();

  return checkError({solNorm}, {4.40937006291}, 1.e-10);
}
