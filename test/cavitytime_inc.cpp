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
  std::get<3>(bcsVel) << zero;
  // std::get<3>(bcsVel) << [](Vec3 const &) { return Vec2{1.0, 0.0}; };
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
  Var vel{"vel", dofU * dim};
  Var dvel{"dvel", dofU * dim};
  Vec resVel{dofU * dim};

  AssemblyTensorStiffness diffusion(nu, feSpaceVel);
  // AssemblyStiffness stiffness(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);
  AssemblyScalarMass timeDer(1. / dt, feSpaceVel);
  AssemblyProjection timeDerRhs(1. / dt, vel.data, feSpaceVel);
  AssemblyProjection residual(1.0, resVel, feSpaceVel);
  AssemblyAdvection advection(1.0, vel.data, feSpaceVel, feSpaceVel);
  // we need this in order to properly apply the pinning bc on the pressure
  AssemblyDummy dummy{feSpaceP};
  t.stop();

  t.start("ic");
  Vec sol = Vec::Zero(numDOFs);
  Vec solFixed;
  auto ic = [](Vec3 const &) { return Vec2{1., 0.}; };
  interpolateAnalyticFunction(ic, feSpaceVel, vel.data);
  t.stop();

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_cavitytime_inc/sol_v"};
  ioVel.print({vel, dvel});
  IOManager ioP{feSpaceP, "output_cavitytime_inc/sol_p"};
  Var p{"p", sol, dofU * dim, dofP};
  ioP.print({p});
  t.stop();

  Builder<StorageType::RowMajor> builder{numDOFs};

  t.start("assembly fixed");
  Builder<StorageType::RowMajor> fixedBuilder{numDOFs};
  fixedBuilder.buildLhs(std::tuple{timeDer, diffusion}, bcsVel);
  fixedBuilder.buildCoupling(grad, bcsVel, bcsP);
  fixedBuilder.buildCoupling(div, bcsP, bcsVel);
  fixedBuilder.buildLhs(std::tuple{dummy}, bcsP);
  fixedBuilder.closeMatrix();
  auto const fixedMat = fixedBuilder.A;
  auto const fixedRhs = fixedBuilder.b;
  t.stop();

  IterSolver solver;
  IterSolver solverFixed;
  auto const ntime = config["ntime"].as<uint>();
  uint const printStep = config["printStep"].as<uint>();
  double const toll = config["toll"].as<double>();
  double time = 0.0;
  MilliTimer timerStep;
  auto const lhs = std::tuple{timeDer, advection, diffusion};
  auto const rhs = std::tuple{residual};
  for (uint itime = 0; itime < ntime; itime++)
  {
    timerStep.start();
    time += dt;
    std::cout << "solving timestep " << itime << ", time = " << time << std::endl;

    t.start("update");
    vel.data += sol.block(0, 0, dofU * dim, 1);
    t.stop();

    t.start("assembly");
    builder.clear();
    builder.buildLhs(lhs, bcsVel);
    builder.buildCoupling(grad, bcsVel, bcsP);
    builder.buildCoupling(div, bcsP, bcsVel);
    builder.closeMatrix();
    resVel = vel.data / dt - builder.A.block(0, 0, dofU * dim, dofU * dim) * vel.data;
    std::cout << "rhs norm: " << resVel.norm() << std::endl;
    builder.buildRhs(rhs, bcsVel);
    builder.buildLhs(std::tuple{dummy}, bcsP);
    t.stop();

    t.start("assembly fixed");
    fixedBuilder.clear();
    fixedBuilder.buildRhs(std::tuple{timeDerRhs}, bcsVel);
    fixedBuilder.buildLhs(std::tuple{advection}, bcsVel);
    fixedBuilder.closeMatrix();
    fixedBuilder.A += fixedMat;
    fixedBuilder.b += fixedRhs;
    t.stop();

    // auto const diffMat = builder.A - fixedBuilder.A;
    // std::cout << "diffMat norm: " << diffMat.norm() << std::endl;

    t.start("solve");
    solver.compute(builder.A);
    sol = solver.solve(builder.b);
    dvel.data = sol.block(0, 0, dofU * dim, 1);
    vel.data += dvel.data;
    t.stop();

    std::cout << "residual norm: " << sol.norm() << std::endl;

    t.start("solve fixed");
    solverFixed.compute(fixedBuilder.A);
    solFixed = solverFixed.solve(fixedBuilder.b);
    t.stop();

    t.start("check");
    Vec tmp{numDOFs};
    tmp = sol;
    tmp.block(0, 0, dofU * dim, 1) += vel.data;
    auto const solDiffNorm = (tmp - solFixed).norm();
    t.stop();

    std::cout << "solution difference norm: " << solDiffNorm << std::endl;
    // if (solDiffNorm > toll)
    // {
    //   std::cerr << "the 2 solutions differ" << std::endl;
    //   return 2;
    // }

    t.start("print");
    if (itime % printStep == 0)
    {
      ioVel.print({vel, dvel}, time);
      p.data = sol.block(dofU * dim, 0, dofP, 1);
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

  return checkError({solNorm}, {4.40937006291}, 1.e-10);
}
