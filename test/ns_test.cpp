#include "def.hpp"

#include "ns.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  uint const dim = Elem_T::dim;
  using Mesh_T = Mesh<Elem_T>;

  MilliTimer t;
  ParameterDict config;

  // default configuration
  config["mesh"]["origin"] = Vec3{0., 0., 0.};
  config["mesh"]["length"] = Vec3{1., 10., 0.};
  config["mesh"]["n"] = std::array{4u, 8u, 0u};
  config["ntime"] = 50u;
  config["ns"]["dt"] = 0.1;
  config["ns"]["nu"] = 0.1;
  config["ns"]["output_dir"] = "output_ns/monolithic";

  if (argc > 1)
  {
    // override from command line
    config.override(argv[1]);
  }

  config.validate({"mesh", "ns", "ntime"});

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new (Mesh_T)};
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.start("eqn monolithic");
  NSSolverMonolithic ns{*mesh, ParameterDict{config["ns"]}};
  t.stop();

  config["ns"]["output_dir"] = "output_ns/split";

  t.start("eqn split");
  NSSolverSplit split{*mesh, ParameterDict{config["ns"]}};
  t.stop();

  t.start("bc");
  auto const zero2d = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  auto const zero = [](Vec3 const &) { return 0.0; };
  auto const inlet = [](Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p(0) * p(0))}; };
  // auto const inlet = [](Vec3 const &) { return Vec2{0.0, 1.0}; };
  auto const bcsVel = std::vector{
      BCEss{ns.feSpaceVel, side::RIGHT, zero2d},
      // BCEss{ns.feSpaceVel, side::TOP, zero2d, Comp::x},
      BCEss{ns.feSpaceVel, side::LEFT, zero2d, Comp::x},
      BCEss{ns.feSpaceVel, side::BOTTOM, inlet},
  };

  // // select the point(s) on the top boundary in the middle
  // auto const o = config["mesh"]["origin"].as<Vec3>();
  // auto const l = config["mesh"]["length"].as<Vec3>();
  // auto const hx = l[0] / config["mesh"]["n"].as<std::array<uint, 3>>()[0];
  // auto const xm = o[0] + 0.5 * l[0];
  // auto const ymax = o[1] + l[1];
  // DOFCoordSet pinSet{
  //     ns.feSpaceP,
  //     [hx, xm, ymax](Vec3 const & p)
  //     { return std::fabs(p[0] - xm) < hx && std::fabs(p[1] - ymax) < 1.e-12; },
  // };
  // auto bcPin = BCEss{ns.feSpaceP, pinSet.ids};
  // bcPin << [](Vec3 const &) { return 0.; };
  // auto const bcsP = std::vector{bcPin};

  auto const bcsPSplit = std::vector{
      BCEss{split.feSpaceP, side::TOP, zero},
  };
  t.stop();

  t.start("monolithic init");
  ns.init(bcsVel);
  // when monolithic, fixed pressure values MUST be imposed weakly via Neumann bcs
  ns.assemblyBCNeumann.marker = side::TOP;
  ns.assemblyBCNeumann.rhs = [](Vec3 const &) { return -1.0; };
  t.stop();

  t.start("split init");
  // no bcs for the velocity, they are derived from velStar
  split.init(NSSolverSplit<Mesh_T>::BCVelList_T{}, bcsPSplit);
  t.stop();

  t.start("ic");
  auto const ic = [](Vec3 const &) { return Vec2{0.0, 1.0}; };
  // auto const ic = [] (Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p(0)*p(0))}; };
  ns.ic(ic);
  for (auto k = 0u; k < ns.feSpaceP.dof.size; k++)
    ns.sol.data[ns.feSpaceVel.dof.size * 2 + k] = 1.0;
  split.ic(ic);
  t.stop();

  t.start("print");
  ns.print();
  split.print();
  t.stop();

  auto const ntime = config["ntime"].as<uint>();
  double time = 0.0;
  for (uint itime = 0; itime < ntime; ++itime)
  {
    MilliTimer stepTimer;
    stepTimer.start();

    time += ns.config["dt"].as<double>();
    std::cout << Utils::separator << "solving timestep " << itime + 1
              << ", time = " << time << std::endl;

    t.start("monolithic assembly");
    ns.assemblyStep(bcsVel);
    t.stop();

    t.start("monolithic solve");
    ns.solve();
    t.stop();

    t.start("split assembly");
    split.assemblyStepVelStar(bcsVel);
    t.stop();

    t.start("split solve");
    split.solveVelStar<VelStarSolverType::MONOLITHIC>();
    t.stop();

    t.start("split assembly");
    split.assemblyStepP(bcsPSplit);
    t.stop();

    t.start("split solve");
    split.solveP();
    t.stop();

    t.start("split assembly");
    split.assemblyStepVel();
    t.stop();

    t.start("split solve");
    split.solveVel();
    t.stop();

    t.start("print");
    ns.print(time);
    split.print(time);
    t.stop();

    std::cout << "time required: " << stepTimer << " ms" << std::endl;
  }

  t.print();

  Var errorVel{"errorU"};
  errorVel.data = ns.sol.data.head(ns.feSpaceVel.dof.size * dim) - split.vel.data;
  Var errorP{"errorP"};
  errorP.data = ns.p.data - split.p.data;

  IOManager ioErrorVel{split.feSpaceVel, "output_ns/eu"};
  ioErrorVel.print({errorVel});
  IOManager ioErrorP{split.feSpaceP, "output_ns/ep"};
  ioErrorP.print({errorP});

  using FESpaceComp_T = Scalar_T<NSSolverMonolithic<Mesh_T>::FESpaceVel_T>;
  FESpaceComp_T const feSpaceComp{*mesh};

  Vec errorU{feSpaceComp.dof.size};
  getComponent(errorU, feSpaceComp, errorVel.data, ns.feSpaceVel, 0);
  Vec errorV{feSpaceComp.dof.size};
  getComponent(errorV, feSpaceComp, errorVel.data, ns.feSpaceVel, 1);

  double const errorNormU = errorU.norm();
  double const errorNormV = errorV.norm();
  double const errorNormP = errorP.data.norm();
  std::cout << "errorU: " << std::setprecision(16) << errorNormU << std::endl;
  std::cout << "errorV: " << std::setprecision(16) << errorNormV << std::endl;
  std::cout << "errorP: " << std::setprecision(16) << errorNormP << std::endl;

  return checkError(
      {errorNormU, errorNormV, errorNormP},
      {5.18103317234e-05, 7.56511221503e-05, 5.25647216673e-05},
      1.e-11);
}
