#include "def.hpp"

#include "ns.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;

  MilliTimer t;
  ParameterDict config;

  // default configuration
  config["mesh"]["origin"] = Vec3{0., 0., 0.};
  config["mesh"]["length"] = Vec3{1., 10., 0.};
  config["mesh"]["n"] = std::array{4U, 8U, 0U};
  config["ntime"] = 50U;
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
  NSSolverSplit2D split{*mesh, ParameterDict{config["ns"]}};
  t.stop();

  t.start("bc monolithic");
  auto const inlet = [](Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p(0) * p(0))}; };
  // auto const inlet = [] (Vec3 const & p) { return Vec2{0.0, 1.0}; };
  auto const zero2d = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  auto const inletX = [&inlet](Vec3 const & p) { return inlet(p)[0]; };
  auto const inletY = [&inlet](Vec3 const & p) { return inlet(p)[1]; };
  auto const zero = [](Vec3 const &) { return 0.0; };
  auto bcsVel = std::tuple{
      BCEss{ns.feSpaceVel, side::BOTTOM},
      BCEss{ns.feSpaceVel, side::RIGHT},
      BCEss{ns.feSpaceVel, side::LEFT, Comp::u},
      BCEss{ns.feSpaceVel, side::TOP, Comp::u},
  };
  std::get<0>(bcsVel) << inlet;
  std::get<1>(bcsVel) << zero2d;
  std::get<2>(bcsVel) << zero2d;
  std::get<3>(bcsVel) << zero2d;
  auto const bcsP = std::tuple{};
  // auto const bcsP = std::tuple{BCEss{feSpaceP, side::TOP, zero}};
  t.stop();

  t.start("bc split");
  auto bcsU = std::tuple{
      BCEss{split.feSpaceU, side::BOTTOM},
      BCEss{split.feSpaceU, side::RIGHT},
      BCEss{split.feSpaceU, side::LEFT},
      BCEss{split.feSpaceU, side::TOP},
  };
  std::get<0>(bcsU) << inletX;
  std::get<1>(bcsU) << zero;
  std::get<2>(bcsU) << zero;
  std::get<3>(bcsU) << zero;
  auto bcsV = std::tuple{
      BCEss{split.feSpaceU, side::BOTTOM},
      BCEss{split.feSpaceU, side::RIGHT},
  };
  std::get<0>(bcsV) << inletY;
  std::get<1>(bcsV) << zero;
  auto bcTopP = BCEss{split.feSpaceP, side::TOP};
  bcTopP << zero;
  auto const bcsPSplit = std::tuple{bcTopP};
  t.stop();

  t.start("monolithic init");
  ns.init(bcsVel, bcsP);
  t.stop();

  t.start("split init");
  split.init(bcsU, bcsV, bcsPSplit);
  t.stop();

  t.start("ic");
  auto const ic = [](Vec3 const &) { return Vec2{0.0, 1.0}; };
  // auto const ic = [] (Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p(0)*p(0))}; };
  ns.ic(ic);
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
    split.assemblyStepVelStar(bcsU, bcsV);
    t.stop();

    t.start("split solve");
    split.solveVelStar();
    t.stop();

    t.start("split assembly");
    split.assemblyStepP(bcsP);
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

  // std::array<Vec, 2> vel = {Vec{split.feSpaceU.dof.size},
  // Vec{split.feSpaceU.dof.size}}; getComponents(vel, ns.sol.data, ns.feSpaceVel);
  Vec u{split.feSpaceU.dof.size};
  Vec v{split.feSpaceU.dof.size};
  getComponent(u, split.feSpaceU, ns.sol.data, ns.feSpaceVel, 0);
  getComponent(v, split.feSpaceU, ns.sol.data, ns.feSpaceVel, 1);

  Var errorU{"errorU"};
  errorU.data = u - split.u.data;
  Var errorV{"errorV"};
  errorV.data = v - split.v.data;
  Var errorP{"errorP"};
  errorP.data = ns.p.data - split.p.data;

  IOManager ioErrorU{split.feSpaceU, "output_ns/eu"};
  ioErrorU.print({errorU, errorV});
  IOManager ioErrorP{split.feSpaceP, "output_ns/ep"};
  ioErrorP.print({errorP});

  double const errorNormU = errorU.data.norm();
  double const errorNormV = errorV.data.norm();
  double const errorNormP = errorP.data.norm();
  std::cout << "errorU: " << std::setprecision(16) << errorNormU << std::endl;
  std::cout << "errorV: " << std::setprecision(16) << errorNormV << std::endl;
  std::cout << "errorP: " << std::setprecision(16) << errorNormP << std::endl;

  return checkError(
      {errorNormU, errorNormV, errorNormP},
      {5.18337674567542e-05, 7.57375946701665e-05, 5.272719445641985e-05},
      1.e-11);
}
