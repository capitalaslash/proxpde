#include "def.hpp"
#include "ns.hpp"
#include "timer.hpp"


int main(int argc, char* argv[])
{
  MilliTimer t;

  ParameterDict config;

  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["nx"] = 4;
    config["ny"] = 8;
    config["dt"] = 0.1;
    config["ntime"]= 50U;
    config["nu"] = 0.1;
  }

  validateConfig(config, {"nx", "ny", "dt", "ntime", "nu"});

  auto const dt = config["dt"].as<double>();
  auto const nu = config["nu"].as<double>();
  NSParameters parMonolithic{dt, nu, "output_ns/monolithic"};
  NSParameters parSplit{dt, nu, "output_ns/split"};

  t.start("mesh");
  std::unique_ptr<Mesh<Quad>> mesh{new(Mesh<Quad>)};
  buildHyperCube(
        *mesh,
        {0.0, 0.0, 0.0},
        {1.0, 10.0, 0.0},
        {config["nx"].as<uint>(), config["ny"].as<uint>(), 0});
  t.stop();

  t.start("monolithic ctor");
  NSSolverMonolithic ns{*mesh, parMonolithic};
  NSSolverSplit2D split{*mesh, parSplit};
  t.stop();

  t.start("monolithic bc");
  auto const inlet = [] (Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p(0)*p(0))}; };
  // auto const inlet = [] (Vec3 const & p) { return Vec2{0.0, 1.0}; };
  auto const zero2d = [] (Vec3 const & ) { return Vec2{0.0, 0.0}; };
  auto const inletX = [&inlet] (Vec3 const & p) { return inlet(p)[0];};
  auto const inletY = [&inlet] (Vec3 const & p) { return inlet(p)[1];};
  auto const zero = [] (Vec3 const & ) { return 0.0; };
  ns.bcsVel.addBC(BCEss{ns.feSpaceVel, side::BOTTOM, inlet});
  ns.bcsVel.addBC(BCEss{ns.feSpaceVel, side::RIGHT, zero2d});
  ns.bcsVel.addBC(BCEss{ns.feSpaceVel, side::TOP, zero2d, uComp});
  ns.bcsVel.addBC(BCEss{ns.feSpaceVel, side::LEFT, inlet, uComp});
  // ns.bcsP.addBC(BCEss{ns.feSpaceP, side::TOP, [](Vec3 const &){return 0.;}});
  t.stop();

  t.start("split bc");
  split.bcsU.addBC(BCEss{split.feSpaceU, side::BOTTOM, inletX});
  split.bcsU.addBC(BCEss{split.feSpaceU, side::RIGHT, zero});
  split.bcsU.addBC(BCEss{split.feSpaceU, side::TOP, zero});
  split.bcsU.addBC(BCEss{split.feSpaceU, side::LEFT, zero});
  split.bcsV.addBC(BCEss{split.feSpaceU, side::BOTTOM, inletY});
  split.bcsV.addBC(BCEss{split.feSpaceU, side::RIGHT, zero});
  split.bcsP.addBC(BCEss{split.feSpaceP, side::TOP, zero});
  t.stop();

  t.start("init");
  ns.init();
  split.init();
  t.stop();

  t.start("ic");
  auto const ic = [] (Vec3 const & ) { return Vec2{0.0, 1.0}; };
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
  for (uint itime=0; itime<ntime; ++itime)
  {
    MilliTimer stepTimer;
    stepTimer.start();

    time += dt;
    std::cout << separator << "solving timestep " << itime
              << ", time = " << time << std::endl;

    t.start("monolithic assembly");
    ns.assemblyStep();
    t.stop();

    t.start("monolithic solve");
    ns.solve();
    t.stop();

    t.start("split assembly");
    split.assemblyStepVelStar();
    t.stop();

    t.start("split solve");
    split.solveVelStar();
    t.stop();

    t.start("split assembly");
    split.assemblyStepP();
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

  // array<Vec, 2> vel = {Vec{split.feSpaceU.dof.size}, Vec{split.feSpaceU.dof.size}};
  // getComponents(vel, ns.sol.data, ns.feSpaceVel);
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

  std::cout << "errorU: " << std::setprecision(16) << errorU.data.norm() << std::endl;
  std::cout << "errorV: " << std::setprecision(16) << errorV.data.norm() << std::endl;
  std::cout << "errorP: " << std::setprecision(16) << errorP.data.norm() << std::endl;

  if (std::fabs(errorU.data.norm() - 5.188514402817515e-05) > 1.e-12 ||
      std::fabs(errorV.data.norm() - 7.581224214010845e-05) > 1.e-12 ||
      std::fabs(errorP.data.norm() - 5.278919058528884e-05) > 1.e-11 )
  {
    std::cerr << "one of the error norms does not coincide with its expected value." << std::endl;
    return 1;
  }
  return 0;
}
