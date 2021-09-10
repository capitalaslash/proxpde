#include "def.hpp"

#include "ns.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  uint constexpr dim = 2;
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpaceVel_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 2>::RefFE_T,
      typename LagrangeFE<Elem_T, 2>::RecommendedQR,
      dim>;
  using FESpaceU_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 2>::RefFE_T,
      typename LagrangeFE<Elem_T, 2>::RecommendedQR,
      1>;
  using FESpaceP_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 1>::RefFE_T,
      typename LagrangeFE<Elem_T, 2>::RecommendedQR>;
  MilliTimer t;

  ParameterDict config;

  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["origin"] = Vec3{0., 0., 0.};
    config["length"] = Vec3{1., 10., 0.};
    config["nx"] = 4;
    config["ny"] = 8;
    config["dt"] = 0.1;
    config["ntime"] = 50U;
    config["nu"] = 0.1;
  }

  config.validate({"origin", "length", "nx", "ny", "dt", "ntime", "nu"});

  auto const dt = config["dt"].as<double>();
  auto const nu = config["nu"].as<double>();
  NSParameters parMonolithic{dt, nu, "output_ns/monolithic"};
  NSParameters parSplit{dt, nu, "output_ns/split"};

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new (Mesh_T)};
  buildHyperCube(
      *mesh,
      config["origin"].as<Vec3>(),
      config["length"].as<Vec3>(),
      {config["nx"].as<uint>(), config["ny"].as<uint>(), 0});
  t.stop();

  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh, feSpaceVel.dof.size * FESpaceVel_T::dim};
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpacePSplit{*mesh};

  t.start("monolithic bc");
  auto const inlet = [](Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p(0) * p(0))}; };
  // auto const inlet = [] (Vec3 const & p) { return Vec2{0.0, 1.0}; };
  auto const zero2d = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  auto const inletX = [&inlet](Vec3 const & p) { return inlet(p)[0]; };
  auto const inletY = [&inlet](Vec3 const & p) { return inlet(p)[1]; };
  auto const zero = [](Vec3 const &) { return 0.0; };
  auto bcsVel = std::make_tuple(
      BCEss{feSpaceVel, side::BOTTOM},
      BCEss{feSpaceVel, side::RIGHT},
      BCEss{feSpaceVel, side::TOP, uComp},
      BCEss{feSpaceVel, side::LEFT, uComp});
  std::get<0>(bcsVel) << inlet;
  std::get<1>(bcsVel) << zero2d;
  std::get<2>(bcsVel) << zero2d;
  std::get<3>(bcsVel) << zero2d;
  auto const bcsP = std::make_tuple();
  // auto const bcsP = std::make_tuple(
  //       BCEss{feSpaceP, side::TOP, zero});
  t.stop();

  t.start("split bc");
  auto bcsU = std::make_tuple(
      BCEss{feSpaceU, side::BOTTOM},
      BCEss{feSpaceU, side::RIGHT},
      BCEss{feSpaceU, side::TOP, {0}},
      BCEss{feSpaceU, side::LEFT, {0}});
  std::get<0>(bcsU) << inletX;
  std::get<1>(bcsU) << zero;
  std::get<2>(bcsU) << zero;
  std::get<3>(bcsU) << zero;
  auto bcsV =
      std::make_tuple(BCEss{feSpaceU, side::BOTTOM}, BCEss{feSpaceU, side::RIGHT});
  std::get<0>(bcsV) << inletY;
  std::get<1>(bcsV) << zero;
  auto bcTopP = BCEss{feSpacePSplit, side::TOP};
  bcTopP << zero;
  auto const bcsPSplit = std::make_tuple(bcTopP);

  t.stop();

  t.start("monolithic ctor");
  NSSolverMonolithic ns{feSpaceVel, feSpaceP, bcsVel, bcsP, parMonolithic};
  NSSolverSplit2D split{feSpaceU, feSpacePSplit, bcsU, bcsV, bcsPSplit, parSplit};
  t.stop();

  t.start("init");
  ns.init();
  split.init();
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

    time += dt;
    std::cout << separator << "solving timestep " << itime << ", time = " << time
              << std::endl;

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
