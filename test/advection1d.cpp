#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "fv.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  // implicit finite element central
  using FESpaceP1_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;
  // explicit finite volume upwind
  using FESpaceP0_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 0>::RefFE_T,
      LagrangeFE<Elem_T, 0>::RecommendedQR>;
  // velocity field
  using FESpaceVel_T = FESpaceP1_T;

  MilliTimer t;

  ParameterDict config;
  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    config["mesh"]["length"] = Vec3{1.0, 0.0, 0.0};
    config["mesh"]["n"] = std::array{10U, 0U, 0U};
    config["mesh"]["flags"] =
        MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS | MeshFlags::FACET_PTRS;
    config["dt"] = 0.1;
    config["final_time"] = 3.0;
    config["velocity"] = 0.5;
    config["threshold"] = 0.37;
  }
  config.validate({"mesh", "dt", "velocity", "threshold", "final_time"});

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, config["mesh"]);
  auto const numElems = config["mesh"]["n"].as<std::array<uint, 3>>()[0];
  // auto const r = 1. / 1.1;
  // auto const starting = (1. - r) / (1. - cepow(r, numElems));
  // auto counter = 0;
  // for (auto & p: mesh->pointList)
  // {
  //   p.coord[0] = starting * (1. - cepow(r, counter)) / (1. - r);
  //   counter++;
  // }
  t.stop();

  t.start("fespace");
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
  t.stop();

  t.start("bcs");
  auto const one = [](Vec3 const &) { return 1.; };
  auto bcLeftP1 = BCEss{feSpaceP1, side::LEFT};
  bcLeftP1 << one;
  auto const bcsP1 = std::vector{bcLeftP1};
  auto bcLeftP0 = BCEss{feSpaceP0, side::LEFT};
  bcLeftP0 << one;
  auto const bcsP0 = std::vector{bcLeftP0};
  t.stop();

  auto const dt = config["dt"].as<double>();
  // the only 1D velocity that is divergence free is a constant one
  auto const velocity = config["velocity"].as<double>();
  FESpaceVel_T feSpaceVel{*mesh};
  FEVar vel{"velocity", feSpaceVel};
  vel << velocity;
  double const cfl = velocity * dt * numElems;
  fmt::print("cfl = {}\n", cfl);
  assert(cfl < 1. - 1.e-8);

  Builder builder{feSpaceP1.dof.size};
  LUSolver solver;
  AssemblyAdvection advection(1.0, vel.data, feSpaceVel, feSpaceP1);
  AssemblyScalarMass timeDer(1. / dt, feSpaceP1);
  Vec concP1Old(feSpaceP1.dof.size);
  AssemblyProjection timeDerRhs(1. / dt, concP1Old, feSpaceP1);

  FEVar concP1{"conc", feSpaceP1};
  auto const threshold = config["threshold"].as<double>();
  scalarFun_T ic = [threshold](Vec3 const & p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
    if (p(0) < threshold)
      return 1.;
    return 0.;
  };
  concP1 << ic;

  t.start("p0 ic");
  FEVar concP0{"concP0", feSpaceP0};
  // we need to use the highest order available QR to integrate discontinuous functions
  // FESpace<Mesh_T, RefLineP0, GaussQR<Line, 4>> feSpaceIC{*mesh};
  FESpace<Mesh_T, RefLineP0, MiniQR<Line, 20>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, concP0.data);
  t.stop();

  FVSolver fv{feSpaceP0, bcsP0, UpwindLimiter{}};

  IOManager ioP1{feSpaceP1, "output_advection1d/solP1"};
  ioP1.print({concP1});
  IOManager ioP0{feSpaceP0, "output_advection1d/solP0"};
  ioP0.print({concP0});

  auto const lhs = std::tuple{timeDer, advection};
  auto const rhs = std::tuple{timeDerRhs};

  auto const ntime =
      static_cast<uint>(std::nearbyint(config["final_time"].as<double>() / dt));
  double time = 0.0;
  for (uint itime = 0; itime < ntime; itime++)
  {
    time += dt;
    fmt::print("solving timestep {:4d}, time = {:.6e}\n", itime, time);

    // central implicit
    t.start("p1 assemby");
    concP1Old = concP1.data;
    builder.buildLhs(lhs, bcsP1);
    builder.buildRhs(rhs, bcsP1);
    builder.closeMatrix();
    t.stop();

    t.start("p1 solve");
    solver.compute(builder.A);
    concP1.data = solver.solve(builder.b);
    builder.clear();
    t.stop();

    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;
    // std::cout << "sol:\n" << c.data << std::endl;

    // explicit upwind
    t.start("p0 update");
    fv.update(concP0.data);
    fv.computeFluxes(vel);
    fv.advance(concP0.data, dt);
    t.stop();

    // print
    t.start("print");
    ioP1.print({concP1}, time);
    ioP0.print({concP0}, time);
    t.stop();
  }

  t.print();

  Vec oneFieldP1;
  interpolateAnalyticFunction(one, feSpaceP1, oneFieldP1);
  Vec oneFieldP0;
  interpolateAnalyticFunction(one, feSpaceP0, oneFieldP0);

  double const errorNormP1 = (concP1.data - oneFieldP1).norm();
  fmt::print("the norm of the P1 error is {:.16e}\n", errorNormP1);
  double const errorNormP0 = (concP0.data - oneFieldP0).norm();
  fmt::print("the norm of the P0 error is {:.16e}\n", errorNormP0);
  return checkError({errorNormP1, errorNormP0}, {1.15355569566e-02, 3.35855289229e-04});
}
