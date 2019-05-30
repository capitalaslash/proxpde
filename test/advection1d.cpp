#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"
#include "fv.hpp"

#include <yaml-cpp/yaml.h>

#include <iostream>

using Elem_T = Line;
using Mesh_T = Mesh<Elem_T>;
// implicit finite element central
using FESpaceP1_T = FESpace<Mesh_T,
                            FEType<Elem_T, 1>::RefFE_T,
                            FEType<Elem_T, 1>::RecommendedQR>;
// explicit finite volume upwind
using FESpaceP0_T = FESpace<Mesh_T,
                            FEType<Elem_T, 0>::RefFE_T,
                            FEType<Elem_T, 0>::RecommendedQR>;
using FVSolver_T = FVSolver<FESpaceP0_T, LimiterType::UPWIND>;

int main(int argc, char* argv[])
{
  MilliTimer t;

  YAML::Node config;
  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config = YAML::LoadFile("advection1d.yaml");
    // n: 10,
    // dt: 0.1,
    // velocity: 0.5,
    // threshold: 0.37,
    // final_time: 3.0,
  }

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  uint const numElems = config["n"].as<uint>();
  buildHyperCube(*mesh, origin, length, {numElems, 0, 0}, INTERNAL_FACETS | NORMALS);
  // auto const r = 1. / 1.1;
  // auto const starting = (1. - r) / (1. - pow(r, numElems));
  // auto counter = 0;
  // for (auto & p: mesh->pointList)
  // {
  //   p.coord[0] = starting * (1. - pow(r, counter)) / (1. - r);
  //   counter++;
  // }
  t.stop();

  t.start("fespace");
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
  t.stop();

  t.start("bcs");
  auto const one = [](Vec3 const & ){return 1.;};
  BCList bcsP1{feSpaceP1};
  bcsP1.addBC(BCEss{feSpaceP1, side::LEFT, one});
  BCList bcsP0{feSpaceP0};
  bcsP0.addBC(BCEss{feSpaceP0, side::LEFT, one});
  t.stop();

  auto const dt = config["dt"].as<double>();
  // the only 1D velocity that is divergence free is a constant one
  auto const velocity = config["velocity"].as<double>();
  Table<double, 1> vel(feSpaceP1.dof.size, 1);
  vel.block(0, 0, feSpaceP1.dof.size, 1) = Vec::Constant(feSpaceP1.dof.size, velocity);
  double const hinv = numElems;
  std::cout << "cfl = " << velocity * dt * hinv << std::endl;

  Builder builder{feSpaceP1.dof.size};
  LUSolver solver;
  AssemblyAdvection advection(1.0, vel, feSpaceP1);
  AssemblyMass timeder(1./dt, feSpaceP1);
  Vec concP1Old(feSpaceP1.dof.size);
  AssemblyProjection timeder_rhs(1./dt, concP1Old, feSpaceP1);

  auto const threshold = config["threshold"].as<double>();
  scalarFun_T ic = [threshold] (Vec3 const& p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
    if(p(0) < threshold) return 1.;
    return 0.;
  };
  Var concP1{"conc"};
  interpolateAnalyticFunction(ic, feSpaceP1, concP1.data);

  t.start("p0 ic");
  Var concP0{"concP0"};
  // we need to use the highest order available QR to integrate discontinuous functions
  // FESpace<Mesh_T, RefLineP0, GaussQR<Line, 4>> feSpaceIC{*mesh};
  FESpace<Mesh_T, RefLineP0, MiniQR<Line, 20>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, concP0.data);
  t.stop();

  FVSolver_T fv{feSpaceP0, bcsP0};

  auto const ntime = static_cast<uint>(std::nearbyint(config["final_time"].as<double>() / dt));
  double time = 0.0;
  IOManager ioP1{feSpaceP1, "output_advection1d/solP1"};
  ioP1.print({concP1});
  IOManager ioP0{feSpaceP0, "output_advection1d/solP0"};
  ioP0.print({concP0});

  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << ", time = " << time << std::endl;

    // central implicit
    t.start("p1 assemby");
    concP1Old = concP1.data;
    builder.buildProblem(timeder, bcsP1);
    builder.buildProblem(timeder_rhs, bcsP1);
    builder.buildProblem(advection, bcsP1);
    builder.closeMatrix();
    t.stop();

    t.start("p1 solve");
    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    concP1.data = solver.solve(builder.b);
    builder.clear();
    t.stop();

    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;
    // std::cout << "sol:\n" << c.data << std::endl;

    // explicit upwind
    t.start("p0 update");
    fv.update(concP0.data);
    fv.computeFluxes(vel, feSpaceP1);
    fv.advance(concP0.data, dt);
    t.stop();

    // print
    t.start("print");
    ioP1.time = time;
    ioP1.iter += 1;
    ioP1.print({concP1});
    ioP0.time = time;
    ioP0.iter += 1;
    ioP0.print({concP0});
    t.stop();
  }

  t.print();

  Vec oneFieldP1;
  interpolateAnalyticFunction(one, feSpaceP1, oneFieldP1);
  Vec oneFieldP0;
  interpolateAnalyticFunction(one, feSpaceP0, oneFieldP0);

  double errorNormP1 = (concP1.data - oneFieldP1).norm();
  std::cout << "the norm of the P1 error is " << std::setprecision(16) << errorNormP1 << std::endl;
  double errorNormP0 = (concP0.data - oneFieldP0).norm();
  std::cout << "the norm of the P0 error is " << std::setprecision(16) << errorNormP0 << std::endl;
  if (std::fabs(errorNormP1 - 0.01153555695665251) > 1.e-10 ||
      std::fabs(errorNormP0 - 0.0003358552892295136) > 1.e-10)
  {
     std::cerr << "the norm of the error is not the prescribed value" << std::endl;
     return 1;
  }

  return 0;
}
