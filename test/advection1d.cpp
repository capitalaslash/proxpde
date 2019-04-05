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

  auto const configFile = (argc > 1) ? argv[1] : "advection1d.yaml";
  auto const config = YAML::LoadFile(configFile);

  t.start();
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  uint const numElems = config["n"].as<uint>();
  buildHyperCube(*mesh, origin, length, {numElems, 0, 0}, INTERNAL_FACETS | NORMALS);
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto const leftBC = [](Vec3 const &){return 1.;};
  BCList bcsP1{feSpaceP1};
  bcsP1.addBC(BCEss{feSpaceP1, side::LEFT, leftBC});
  BCList bcsP0{feSpaceP0};
  bcsP0.addBC(BCEss{feSpaceP0, side::LEFT, leftBC});
  std::cout << "bcs: " << t << " ms" << std::endl;

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
  Var concP0{"concP0"};
  interpolateAnalyticFunction(ic, feSpaceP0, concP0.data);

  FVSolver_T fv{feSpaceP0, bcsP0};

  uint const ntime = config["final_time"].as<double>() / dt;
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
    concP1Old = concP1.data;

    builder.buildProblem(timeder, bcsP1);
    builder.buildProblem(timeder_rhs, bcsP1);
    builder.buildProblem(advection, bcsP1);
    builder.closeMatrix();

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    concP1.data = solver.solve(builder.b);
    builder.clear();

    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;
    // std::cout << "sol:\n" << c.data << std::endl;

    // explicit upwind
    fv.update(concP0.data);
    fv.computeFluxes(vel, feSpaceP1);
    fv.advance(concP0.data, dt);

    // print
    ioP1.time = time;
    ioP1.iter += 1;
    ioP1.print({concP1});
    ioP0.time = time;
    ioP0.iter += 1;
    ioP0.print({concP0});
  }

  double normP1 = concP1.data.norm();
  double normP0 = concP0.data.norm();
  std::cout << "the norm of the P1 solution is " << std::setprecision(12) << normP1 << std::endl;
  std::cout << "the norm of the P0 solution is " << std::setprecision(12) << normP0 << std::endl;
  if(std::fabs(normP1 - 3.32129054498) > 1.e-10 || std::fabs(normP0 - 3.16198922903) > 1.e-10)
  {
     std::cerr << "the norm of the solution is not the prescribed value" << std::endl;
     return 1;
  }

  return 0;
}
