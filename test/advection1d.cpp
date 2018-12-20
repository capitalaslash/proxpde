#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

#include <yaml-cpp/yaml.h>

#include <iostream>

using Elem_T = Line;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          FEType<Elem_T,1>::RecommendedQR>;

int main(int argc, char* argv[])
{
  MilliTimer t;

  auto const configFile = (argc > 1) ? argv[1] : "advection1d.yaml";

  auto const config = YAML::LoadFile(configFile);
  uint const numPts = config["n"].as<uint>() + 1;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start();
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, origin, length, {{numPts, 0, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace(*mesh);
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [](Vec3 const &){return 1.;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  double const dt = config["dt"].as<double>();

  Vec vel = Vec::Constant(feSpace.dof.size, config["velocity"].as<double>());
  AssemblyAdvection advection(1.0, vel, feSpace);
  AssemblyMass timeder(1./dt, feSpace);
  Vec cOld(feSpace.dof.size);
  AssemblyProjection timeder_rhs(1./dt, cOld, feSpace);

  Var c{"conc"};
  double const threshold = config["threshold"].as<double>();
  scalarFun_T ic = [threshold] (Vec3 const& p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
    if(p(0) < threshold) return 1.;
    return 0.;
  };
  interpolateAnalyticFunction(ic, feSpace, c.data);
  IOManager io{feSpace, "output_advection1d/sol"};

  Builder builder{feSpace.dof.size};
  LUSolver solver;
  uint const ntime = config["final_time"].as<double>() / dt;
  double time = 0.0;
  io.time = time;
  io.print({c});
  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << ", time = " << time << std::endl;

    cOld = c.data;

    builder.buildProblem(timeder, bcs);
    builder.buildProblem(timeder_rhs, bcs);
    builder.buildProblem(advection, bcs);
    builder.closeMatrix();

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    c.data = solver.solve(builder.b);
    builder.clear();

    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;
    // std::cout << "sol:\n" << c.data << std::endl;

    io.time = time;
    io.iter += 1;
    io.print({c});
  }

  double norm = c.data.norm();
  std::cout << "the norm of the solution is " << std::setprecision(12) << norm << std::endl;
  if(std::fabs(norm - 3.31662383156) > 1.e-10)
  {
     std::cerr << "the norm of the solution is not the prescribed value" << std::endl;
     return 1;
  }

  return 0;
}
