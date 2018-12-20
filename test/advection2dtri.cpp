#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

#include <iostream>

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          FEType<Elem_T,1>::RecommendedQR>;

static scalarFun_T ic = [] (Vec3 const& /*p*/)
{
  // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
  // if(p(0) < .4) return 1.;
  return 0.;
};

template <typename Mesh>
double computeMaxCFL(Mesh const & mesh, Vec2 const & vel, double const dt)
{
  double cfl = 0.;
  for (auto const & e: mesh.elementList)
  {
    cfl = std::max(cfl, vel.norm() * dt / e.h_min());
  }
  return cfl;
}

int main()
{
  MilliTimer t;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start();
  // MeshBuilder<Elem_T> meshBuilder;
  // uint const numPts_x = (argc < 3)? 5 : std::stoi(argv[1]);
  // uint const numPts_y = (argc < 3)? 5 : std::stoi(argv[2]);
  // meshBuilder.build(*mesh, origin, length, {{numPts_x, numPts_y, 0}});
  readGMSH(*mesh, "square_uns.msh");
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [](Vec3 const &){return 1.;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  auto const velocity = Vec2(0.1, 0.0);
  double const dt = 0.4;
  auto const cfl = computeMaxCFL(*mesh, velocity, dt);
  std::cout << "max cfl = " << cfl << std::endl;

  auto const & size = feSpace.dof.size;
  Vec vel = Vec::Zero(2*size);
  vel.block(   0, 0, size, 1) = Vec::Constant(size, velocity[0]);
  vel.block(size, 0, size, 1) = Vec::Constant(size, velocity[1]);
  AssemblyAdvection advection(1.0, vel, feSpace);
  AssemblyMass timeder(1./dt, feSpace);
  Vec cOld(size);
  AssemblyProjection timeder_rhs(1./dt, cOld, feSpace);

  Var c{"conc"};
  interpolateAnalyticFunction(ic, feSpace, c.data);
  IOManager io{feSpace, "output_advection2dtri/sol"};

  Builder builder{size};
  LUSolver solver;

  uint const ntime = 50;
  double time = 0.0;
  io.time = time;
  io.iter = 0;
  io.print({c});

  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << std::endl;

    cOld = c.data;

    builder.buildProblem(timeder, bcs);
    builder.buildProblem(timeder_rhs, bcs);
    builder.buildProblem(advection, bcs);
    builder.closeMatrix();
    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    c.data = solver.solve(builder.b);
    builder.clear();

    // std::cout << "sol:\n" << c.data << std::endl;

    io.time = time;
    io.iter += 1;
    io.print({c});
  }

  double norm = c.data.norm();
  std::cout << "the norm of the solution is " << std::setprecision(12) << norm << std::endl;
  if(std::fabs(norm - 11.9149220338) > 1.e-10)
  {
    std::cerr << "the norm of the solution is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
