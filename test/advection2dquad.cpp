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

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          FEType<Elem_T,1>::RecommendedQR>;

static scalarFun_T ic = [] (Vec3 const& p)
{
  // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
  if(p(0) < .4) return 1.;
  return 0.;
};

int main(int argc, char* argv[])
{
  MilliTimer t;
  uint const numPts_x = (argc < 3)? 5 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 5 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start();
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, origin, length, {{numPts_x, numPts_y, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace(*mesh);
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [](Vec3 const &){return 1.;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  double const dt = 0.1;

  Vec vel = Vec::Zero(2*feSpace.dof.size);
  vel.block(0, 0, feSpace.dof.size, 1) = Vec::Constant(feSpace.dof.size, 0.1);
  AssemblyAdvection advection(1.0, vel, feSpace);
  AssemblyMass timeder(1./dt, feSpace);
  Vec cOld(feSpace.dof.size);
  AssemblyProjection timeder_rhs(1./dt, cOld, feSpace);

  Var c{"conc"};
  interpolateAnalyticFunction(ic, feSpace, c.data);
  IOManager io{feSpace, "output_advection2dquad/sol"};

  Builder builder{feSpace.dof.size};
  LUSolver solver;
  uint const ntime = 200;
  double time = 0.0;
  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << std::endl;

    cOld = c.data;

    builder.buildProblem(timeder, bcs);
    builder.buildProblem(timeder_rhs, bcs);
    builder.buildProblem(advection, bcs);
    builder.closeMatrix();
    std::cout << "A:\n" << builder.A << std::endl;
    std::cout << "b:\n" << builder.b << std::endl;

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    c.data = solver.solve(builder.b);
    builder.clear();

    // std::cout << "sol:\n" << c.data << std::endl;

    io.time = time;
    io.iter += 1;
    io.print({c});
  }

  // double norm = error.data.norm();
  // std::cout << "the norm of the error is " << norm << std::endl;
  // if(std::fabs(norm - 2.61664e-11) > 1.e-10)
  // {
  //   std::cerr << "the norm of the error is not the prescribed value" << std::endl;
  //   return 1;
  // }

  return 0;
}
