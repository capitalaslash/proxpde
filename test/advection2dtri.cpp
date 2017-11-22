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

scalarFun_T ic = [] (Vec3 const& p)
{
  // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
  if(p(0) < .4) return 1.;
  return 0.;
};

scalarFun_T rhs = [] (Vec3 const& p)
{
  return M_PI*std::sin(M_PI*p(0));
};
scalarFun_T exact_sol = [] (Vec3 const& p)
{
  return std::sin(M_PI*p(0))/M_PI + p(0);
};

int main(int argc, char* argv[])
{
  MilliTimer t;
  uint const numPts_x = (argc < 3)? 5 : std::stoi(argv[1]);
  uint const numPts_y = (argc < 3)? 5 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr(new Mesh_T);

  t.start();
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace(meshPtr);
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList<FESpace_T> bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [](Vec3 const &){return 1.;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  double const dt = 0.1;

  Field3 vel = Field3::Zero(feSpace.dof.totalNum, 3);
  vel.col(0) = Vec::Constant(feSpace.dof.totalNum, 0.1);
  AssemblyAdvection<FESpace_T> advection(vel, feSpace);
  AssemblyMass<FESpace_T> timeder(1./dt, feSpace);
  Vec c_old(feSpace.dof.totalNum);
  AssemblyVecRhs<FESpace_T> timeder_rhs(c_old, feSpace);

  Var c{"conc"};
  c.data = Vec::Zero(feSpace.dof.totalNum);
  interpolateAnalyticFunction(ic, feSpace, c.data);
  IOManager<FESpace_T> io{feSpace, "sol_advection2dtri.xmf", 0.0};

  Builder builder{feSpace.dof.totalNum};
  LUSolver solver;
  uint const ntime = 200;
  for(uint itime=0; itime<ntime; itime++)
  {
    std::cout << "solving timestep " << itime << std::endl;

    c_old = c.data / dt;

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

    io.fileName = "output/sol_advection2dtri_" + std::to_string(itime) + ".xmf";
    io.time = (itime+1) * dt;
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
