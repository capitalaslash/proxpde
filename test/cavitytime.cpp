#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "assembler.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

// #include <yaml-cpp/yaml.h>

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T,2>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T,2>::RecommendedQR;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR,2>;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;

int main(int argc, char* argv[])
{
  // YAML::Node config = YAML::LoadFile("cavitytime.yaml");

  MilliTimer t;
  uint const numPts_x = 3; // config["nx"].as<uint>()+1;
  uint const numPts_y = 3; // config["ny"].as<uint>()+1;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::shared_ptr<Mesh_T> meshPtr{new Mesh_T};

  t.start();
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(meshPtr, origin, length, {{numPts_x, numPts_y, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpaceVel_T feSpaceVel{meshPtr};
  FESpaceP_T feSpaceP{meshPtr};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  BCList bcsVel{feSpaceVel};
  bcsVel.addEssentialBC(side::RIGHT, zero);
  bcsVel.addEssentialBC(side::LEFT, zero);
  bcsVel.addEssentialBC(side::BOTTOM, zero);
  bcsVel.addEssentialBC(side::TOP, [] (Vec3 const &) {return Vec2(1.0, 0.0);});
  BCList bcsP{feSpaceP};
  // DofSet_T pinSet = {1};
  // bcsP.addEssentialBC(pinSet, [] (Vec3 const &) {return 0.;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  // t.start();
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU*FESpaceVel_T::dim + dofP;

  double const mu = 0.1; // config["mu"].as<double>();
  // AssemblyTensorStiffness stiffness(mu, feSpaceVel);
  AssemblyStiffness stiffness(mu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP, {0,1}, 0, 2*dofU);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel, {0,1}, 2*dofU, 0);

  double const dt = 0.01; // config["timestep"].as<double>();
  AssemblyMass timeder(1./dt, feSpaceVel);
  Vec velOld{2*dofU};
  AssemblyProjection timeder_rhs(1./dt, velOld, feSpaceVel);
  AssemblyAdvection advection(1.0, velOld, feSpaceVel);

  Var sol{"vel", numDOFs};
  auto ic = [](Vec3 const &) {return Vec2(1., 0.);};
  interpolateAnalyticFunction(ic, feSpaceVel, sol.data);

  IOManager ioVel{feSpaceVel, "output_cavitytime/sol_v"};
  ioVel.print({sol});
  IOManager ioP{feSpaceP, "output_cavitytime/sol_p"};
  Var p{"p", sol.data, 2*dofU, dofP};
  ioP.print({p});

  Builder builder{numDOFs};
  Eigen::UmfPackLU<Mat> solver;
  // GMRESSolver solver;
  uint const ntime = 100; // config["numsteps"].as<uint>();
  double time = 0.0;
  for (uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime
              << ", time = " << time << std::endl;

    velOld = sol.data;

    builder.buildProblem(timeder, bcsVel);
    builder.buildProblem(timeder_rhs, bcsVel);
    builder.buildProblem(advection, bcsVel);
    builder.buildProblem(stiffness, bcsVel);
    builder.buildProblem(grad, bcsVel, bcsP);
    builder.buildProblem(div, bcsP, bcsVel);
    builder.closeMatrix();

    solver.compute(builder.A);
    sol.data = solver.solve(builder.b);
    auto res = builder.A*sol.data-builder.b;
    std::cout << "residual norm: " << res.norm() << std::endl;

    builder.clear();

    ioVel.time = time;
    ioVel.iter += 1;
    ioVel.print({sol});
    p.data = sol.data.block(2*dofU,0,dofP,1);
    ioP.time = time;
    ioP.iter += 1;
    ioP.print({p});
  }

  return 0;
}
