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

#include <yaml-cpp/yaml.h>

static constexpr uint dim = 2;
using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = FEType<Elem_T, dim>::RefFE_T;
using LinearRefFE = FEType<Elem_T,1>::RefFE_T;
using QuadraticQR = FEType<Elem_T, dim>::RecommendedQR;
using FESpaceVel_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR, dim>;
using FESpaceP_T = FESpace<Mesh_T,LinearRefFE,QuadraticQR>;

int main(int argc, char* argv[])
{
  MilliTimer t;

  t.start();
  auto const configFile = (argc > 1) ? argv[1] : "cavitytime.yaml";
  YAML::Node const config = YAML::LoadFile(configFile);
  std::cout << "read config file: " << t << " ms" << std::endl;

  t.start();
  auto const numElemsX = config["nx"].as<uint>();
  auto const numElemsY = config["ny"].as<uint>();

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {{numElemsX, numElemsY, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  BCList bcsVel{feSpaceVel};
  bcsVel.addBC(BCEss{feSpaceVel, side::RIGHT, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::LEFT, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::BOTTOM, zero});
  bcsVel.addBC(BCEss{feSpaceVel, side::TOP, [] (Vec3 const &) {return Vec2(1.0, 0.0);}});
  BCList bcsP{feSpaceP};
  // select the point on the bottom boundary in the middle
  DOFCoordSet pinSet{
      feSpaceP,
      [](Vec3 const & p){return std::fabs(p[0] - 0.5) < 1e-12 && std::fabs(p[1]) < 1e-12;}};
  bcsP.addBC(BCEss{feSpaceP, pinSet.ids, [] (Vec3 const &) {return 0.;}});
  std::cout << "bcs: " << t << " ms" << std::endl;

  // t.start();
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU*dim + dofP;

  auto const mu = config["mu"].as<double>();
  auto const dt = config["timestep"].as<double>();
  Vec velOld{dofU*dim};

  AssemblyTensorStiffness stiffness(mu, feSpaceVel);
  // AssemblyStiffness stiffness(mu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP, {0,1}, 0, dofU*dim);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel, {0,1}, dofU*dim, 0);
  AssemblyMass timeder(1./dt, feSpaceVel);
  AssemblyProjection timeder_rhs(1./dt, velOld, feSpaceVel);
  AssemblyAdvection advection(1.0, velOld, feSpaceVel);
  // we need this in order to properly apply the pinning bc on the pressure
  AssemblyMass dummy(0.0, feSpaceP, {0}, dofU*dim, dofU*dim);

  Var sol{"vel", numDOFs};
  Var fixedSol{"vel", numDOFs};
  auto ic = [](Vec3 const &) {return Vec2(1., 0.);};
  interpolateAnalyticFunction(ic, feSpaceVel, sol.data);
  fixedSol.data = sol.data;

  IOManager ioVel{feSpaceVel, "output_cavitytime/sol_v"};
  ioVel.print({sol});
  IOManager ioP{feSpaceP, "output_cavitytime/sol_p"};
  Var p{"p", sol.data, dofU*dim, dofP};
  ioP.print({p});

  Builder<StorageType::RowMajor> builder{numDOFs};

  Builder<StorageType::RowMajor> fixedBuilder{numDOFs};
  fixedBuilder.buildProblem(timeder, bcsVel);
  fixedBuilder.buildProblem(stiffness, bcsVel);
  fixedBuilder.buildProblem(grad, bcsVel, bcsP);
  fixedBuilder.buildProblem(div, bcsP, bcsVel);
  fixedBuilder.buildProblem(dummy, bcsP);
  fixedBuilder.closeMatrix();
  auto const fixedMat = fixedBuilder.A;
  auto const fixedRhs = fixedBuilder.b;

  IterSolver solver;
  IterSolver fixedSolver;
  auto const ntime = config["numsteps"].as<uint>();
  double time = 0.0;
  for (uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime
              << ", time = " << time << std::endl;

    velOld = sol.data;

    builder.clear();
    builder.buildProblem(timeder_rhs, bcsVel);
    builder.buildProblem(advection, bcsVel);
    builder.buildProblem(timeder, bcsVel);
    builder.buildProblem(stiffness, bcsVel);
    builder.buildProblem(grad, bcsVel, bcsP);
    builder.buildProblem(div, bcsP, bcsVel);
    builder.buildProblem(dummy, bcsP);
    builder.closeMatrix();

    fixedBuilder.clear();
    fixedBuilder.buildProblem(timeder_rhs, bcsVel);
    fixedBuilder.buildProblem(advection, bcsVel);
    fixedBuilder.closeMatrix();
    fixedBuilder.A += fixedMat;
    fixedBuilder.b += fixedRhs;

    // auto const diffMat = builder.A - fixedBuilder.A;
    // std::cout << "diffMat norm: " << diffMat.norm() << std::endl;

    solver.compute(builder.A);
    sol.data = solver.solve(builder.b);
    fixedSolver.compute(fixedBuilder.A);
    fixedSol.data = fixedSolver.solve(fixedBuilder.b);
    auto res = fixedBuilder.A*sol.data - fixedBuilder.b;
    std::cout << "residual norm: " << res.norm() << std::endl;
    auto const solDiffNorm = (sol.data - fixedSol.data).norm();
    std::cout << "solution difference norm: " << solDiffNorm << std::endl;
    if (solDiffNorm > 1.e-11)
    {
      std::cerr << "the 2 solutions differ" << std::endl;
      return 2;
    }

    ioVel.time = time;
    ioVel.iter += 1;
    ioVel.print({sol});
    p.data = sol.data.block(dofU*dim,0,dofP,1);
    ioP.time = time;
    ioP.iter += 1;
    ioP.print({p});
  }

  auto const solNorm = sol.data.norm();
  auto const targetNorm = 4.40937006291;
  std::cout << "solution norm: " << std::setprecision(12) << solNorm << std::endl;
  if (std::fabs(solNorm - targetNorm) > 1.e-10 * targetNorm)
  {
    std::cerr << "the solution norm is not the prescribed value, distance: "
              << std::fabs(solNorm - targetNorm) << std::endl;
    return 1;
  }

  return 0;
}
