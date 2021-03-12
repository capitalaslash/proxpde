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

static constexpr uint dim = 2;
using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = LagrangeFE<Elem_T, dim>::RefFE_T;
using LinearRefFE = LagrangeFE<Elem_T,1>::RefFE_T;
using QuadraticQR = LagrangeFE<Elem_T, dim>::RecommendedQR;
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
  FESpaceP_T feSpaceP{*mesh, feSpaceVel.dof.size * dim};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto zero = [] (Vec3 const &) {return Vec2::Constant(0.);};
  auto bcsVel = std::make_tuple(
        BCEss{feSpaceVel, side::RIGHT},
        BCEss{feSpaceVel, side::LEFT},
        BCEss{feSpaceVel, side::BOTTOM},
        BCEss{feSpaceVel, side::TOP});
  std::get<0>(bcsVel) << zero;
  std::get<1>(bcsVel) << zero;
  std::get<2>(bcsVel) << zero;
  std::get<3>(bcsVel) << [] (Vec3 const &) { return Vec2(1.0, 0.0); };
  // select the point on the bottom boundary in the middle
  DOFCoordSet pinSet{
      feSpaceP,
      [](Vec3 const & p){return std::fabs(p[0] - 0.5) < 1e-12 && std::fabs(p[1]) < 1e-12;}
  };
  auto bcPin = BCEss{feSpaceP, pinSet.ids};
  bcPin << [] (Vec3 const &) { return 0.; };
  auto const bcsP = std::tuple{bcPin};
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
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);
  AssemblyScalarMass timeder(1./dt, feSpaceVel);
  AssemblyProjection timeder_rhs(1./dt, velOld, feSpaceVel);
  AssemblyAdvection advection(1.0, velOld, feSpaceVel, feSpaceVel);
  // we need this in order to properly apply the pinning bc on the pressure
  AssemblyDummy dummy{feSpaceP};

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
  fixedBuilder.buildLhs(std::tuple{timeder, stiffness}, bcsVel);
  fixedBuilder.buildCoupling(grad, bcsVel, bcsP);
  fixedBuilder.buildCoupling(div, bcsP, bcsVel);
  fixedBuilder.buildLhs(std::tuple{dummy}, bcsP);
  fixedBuilder.closeMatrix();
  auto const fixedMat = fixedBuilder.A;
  auto const fixedRhs = fixedBuilder.b;

  IterSolver solver;
  IterSolver fixedSolver;
  auto const ntime = config["numsteps"].as<uint>();
  double time = 0.0;
  auto const lhs = std::tuple{advection, timeder, stiffness};
  auto const rhs = std::tuple{timeder_rhs};
  for (uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime
              << ", time = " << time << std::endl;

    velOld = sol.data;

    builder.clear();
    builder.buildRhs(rhs, bcsVel);
    builder.buildLhs(lhs, bcsVel);
    builder.buildCoupling(grad, bcsVel, bcsP);
    builder.buildCoupling(div, bcsP, bcsVel);
    builder.buildLhs(std::tuple{dummy}, bcsP);
    builder.closeMatrix();

    fixedBuilder.clear();
    fixedBuilder.buildRhs(std::tuple{timeder_rhs}, bcsVel);
    fixedBuilder.buildLhs(std::tuple{advection}, bcsVel);
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

    ioVel.print({sol}, time);
    p.data = sol.data.block(dofU*dim,0,dofP,1);
    ioP.print({p}, time);
  }

  auto const solNorm = sol.data.norm();
  std::cout << "solution norm: " << std::setprecision(16) << solNorm << std::endl;
  return checkError({solNorm}, {4.40937006291}, 1.e-10);
}
