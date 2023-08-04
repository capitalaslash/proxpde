#include "def.hpp"

#include "assembler.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"
#include "var.hpp"

int main()
{
  using namespace proxpde;

  ParameterDict config = YAML::LoadFile("steptime.yaml");

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, 2>;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;

  MilliTimer t;
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start();
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh, feSpaceVel.dof.size * FESpaceVel_T::dim};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto zero = [](Vec3 const &) { return Vec2::Constant(0.); };
  auto bcLeft = BCEss{feSpaceVel, side::LEFT};
  bcLeft << [](Vec3 const & p) { return p[1] > .05 ? Vec2(1.0, 0.0) : Vec2(0.0, 0.0); };
  auto const bcBottom = BCEss{feSpaceVel, side::BOTTOM, zero};
  auto const bcTop = BCEss{feSpaceVel, side::TOP, zero};
  auto const bcsVel = std::vector{bcLeft, bcBottom, bcTop};

  auto const bcsP = std::vector<BCEss<FESpaceP_T>>{};
  // DofSet_T pinSet = {1};
  // auto const bcsP =
  //     std::vector{BCEss{feSpaceP, pinSet, [](Vec3 const &) { return 0.; }}};
  std::cout << "bcs: " << t << " ms" << std::endl;

  // t.start();
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU * FESpaceVel_T::dim + dofP;

  auto const mu = config["mu"].as<double>();
  AssemblyTensorStiffness stiffness{mu, feSpaceVel};
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);

  auto const dt = config["timestep"].as<double>();
  AssemblyScalarMass timeder(1. / dt, feSpaceVel);
  Vec velOld{2 * dofU};
  AssemblyProjection timeder_rhs(1. / dt, velOld, feSpaceVel);
  AssemblyAdvection advection(1.0, velOld, feSpaceVel, feSpaceVel);

  Var sol{"sol", numDOFs};
  auto ic = [](Vec3 const &) { return Vec2(1., 0.); };
  interpolateAnalyticFunction(ic, feSpaceVel, sol.data);

  IOManager ioVel{feSpaceVel, "output_steptime/sol_v"};
  ioVel.print({sol});
  IOManager ioP{feSpaceP, "output_steptime/sol_p"};
  Var p{"p", sol.data, 2 * dofU, dofP};
  ioP.print({p});

  Builder builder{numDOFs};
  IterSolver solver;
  auto const ntime = config["numsteps"].as<uint>();
  double time = 0.0;
  for (uint itime = 0; itime < ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime + 1 << ", time = " << time << std::endl;

    velOld = sol.data;

    builder.buildLhs(std::tuple{timeder, advection, stiffness}, bcsVel);
    builder.buildRhs(std::tuple{timeder_rhs}, bcsVel);
    builder.buildCoupling(grad, bcsVel, bcsP);
    builder.buildCoupling(div, bcsP, bcsVel);
    builder.closeMatrix();

    solver.compute(builder.A);
    sol.data = solver.solve(builder.b);
    auto res = builder.A * sol.data - builder.b;
    std::cout << "residual norm: " << res.norm() << std::endl;

    builder.clear();

    ioVel.print({sol}, time);
    p.data = sol.data.tail(dofP);
    ioP.print({p}, time);

    if ((sol.data - dt * velOld).norm() < 1.e-12)
    {
      std::cout << "stationary state has been reached" << std::endl;
      break;
    }
  }

  return 0;
}
