#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"
#include "var.hpp"

#include <yaml-cpp/yaml.h>

static constexpr uint dim = 3;
using Elem_T = Tetrahedron;
using Mesh_T = Mesh<Elem_T>;
using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;
using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, dim>;

template <int n>
static constexpr std::vector<uint> vecIota()
{
  std::vector<uint> comp(n);
  std::iota(comp.begin(), comp.end(), 0);
  return comp;
}
static const auto dimIota = vecIota<dim>();

static auto zero = [](Vec3 const &) { return Vec2::Constant(0.); };
static std::map<std::string, std::pair<Fun<dim, 3>, std::vector<uint>>> bcTypes = {
    {"outlet", {zero, {0}}}, // should clamp directions tangent to the face
    {"wall", {zero, dimIota}},
    {"sym", {zero, {0}}}, // should clamp direction normal to the face
};
int main(int argc, char * argv[])
{
  MilliTimer t;
  auto const configFile = (argc > 1) ? argv[1] : "config.yaml";
  YAML::Node const config = YAML::LoadFile(configFile);

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readGMSH(*mesh, config["meshFile"].as<std::string>());

  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh};

  FVec<dim> const vInlet =
      FVec<dim>(config["v_inlet"].as<std::vector<double>>().data());
  auto inlet = [&vInlet](Vec3 const &) { return vInlet; };
  bcTypes["inlet"] = std::pair(inlet, dimIota);
  BCList bcsVel{feSpaceVel};
  // bcsVel.addEssentialBC(1, inlet);     // inlet
  // bcsVel.addEssentialBC(2, zero, {0}); // outlet
  // bcsVel.addEssentialBC(3, zero);      // heat
  // bcsVel.addEssentialBC(4, zero, {0}); // sym
  // bcsVel.addEssentialBC(5, zero);      // wall
  for (auto const & [boundaryId, boundaryType]:
       config["bcs"].as<std::map<marker_T, std::string>>())
  {
    std::cout << boundaryId << " -> " << boundaryType << std::endl;
    bcsVel.addEssentialBC(
        boundaryId, bcTypes[boundaryType].first, bcTypes[boundaryType].second);
  }
  BCList bcsP{feSpaceP};

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;
  uint const numDOFs = dofU * dim + dofP;

  Vec velOld{dofU * dim};
  double const dt = config["dt"].as<double>();
  double const mu = config["mu"].as<double>();
  AssemblyMass timeder{1. / dt, feSpaceVel};
  AssemblyAdvection advection{1.0, velOld, feSpaceVel};
  AssemblyTensorStiffness stiffness{mu, feSpaceVel};
  AssemblyGrad grad{-1.0, feSpaceVel, feSpaceP, {0, 1}, 0, dofU * dim};
  AssemblyDiv div{-1.0, feSpaceP, feSpaceVel, {0, 1}, dofU * dim, 0};
  AssemblyProjection timederRhs(1. / dt, velOld, feSpaceVel);

  Var sol{"vel", numDOFs};
  Var p{"p", sol.data, dofU * dim, dofP};

  GMRESSolver solver;

  IOManager ioVel{feSpaceVel, "output/sol_v"};
  ioVel.print({sol});
  IOManager ioP{feSpaceP, "output/sol_p"};
  ioP.print({p});

  Builder builder{numDOFs};
  builder.buildProblem(timeder, bcsVel);
  builder.buildProblem(stiffness, bcsVel);
  builder.buildProblem(grad, bcsVel, bcsP);
  builder.buildProblem(div, bcsP, bcsVel);
  builder.closeMatrix();

  Mat fixedMat = builder.A;
  Vec fixedRhs = builder.b;

  double time = 0.0;
  uint const ntime = config["final_time"].as<double>() / dt;
  uint const print_step = config["print_step"].as<uint>();
  for (uint itime = 0; itime < ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime + 1 << ", time = " << time << std::endl;

    velOld = sol.data;

    builder.clear();
    builder.buildProblem(advection, bcsVel);
    builder.buildProblem(timederRhs, bcsVel);
    builder.closeMatrix();
    builder.A += fixedMat;
    builder.b += fixedRhs;
    // std::cout << "rhs norm: " << builder.b.norm() << std::endl;

    solver.compute(builder.A);
    sol.data = solver.solve(builder.b);
    auto const res = builder.A * sol.data - builder.b;
    std::cout << "residual norm: " << res.norm() << std::endl;

    if (itime % print_step == 0)
    {
      ioVel.iter = itime + 1;
      ioVel.time = time;
      ioVel.print({sol});
      p.data = sol.data.tail(dofP);
      ioP.iter = itime + 1;
      ioP.time = time;
      ioP.print({p});
    }
  }

  return 0;
}
