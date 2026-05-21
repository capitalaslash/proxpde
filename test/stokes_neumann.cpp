// proxpde configuration - always first
#include "def.hpp"

// external libs
#include <yaml-cpp/yaml.h>

// proxpde
#include "assembly_bc.hpp"
#include "assembly_lhs.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"
#include "ns.hpp"
#include "timer.hpp"
#include "var.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  constexpr uint dim = Elem_T::dim;
  using Mesh_T = Mesh<Elem_T>;
  using RefFEQuadratic_T = LagrangeFE<Elem_T, 2u>::RefFE_T;
  using RefFELinear_T = LagrangeFE<Elem_T, 1u>::RefFE_T;
  using RefFEConstant_T = LagrangeFE<Elem_T, 0u>::RefFE_T;
  using QRQuadratic_T = LagrangeFE<Elem_T, 2u>::RecommendedQR;
  using QRConstant_T = LagrangeFE<Elem_T, 2u>::RecommendedQR;
  using FESpaceP2V_T = FESpace<Mesh_T, RefFEQuadratic_T, QRQuadratic_T, dim>;
  using FESpaceP2_T = FESpace<Mesh_T, RefFEQuadratic_T, QRQuadratic_T>;
  using FESpaceP1_T = FESpace<Mesh_T, RefFELinear_T, QRQuadratic_T>;
  using FESpaceP0_T = FESpace<Mesh_T, RefFEConstant_T, QRConstant_T>;
  using QRQuadraticSide_T = SideGaussQR<Elem_T, SideQR<QRQuadratic_T>::type::numPts>;
  using FESpaceSide_T = FESpace<Mesh_T, RefFEQuadratic_T, QRQuadraticSide_T>;

  MilliTimer t;

  ParameterDict config;
  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    // config["mesh"]["type"] = MeshType::STRUCTURED;
    // config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    // config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
    // config["mesh"]["n"] = std::array{10u, 10u, 0u};
    // config["mesh"]["flags"] = Bitmask{MeshFlags::BOUNDARY_FACETS};

    config["mesh"]["type"] = MeshType::GMSH;
    auto const pathExe = std::filesystem::path{argv[0]}.parent_path();
    config["mesh"]["filename"] = (pathExe / "square_q.msh").string();

    config["nu"] = 0.1;
  }
  config.validate({"mesh", "nu"});

  t.start("mesh");
  std::unique_ptr<Mesh_T> meshCoarse{new Mesh_T};
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*meshCoarse, config["mesh"]);
  uniformRefine(*meshCoarse, *mesh);

  // std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // readMesh(*mesh, config["mesh"]);
  t.stop();

  t.start("fespace");
  auto feSpaceVel = FESpaceP2V_T{*mesh};
  auto feSpaceComponent = FESpaceP2_T{*mesh};
  auto const dofU = feSpaceVel.dof.size;
  auto feSpaceP = FESpaceP1_T{*mesh, dim * dofU};
  auto const dofP = feSpaceP.dof.size;
  FESpaceSide_T feSpaceSide{*mesh};
  t.stop();

  t.start("bc");
  auto const zero = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  // fully developed
  auto const inlet = [](Vec3 const & p) { return Vec2{0.0, 0.5 * (1. - p[0] * p[0])}; };
  auto const bcsVel = std::vector{
      BCEss{feSpaceVel, side::RIGHT, zero},
      // BCEss{feSpaceVel, side::TOP, zero, {0}},
      BCEss{feSpaceVel, side::LEFT, zero, {0}},
      BCEss{feSpaceVel, side::BOTTOM, inlet},
  };
  auto const pOutlet = 1.0;
  auto pOutletVar = FEVar{"pOutlet", feSpaceSide};
  pOutletVar.data = Vec::Constant(feSpaceSide.dof.size, -pOutlet);
  auto const bcTop = AssemblyBCNormalFE{pOutletVar, side::TOP, feSpaceVel, Comp::y};
  t.stop();

  t.start("build");
  auto const nu = config["nu"].as<double>();
  auto const diffusion = AssemblyTensorStiffness{nu, feSpaceVel};
  auto const grad = AssemblyGrad{-1.0, feSpaceVel, feSpaceP};
  auto const div = AssemblyDiv{-1.0, feSpaceP, feSpaceVel};

  Builder builder{dofU * dim + dofP};
  builder.buildLhs(std::tuple{diffusion}, bcsVel);
  builder.buildCoupling(grad, bcsVel, {});
  builder.buildCoupling(div, {}, bcsVel);
  builder.closeMatrix();
  builder.buildRhs(std::tuple{bcTop}, bcsVel);
  t.stop();

  t.start("solve");
  auto solver = LUSolver{builder.A};
  auto const sol = solver.solve(builder.b);
  auto vel = Var{"vel"};
  vel.data = sol.head(dim * dofU);
  auto u = Var{"u"};
  auto v = Var{"v"};
  getComponent(u.data, feSpaceComponent, vel.data, feSpaceVel, 0u);
  getComponent(v.data, feSpaceComponent, vel.data, feSpaceVel, 1u);
  Var p{"p"};
  p.data = sol.tail(dofP);
  t.stop();

  t.start("exact");
  auto exact = Vec{sol.size()};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact);
  interpolateAnalyticFunction(
      [nu, pOutlet](Vec3 const & p) { return pOutlet + nu * (1.0 - p(1)); },
      feSpaceP,
      exact);
  auto velExact = Var{"velExact"};
  velExact.data = exact.head(dim * dofU);
  auto uExact = Var{"ue"};
  auto vExact = Var{"ve"};
  getComponent(uExact.data, feSpaceComponent, velExact.data, feSpaceVel, 0u);
  getComponent(vExact.data, feSpaceComponent, velExact.data, feSpaceVel, 1u);
  auto pExact = Var{"pe"};
  pExact.data = exact.tail(dofP);

  // wall shear stress
  t.start("wssCell");
  FESpaceP0_T feSpaceP0{*mesh};
  // mean cell value
  Var wssCell{"wssCell"};
  computeElemWSS(wssCell.data, feSpaceP0, vel.data, feSpaceVel, {side::RIGHT}, nu);
  t.stop();

  Eigen::Index wssCellMaxPos;
  auto const wssCellMaxVal = wssCell.data.maxCoeff(&wssCellMaxPos);
  fmt::println("wssCell min: {:.16e}", wssCell.data.minCoeff());
  fmt::println("wssCell max: {:.16e}", wssCellMaxVal);

  t.start("wssFacet");
  Var wssFacet{"wssFacet"};
  computeFEWSS(wssFacet.data, feSpaceP0, sol, feSpaceVel, {side::RIGHT}, nu);
  Var dist{"distance"};
  computeDistanceCell(dist.data, feSpaceP0);
  t.stop();

  auto const wssFacetMaxVal = wssFacet.data.maxCoeff();
  fmt::println("wssFacet min: {:.16e}", wssFacet.data.minCoeff());
  fmt::println("wssFacet max: {:.16e}", wssFacetMaxVal);

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_stokes_neumann/vel"};
  ioVel.print({vel, velExact});
  IOManager ioP{feSpaceP, "output_stokes_neumann/p"};
  ioP.print({p, pExact});
  IOManager ioComponent{feSpaceComponent, "output_stokes_neumann/vel"};
  ioComponent.print({u, v, uExact, vExact});
  IOManager ioWSS{feSpaceP0, "output_stokes_neumann/wss"};
  ioWSS.print({wssCell, wssFacet, dist});
  t.stop();

  t.print();

  auto const uError = (u.data - uExact.data).norm();
  auto const vError = (v.data - vExact.data).norm();
  auto const pError = (p.data - pExact.data).norm();

  auto const wssCellExact = nu * (1.0 - dist.data[wssCellMaxPos]);
  auto const wssFacetExact = nu;

  bool const res =
      checkError(
          {uError, vError, pError},
          {6.8227762276597068e-15, 1.3084695251378470e-14, 2.1939198064224178e-14}) ||
      checkError({wssCellMaxVal}, {wssCellExact}, 3e-04) ||
      checkError({wssFacetMaxVal}, {wssFacetExact});

  return res;
}
