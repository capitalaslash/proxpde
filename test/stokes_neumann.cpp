#include "def.hpp"

#include <yaml-cpp/yaml.h>

#include "assembly_bc.hpp"
#include "assembly_lhs.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "ns.hpp"
#include "timer.hpp"
#include "var.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  constexpr uint dim = 2u;
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using QRQuadratic_T = LagrangeFE<Elem_T, 2u>::RecommendedQR;
  using QRConstant_T = LagrangeFE<Elem_T, 2u>::RecommendedQR;
  using RefFEQuadratic_T = LagrangeFE<Elem_T, 2u>::RefFE_T;
  using RefFELinear_T = LagrangeFE<Elem_T, 1u>::RefFE_T;
  using RefFEConstant_T = LagrangeFE<Elem_T, 0u>::RefFE_T;
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
    config["mesh"]["type"] = MeshType::STRUCTURED;
    config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    config["mesh"]["length"] = Vec3{1.0, 1.0, 0.0};
    config["mesh"]["n"] = std::array{10u, 10u, 0u};
    config["mesh"]["flags"] = Bitmask{MeshFlags::BOUNDARY_FACETS};
    config["nu"] = 0.1;
  }
  config.validate({"mesh", "nu"});

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, config["mesh"]);
  t.stop();

  t.start("fespace");
  FESpaceP2V_T feSpaceVel{*mesh};
  FESpaceP2_T feSpaceComponent{*mesh};
  auto const dofU = feSpaceVel.dof.size;
  FESpaceP1_T feSpaceP{*mesh, dim * dofU};
  auto const dofP = feSpaceP.dof.size;
  FESpaceSide_T feSpaceSide{*mesh};
  t.stop();

  t.start("bc");
  auto zero = [](Vec3 const &) { return Vec2::Constant(0.0); };
  auto inlet = [](Vec3 const & p) { return Vec2(0., 0.5 * (1. - p(0) * p(0))); };
  auto const bcsVel = std::vector{
      BCEss{feSpaceVel, side::BOTTOM, inlet},
      BCEss{feSpaceVel, side::RIGHT, zero},
      // BCEss{feSPaceVel, side::TOP, zero, {0}},
      BCEss{feSpaceVel, side::LEFT, zero, {0}},
  };
  auto const pOutlet = 1.0;
  auto pOutletVar = FEVar{"pOutlet", feSpaceSide};
  pOutletVar.data = Vec::Constant(feSpaceSide.dof.size, -pOutlet);
  auto const bcTop = AssemblyBCNormalFE{pOutletVar, side::TOP, feSpaceVel};
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
  LUSolver solver(builder.A);
  Vec const sol = solver.solve(builder.b);
  t.stop();

  t.start("exact");
  Vec exact{sol.size()};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact);
  interpolateAnalyticFunction(
      [nu, pOutlet](Vec3 const & p) { return pOutlet + nu * (1.0 - p(1)); },
      feSpaceP,
      exact);

  Var u{"u"};
  Var v{"v"};
  getComponent(u.data, feSpaceComponent, sol, feSpaceVel, 0u);
  getComponent(v.data, feSpaceComponent, sol, feSpaceVel, 1u);
  Var p{"p"};
  p.data = sol.tail(dofP);

  Var uExact{"ue"};
  Var vExact{"ve"};
  getComponent(uExact.data, feSpaceComponent, exact, feSpaceVel, 0u);
  getComponent(vExact.data, feSpaceComponent, exact, feSpaceVel, 1u);
  Var pExact{"pe"};
  pExact.data = exact.tail(dofP);

  // wall shear stress
  t.start("wssCell");
  FESpaceP0_T feSpaceP0{*mesh};
  // mean cell value
  Var wssCell{"wssCell"};
  computeElemWSS(wssCell.data, feSpaceP0, sol, feSpaceVel, {side::RIGHT}, nu);
  t.stop();
  Eigen::Index wssCellMaxPos;
  auto const wssCellMaxVal = wssCell.data.maxCoeff(&wssCellMaxPos);
  fmt::println("wssCell min: {:.16e}", wssCell.data.minCoeff());
  fmt::println("wssCell max: {:.16e}", wssCellMaxVal);

  t.start("wssFacet");
  Var wssFacet{"wssFacet"};
  computeFEWSS(wssFacet.data, feSpaceP0, sol, feSpaceVel, {side::RIGHT}, nu);
  t.stop();
  fmt::println("wssFacet min: {:.16e}", wssFacet.data.minCoeff());
  fmt::println("wssFacet max: {:.16e}", wssFacet.data.maxCoeff());

  t.start("print");
  IOManager ioComponent{feSpaceComponent, "output_stokes_neumann/vel"};
  ioComponent.print({u, v, uExact, vExact});
  IOManager ioP{feSpaceP, "output_stokes_neumann/p"};
  ioP.print({p, pExact});
  IOManager ioWSS{feSpaceP0, "output_stokes_neumann/wss"};
  ioWSS.print({wssCell, wssFacet});
  t.stop();

  t.print();

  auto const uError = (u.data - uExact.data).norm();
  auto const vError = (v.data - vExact.data).norm();
  auto const pError = (p.data - pExact.data).norm();

  Var distCell{"distCell"};
  auto const bigNum = 1e+20;
  distCell.data = Vec::Ones(feSpaceP0.dof.size) * bigNum;

  for (auto & e: mesh->elementList)
  {
    auto const elemDof = feSpaceP0.dof.getId(e.id);
    for (auto s = 0u; s < Elem_T::numFacets; s++)
    {
      auto const facetId = mesh->elemToFacet[e.id][s];
      if (facetId != idNotSet)
      {
        auto const & facet = mesh->facetList[facetId];
        if (facet.onBoundary())
          distCell.data[elemDof] = std::min(
              distCell.data[elemDof], (e.midpoint() - facet.midpoint()).norm());
      }
    }
    if (std::fabs(distCell.data[elemDof] - bigNum) < 1e-6)
      distCell.data[elemDof] = 0.0;
  }

  return checkError(
      {
          uError,
          vError,
          pError,
          wssCellMaxVal,
          wssFacet.data.maxCoeff(),
      },
      {
          1.2253585539112627e-14,
          1.7077476426602275e-14,
          2.5919844141354946e-14,
          nu * (1.0 - distCell.data[wssCellMaxPos]),
          nu,
      },
      1e-15);
}
