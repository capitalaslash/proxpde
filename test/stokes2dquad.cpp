#include <algorithm>     // for min
#include <array>         // for array
#include <cmath>         // for sqrt, fabs
#include <cstdlib>       // for uint
#include <filesystem>    // for path
#include <memory>        // for unique_ptr, cons...
#include <optional>      // for optional
#include <string>        // for basic_string
#include <tuple>         // for tuple, get
#include <type_traits>   // for is_base_of_v
#include <unordered_set> // for unordered_set
#include <vector>        // for vector

#include <Eigen/Core>      // for DenseBase::row
#include <fmt/base.h>      // for println
#include <yaml-cpp/yaml.h> // for Node, LoadFile

#include <proxpde/assembler.hpp>         // for make_assembler
#include <proxpde/assembly_coupling.hpp> // for AssemblyDiv, Ass...
#include <proxpde/assembly_lhs.hpp>      // for AssemblyTensorSt...
#include <proxpde/bc.hpp>                // for BCEss
#include <proxpde/builder.hpp>           // for Builder
#include <proxpde/def.hpp>               // for ParameterDict, Vec3
#include <proxpde/fe.hpp>                // for LagrangeFE
#include <proxpde/fespace.hpp>           // for FESpace, getComp...
#include <proxpde/geo.hpp>               // for Quad
#include <proxpde/iomanager.hpp>         // for IOManager
#include <proxpde/mesh.hpp>              // for Mesh, side, Mesh...
#include <proxpde/ns.hpp>                // for computeDistanceCell
#include <proxpde/reffe.hpp>             // for RefQuadQ1, RefQu...
#include <proxpde/sparse_matrix.hpp>     // for LUSolver
#include <proxpde/timer.hpp>             // for MilliTimer
#include <proxpde/var.hpp>               // for Var

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
  using FESpaceVel_T = FESpace<Mesh_T, RefFEQuadratic_T, QRQuadratic_T, dim>;
  using FESpaceComponent_T = FESpace<Mesh_T, RefFEQuadratic_T, QRQuadratic_T>;
  using FESpaceP_T = FESpace<Mesh_T, RefFELinear_T, QRQuadratic_T>;
  using FESpaceP0_T = FESpace<Mesh_T, RefFEConstant_T, QRConstant_T>;

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
  auto feSpaceVel = FESpaceVel_T{*mesh};
  auto feSpaceComponent = FESpaceComponent_T{*mesh};
  auto const dofU = feSpaceVel.dof.size;
  auto feSpaceP = FESpaceP_T{*mesh, dim * dofU};
  auto const dofP = feSpaceP.dof.size;
  t.stop();

  auto feList = std::tuple{feSpaceVel, feSpaceP};
  auto assembler = make_assembler(feList);
  // auto assembler = make_assembler(std::forward_as_tuple(feSpaceU, feSpaceU,
  // feSpaceP));

  t.start("bc");
  auto const zero = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  // fully developed
  auto const inlet = [](Vec3 const & p) { return Vec2{0.0, 0.5 * (1. - p[0] * p[0])}; };
  // // constant
  // auto const inlet = [] (Vec3 const & ) { return Vec2{0.0, 1.0}; };
  // // step
  // auto const inlet = [] (Vec3 const &p) {
  //   return p[0] < .5 ? Vec2{0., 1.} : Vec2{0., 0.};
  // };
  auto const bcsVel = std::vector{
      BCEss{feSpaceVel, side::RIGHT, zero},
      BCEss{feSpaceVel, side::TOP, zero, {0}},
      BCEss{feSpaceVel, side::LEFT, zero, {0}},
      // BCEss{feSpaceVel, side::LEFT, zero},
      BCEss{feSpaceVel, side::BOTTOM, inlet},
  };
  // auto const bcsVel = std::vector{
  //     BCEss{feSpaceVel, mesh->physicalNames["2"], side::RIGHT, zero},
  //     BCEss{feSpaceVel, mesh->physicalNames["3"], zero, {0}},
  //     // BCEss{feSpaceVel, side::LEFT, zero, {0}},
  //     BCEss{feSpaceVel, mesh->physicalNames["4"], zero},
  //     BCEss{feSpaceVel, mesh->physicalNames["1"], inlet},
  // };

  auto const bcsP = std::vector<BCEss<FESpaceP_T>>{
      // BCEss{feSpaceP, side::TOP, [](Vec3 const &){ return 0.0; }},
  };
  t.stop();

  t.start("build");
  auto const nu = config["nu"].as<double>();
  auto const stiffness = AssemblyTensorStiffness{nu, feSpaceVel};
  auto const grad = AssemblyGrad{-1.0, feSpaceVel, feSpaceP};
  auto const div = AssemblyDiv{-1.0, feSpaceP, feSpaceVel};

  Builder builder{dofU * dim + dofP};
  builder.buildLhs(std::tuple{stiffness}, bcsVel);
  builder.buildCoupling(grad, bcsVel, bcsP);
  builder.buildCoupling(div, bcsP, bcsVel);
  builder.closeMatrix();
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

  // Utils::filelog << "A:\n" << builder.A << std::endl;
  // Utils::filelog << "b:\n" << builder.b << std::endl;
  // Utils::filelog << "sol:\n" << sol << std::endl;

  t.start("exact");
  auto exact = Vec{dim * dofU + dofP};
  interpolateAnalyticFunction(inlet, feSpaceVel, exact);
  interpolateAnalyticFunction(
      [nu](Vec3 const & p) { return nu * (1. - p(1)); }, feSpaceP, exact);
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
  IOManager ioVel{feSpaceVel, "output_stokes2dquad/vel"};
  ioVel.print({vel, velExact});
  IOManager ioP{feSpaceP, "output_stokes2dquad/p"};
  ioP.print({p, pExact});
  IOManager ioComponent{feSpaceComponent, "output_stokes2dquad/vel"};
  ioComponent.print({u, v, uExact, vExact});
  IOManager ioWSS{feSpaceP0, "output_stokes2dquad/wss"};
  ioWSS.print({wssCell, wssFacet, dist});
  t.stop();

  t.print();

  auto const uError = (u.data - uExact.data).norm();
  auto const vError = (v.data - vExact.data).norm();
  auto const pError = (p.data - pExact.data).norm();

  auto const wssCellExact = nu * (1.0 - dist.data[wssCellMaxPos]);
  auto const wssFacetExact = nu;

  return checkError(
      {uError, vError, pError, wssCellMaxVal, wssFacetMaxVal},
      {8.0423101978505405e-14,
       1.6859789996431204e-13,
       3e-14,
       wssCellExact,
       wssFacetExact});
}
