#include <cmath>         // for sqrt
#include <cstdlib>       // for uint
#include <filesystem>    // for path, operator/
#include <memory>        // for unique_ptr, construct_at
#include <optional>      // for optional
#include <string>        // for basic_string
#include <tuple>         // for tuple, get
#include <type_traits>   // for is_base_of_v
#include <unordered_set> // for unordered_set
#include <vector>        // for vector

#include <fmt/base.h>      // for println
#include <yaml-cpp/yaml.h> // for Node, LoadFile

#include <proxpde/assembler.hpp>         // for make_assembler
#include <proxpde/assembly_coupling.hpp> // for AssemblyDiv, AssemblyGrad
#include <proxpde/assembly_lhs.hpp>      // for AssemblyTensorStiffness
#include <proxpde/bc.hpp>                // for BCEss
#include <proxpde/builder.hpp>           // for Builder
#include <proxpde/def.hpp>               // for Vec3, ParameterDict, Vec2
#include <proxpde/fe.hpp>                // for LagrangeFE
#include <proxpde/fespace.hpp>           // for FESpace, interpolateAnalyti...
#include <proxpde/geo.hpp>               // for Triangle
#include <proxpde/iomanager.hpp>         // for IOManager
#include <proxpde/mesh.hpp>              // for Mesh, side, MeshType, readMesh
#include <proxpde/ns.hpp>                // for computeElemWSS
#include <proxpde/reffe.hpp>             // for RefTriangleP1, RefTriangleP2
#include <proxpde/sparse_matrix.hpp>     // for LUSolver
#include <proxpde/timer.hpp>             // for MilliTimer
#include <proxpde/var.hpp>               // for Var

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Triangle;
  constexpr uint dim = Elem_T::dim;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE_T = LagrangeFE<Elem_T, 2u>::RefFE_T;
  using LinearRefFE_T = LagrangeFE<Elem_T, 1u>::RefFE_T;
  using QuadraticQR_T = LagrangeFE<Elem_T, 2u>::RecommendedQR;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE_T, QuadraticQR_T, dim>;
  // using FESpaceComponent_T = FESpace<Mesh_T, QuadraticRefFE_T, QuadraticQR_T>;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE_T, QuadraticQR_T>;

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
    config["mesh"]["filename"] = (pathExe / "square_uns.msh").string();
  }

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, config["mesh"]);
  t.stop();

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  auto const dofU = feSpaceVel.dof.size;
  FESpaceP_T feSpaceP{*mesh, dim * dofU};
  auto const dofP = feSpaceP.dof.size;
  t.stop();

  auto feList = std::tuple{feSpaceVel, feSpaceP};
  auto assembler = make_assembler(feList);
  // auto assembler = make_assembler(std::forward_as_tuple(feSpaceU, feSpaceU,
  // feSpaceP));

  t.start("bc");
  auto const zero = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  // fully developed
  auto const inlet = [](Vec3 const & p) { return Vec2{0.0, 0.5 * (1. - p(0) * p(0))}; };
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
  auto const bcsP = std::vector<BCEss<FESpaceP_T>>{};
  t.stop();

  t.start("build");
  double const nu = 1.0;
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
  auto pExact = Var{"pExact"};
  pExact.data = exact.tail(dofP);

  // wall shear stress
  t.start("wssCell");
  FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, LagrangeFE<Elem_T, 0>::RecommendedQR>
      feSpaceP0{*mesh};
  // mean cell value
  Var wssCell{"wssCell"};
  computeElemWSS(wssCell.data, feSpaceP0, vel.data, feSpaceVel, {side::RIGHT}, nu);
  t.stop();
  fmt::println("wssCell min: {:.16e}", wssCell.data.minCoeff());
  fmt::println("wssCell max: {:.16e}", wssCell.data.maxCoeff());

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_stokes2dtri/vel"};
  ioVel.print({vel, velExact});
  IOManager ioP{feSpaceP, "output_stokes2dtri/p"};
  ioP.print({p, pExact});
  t.stop();

  t.print();

  auto const velError = (vel.data - velExact.data).norm();
  auto const pError = (p.data - pExact.data).norm();

  return checkError(
      {velError, pError}, {1.2722955943638614e-14, 1.9053041663391855e-13}, 1e-16);
}
