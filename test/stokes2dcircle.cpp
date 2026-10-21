#include <cmath>         // for sqrt, fabs
#include <cstdlib>       // for uint
#include <filesystem>    // for path
#include <memory>        // for unique_ptr, construct_at
#include <optional>      // for optional
#include <string>        // for stoi, basic_string
#include <tuple>         // for tuple, get
#include <type_traits>   // for is_base_of_v
#include <unordered_set> // for unordered_set
#include <vector>        // for vector

#include <fmt/base.h> // for println

#include <proxpde/assembly_bc.hpp>       // for AssemblyBCNaturalAnalytic
#include <proxpde/assembly_coupling.hpp> // for AssemblyDiv, AssemblyGrad
#include <proxpde/assembly_lhs.hpp>      // for AssemblyDummy, AssemblyStif...
#include <proxpde/bc.hpp>                // for BCEss, DOFCoordSet
#include <proxpde/builder.hpp>           // for Builder
#include <proxpde/def.hpp>               // for Vec3, Vec2, checkError
#include <proxpde/fe.hpp>                // for LagrangeFE
#include <proxpde/fespace.hpp>           // for FESpace, integrateOnBoundary
#include <proxpde/geo.hpp>               // for Quad
#include <proxpde/iomanager.hpp>         // for IOManager
#include <proxpde/mesh.hpp>              // for side, Mesh, buildQuarterCircle
#include <proxpde/ns.hpp>                // for computeElemWSS
#include <proxpde/reffe.hpp>             // for RefQuadQ2, RefQuadQ1
#include <proxpde/sparse_matrix.hpp>     // for LUSolver
#include <proxpde/timer.hpp>             // for MilliTimer
#include <proxpde/var.hpp>               // for Var

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpaceU_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 2>::RefFE_T,
      LagrangeFE<Elem_T, 2>::RecommendedQR,
      2>;
  using FESpaceP_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 2>::RecommendedQR>;

  MilliTimer t;

  t.start("mesh");
  uint const numElemsX = (argc < 4) ? 6 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 4) ? 6 : std::stoi(argv[2]);
  uint const numElemsR = (argc < 4) ? 6 : std::stoi(argv[3]);

  Vec3 const origin{0., 0., 0.};
  double const radius = 1.;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildQuarterCircle(*mesh, origin, radius, {{numElemsX, numElemsY, numElemsR}});
  t.stop();

  t.start("feSpace");
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpaceP{*mesh, FESpaceU_T::dim * feSpaceU.dof.size};
  t.stop();

  t.start("bcs");
  auto const zeroFun = [](Vec3 const &) { return Vec2{0., 0.}; };
  auto const bcsU = std::vector{
      BCEss{feSpaceU, side::CIRCLE, zeroFun},
      BCEss{feSpaceU, side::TOP, zeroFun, {0}},
      BCEss{feSpaceU, side::LEFT, zeroFun, {1}}};
  DOFCoordSet pinSet{
      feSpaceP,
      [](Vec3 const & p)
      { return std::fabs(p[0] - 1.0) < 1.e-6 && std::fabs(p[1] - 0.0) < 1.e-6; }};
  auto const bcPin = BCEss{feSpaceP, pinSet.ids, [](Vec3 const &) { return 0.; }};
  auto const bcsP = std::vector{bcPin};
  t.stop();

  t.start("assembly");
  double const nu = 0.1;
  AssemblyStiffness stiffness(nu, feSpaceU);
  AssemblyGrad grad(-1.0, feSpaceU, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceU);
  // this is required to properly apply the pinning on the pressure
  AssemblyDummy dummy{feSpaceP};

  uint const numDOFs = 2 * feSpaceU.dof.size + feSpaceP.dof.size;
  Builder builder{numDOFs};
  builder.buildLhs(std::tuple{stiffness}, bcsU);
  builder.buildLhs(std::tuple{dummy}, bcsP);
  builder.buildCoupling(grad, bcsU, bcsP);
  builder.buildCoupling(div, bcsP, bcsU);
  builder.buildRhs(
      std::tuple{AssemblyBCNaturalAnalytic{
          [](Vec3 const &) { return Vec2{1., 0.}; }, side::LEFT, feSpaceU}},
      bcsU);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol{"vel", numDOFs};
  LUSolver solver(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol.data << std::endl;

  t.start("postpro");
  Var p{"p", sol.data, 2 * feSpaceU.dof.size, feSpaceP.dof.size};
  t.stop();

  // integral on boundary
  t.start("boundary int");
  double const pIntegral = integrateOnBoundary(p.data, feSpaceP, side::LEFT);
  fmt::println("integral of pressure on left face: {:.16e}", pIntegral);
  t.stop();

  // wall shear stress on circle
  t.start("wss");
  FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, LagrangeFE<Elem_T, 0>::RecommendedQR>
      feSpaceWSS{*mesh};
  // mean cell value
  Var wss{"wss"};
  computeElemWSS(wss.data, feSpaceWSS, sol.data, feSpaceU, {side::CIRCLE}, nu);
  t.stop();
  fmt::println("wss min: {:.16e}", wss.data.minCoeff());
  fmt::println("wss max: {:.16e}", wss.data.maxCoeff());

  // Var exact{"exact"};
  // interpolateAnalyticFunction(exact_sol, feSpaceU, exact.data);
  // Var error{"e"};
  // error.data = sol.data /*- exact.data*/;

  t.start("print");
  IOManager ioU{feSpaceU, "output_stokes2dcircle/vel"};
  ioU.print({sol});
  IOManager ioP{feSpaceP, "output_stokes2dcircle/p"};
  ioP.print({p});
  IOManager ioWSS{feSpaceWSS, "output_stokes2dcircle/wss"};
  ioWSS.print({wss});
  t.stop();

  t.print();

  // double const norm = error.data.norm();
  return checkError(
      {pIntegral, wss.data.maxCoeff()},
      {9.800943987495e-01, 3.068317649010e-01},
      1.e-12);
}
