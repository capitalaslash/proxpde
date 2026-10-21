#include <cmath>       // for sin, M_PI, sqrt
#include <cstdlib>     // for uint
#include <filesystem>  // for path
#include <memory>      // for unique_ptr, construct_at, destr...
#include <optional>    // for optional
#include <string>      // for char_traits, stoi, basic_string
#include <tuple>       // for tuple, get
#include <type_traits> // for is_base_of_v
#include <vector>      // for vector

#include <fmt/base.h> // for println

#include <proxpde/assembly_lhs.hpp>  // for AssemblyStiffness
#include <proxpde/assembly_rhs.hpp>  // for AssemblyRhsAnalytic
#include <proxpde/bc.hpp>            // for BCEss
#include <proxpde/builder.hpp>       // for Builder
#include <proxpde/def.hpp>           // for Vec3, checkError, scalarFun_T
#include <proxpde/fe.hpp>            // for LagrangeFE
#include <proxpde/fespace.hpp>       // for FESpace, interpolateAnalyticFun...
#include <proxpde/geo.hpp>           // for Triangle
#include <proxpde/iomanager.hpp>     // for IOManager
#include <proxpde/mesh.hpp>          // for side, Mesh, buildHyperCube
#include <proxpde/reffe.hpp>         // for RefTriangleP1
#include <proxpde/sparse_matrix.hpp> // for LUSolver
#include <proxpde/var.hpp>           // for Var

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Triangle;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  scalarFun_T const rhs = [](Vec3 const & p)
  {
    return 2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
  };

  scalarFun_T const exactSol = [](Vec3 const & p)
  { return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1)); };

  uint const numElemsX = (argc < 3) ? 10 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3) ? 10 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});

  FESpace_T feSpace{*mesh};

  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::vector{bcLeft, bcBottom};

  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  builder.buildRhs(std::tuple{AssemblyRhsAnalytic(rhs, feSpace)}, bcs);
  // using an interpolated rhs makes its quality independent of the chosen qr
  // Vec rhsProj;
  // interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  // builder.buildRhs(AssemblyProjection(1.0, rhsProj, feSpace), bcs);
  builder.closeMatrix();

  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  IOManager io{feSpace, "output_poisson2dtri/sol"};
  io.print({sol, exact, error});

  double norm = error.data.norm();
  fmt::println("the norm of the error is {:.16e}", norm);
  return checkError({norm}, {0.05950336122297865});
}
