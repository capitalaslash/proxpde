#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "def.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

int main(/*int argc, char * argv[]*/)
{
  using namespace proxpde;

  using Elem_T = Triangle;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T =
      FESpace<Mesh_T, RefTriangleCR1, LagrangeFE<Elem_T, 1>::RecommendedQR>;
  using FESpacePP_T = FESpace<
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

  // uint const numElemsX = (argc < 3) ? 10 : std::stoi(argv[1]);
  // uint const numElemsY = (argc < 3) ? 10 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});
  readGMSH(*mesh, "square_uns.msh");

  FESpace_T feSpace{*mesh};

  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::vector{bcLeft, bcBottom};

  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  // builder.buildRhs(std::tuple{AssemblyAnalyticRhs(rhs, feSpace)}, bcs);
  // using an interpolated rhs makes its quality independent of the chosen qr
  Vec rhsProj;
  // interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  projectAnalyticFunction(rhs, feSpace, rhsProj);
  builder.buildRhs(std::tuple{AssemblyProjection(1.0, rhsProj, feSpace)}, bcs);
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

  FESpacePP_T feSpacePP{*mesh};
  Var solPP{"u"};
  Var exactPP{"exact"};
  l2Projection(solPP.data, feSpacePP, sol.data, feSpace);
  l2Projection(exactPP.data, feSpacePP, exact.data, feSpace);
  IOManager io{feSpacePP, "output_poisson2dcr/sol"};
  io.print({solPP, exactPP});

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm
            << std::endl;
  return checkError({norm}, {0.05950336122297865});
}
