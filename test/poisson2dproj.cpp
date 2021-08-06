#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

int main(int argc, char* argv[])
{
  using Elem_T = Triangle;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<Mesh_T,
                            LagrangeFE<Elem_T,1>::RefFE_T,
                            LagrangeFE<Elem_T,1>::RecommendedQR>;

  scalarFun_T const rhs = [] (Vec3 const& p)
  {
    return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  };

  scalarFun_T const exactSol = [] (Vec3 const& p)
  {
    return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  };

  uint const numElemsX = (argc < 3)? 10 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3)? 10 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});

  FESpace_T feSpace{*mesh};

  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [] (Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [] (Vec3 const &) { return 0.; };
  auto const bcs = std::tuple{bcLeft, bcBottom};

  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  Var rhsVec{"rhs", feSpace.dof.size};
  interpolateAnalyticFunction(rhs, feSpace, rhsVec.data);
  builder.buildRhs(std::tuple{AssemblyProjection(1.0, rhsVec.data, feSpace)}, bcs);
  builder.closeMatrix();

  Builder builderTest{feSpace.dof.size};
  builderTest.buildLhs(std::tuple{AssemblyStiffness(1.0, feSpace)}, bcs);
  builderTest.buildRhs(std::tuple{AssemblyProjection(1.0, rhsVec.data, feSpace)}, bcs);
  builder.closeMatrix();
  // std::cout << "b VecRhs:\n" << builder.b << "\nb Projection:\n" << builderTest.b << std::endl;

  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = exact.data - sol.data;

  IOManager io{feSpace, "output_poisson2dproj/sol"};
  io.print({sol, exact, error, rhsVec});

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  return checkError({norm}, {0.141867443977535});
}
