#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
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

  MilliTimer t;

  t.start("mesh build");
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
    config["mesh"]["n"] = std::array{10U, 10U, 0U};
    // config["mesh"]["type"] = MeshType::GMSH;
    // config["mesh"]["filename"] = "square_q.msh";
  }

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, config["mesh"]);
  t.stop();

  t.start("fespace");
  // FESpace_T feSpace{*mesh};
  FESpace_T feSpace;
  feSpace.init(*mesh);
  t.stop();

  t.start("bcs");
  // auto bcLeft = BCEss{feSpace, side::LEFT};
  auto bcLeft = BCEss<FESpace_T>{};
  bcLeft.init(feSpace, side::LEFT);
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::vector{bcLeft, bcBottom};
  t.stop();

  t.start("fe space");
  AssemblyStiffness stiffness(1.0, feSpace);
  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{stiffness}, bcs);
  builder.buildRhs(std::tuple{AssemblyRhsAnalytic(rhs, feSpace)}, bcs);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output_poisson2dquad/sol"};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm
            << std::endl;
  return checkError({norm}, {0.02049777877937642});
}
