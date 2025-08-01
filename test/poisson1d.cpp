#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1u>::RefFE_T,
      LagrangeFE<Elem_T, 1u>::RecommendedQR>;
  using FESpaceGrad_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 0u>::RefFE_T,
      LagrangeFE<Elem_T, 1u>::RecommendedQR>;
  using FESpaceSide_T =
      FESpace<Mesh_T, LagrangeFE<Elem_T, 1u>::RefFE_T, SideGaussQR<Line, 1u>>;

  scalarFun_T const rhs = [](Vec3 const & p) { return M_PI * std::sin(M_PI * p(0)); };

  scalarFun_T const exactSol = [](Vec3 const & p)
  { return std::sin(M_PI * p(0)) / M_PI + p(0); };

  MilliTimer t;

  uint const numElems = (argc < 2) ? 20u : std::stoi(argv[1]);

  auto const origin = Vec3{0.0, 0.0, 0.0};
  auto const length = Vec3{1.0, 0.0, 0.0};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start("mesh build");
  buildHyperCube(*mesh, origin, length, {{numElems, 0, 0}});

  // rotation matrix
  double theta = M_PI / 3.;
  FMat<3, 3> R;
  R << std::cos(theta), std::sin(theta), 0.0, -std::sin(theta), std::cos(theta), 0.0,
      0.0, 0.0, 1.0;
  auto Rt = R.transpose();

  // rotate mesh
  for (auto & p: mesh->pointList)
  {
    p.coord = R * p.coord;
  }
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  auto const bcs = std::vector{
      BCEss{feSpace, side::LEFT, [](Vec3 const &) { return 0.; }},
  };
  t.stop();

  t.start("fe build");
  AssemblyStiffness stiffness(1.0, feSpace);
  auto rotatedRhs = [&Rt, &rhs](Vec3 const & p) { return rhs(Rt * p); };
  AssemblyRhsAnalytic f(rotatedRhs, feSpace);
  // // using an interpolated rhs makes its quality independent of the chosen qr
  // Vec rhsProj;
  // interpolateAnalyticFunction(rotatedRhs, feSpace, rhsProj);
  // AssemblyProjection f(1.0, rhsProj, feSpace);

  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{stiffness}, bcs);
  builder.buildRhs(std::tuple{f}, bcs);
  builder.closeMatrix();
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  t.start("solve");
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "u:\n" << sol.data << std::endl;

  Var exact{"exact"};
  auto rotatedESol = [&Rt, &exactSol](Vec3 const & p) { return exactSol(Rt * p); };
  interpolateAnalyticFunction(rotatedESol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output/sol_poisson1d"};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  double const norm = error.data.norm();
  fmt::print("the norm of the error is {:.16e}\n", norm);
  return checkError({norm}, {2.87785419773588e-07});
}
