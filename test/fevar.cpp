#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

int main(int argc, char* argv[])
{
  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<Mesh_T,
                            FEType<Elem_T, 1>::RefFE_T,
                            FEType<Elem_T, 1>::RecommendedQR>;

  static scalarFun_T rhs = [] (Vec3 const & p)
  {
    // return M_PI*std::sin(M_PI*p(0));
    // return 0.25 * M_PI * M_PI * std::sin(0.5 * M_PI * p[0]);
    return std::sin(0.5 * M_PI * p[0]) - 0.25 * M_PI * M_PI * cos(M_PI * p[0]);
  };

  static scalarFun_T exactSol = [] (Vec3 const & p)
  {
    // return std::sin(M_PI*p(0))/M_PI + p(0);
    return std::sin(0.5 * M_PI * p[0]);
  };

  MilliTimer t;

  t.start("mesh build");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  uint const numElems = (argc < 2)? 20 : std::stoi(argv[1]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  buildHyperCube(*mesh, origin, length, {{numElems, 0, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  BCList bcs{feSpace};
  bcs.addBC(BCEss{feSpace, side::LEFT, [] (Vec3 const&) {return 0.;}});
  t.stop();

  t.start("fe build");
  // AssemblyStiffness stiffness(1.0, feSpace);
  FEVar nu{"nu", feSpace};
  nu << [] (Vec3 const & p) { return std::sin(0.5 * M_PI * p[0]); };
  // nu << [] (Vec3 const & ) { return 1.; };
  ScalarCoef one{1.};
  AssemblyMassFE mass{1., one, feSpace};
  AssemblyStiffnessFE stiffness{nu, feSpace};
  AssemblyAnalyticRhs f(rhs, feSpace);
  // // using an interpolated rhs makes its quality independent of the chosen qr
  // Vec rhsProj;
  // interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  // AssemblyProjection f(1.0, rhsProj, feSpace);

  Builder builder{feSpace.dof.size};
  builder.buildLhs(mass, bcs);
  builder.buildLhs(stiffness, bcs);
  builder.buildRhs(f, bcs);
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
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output/sol_fevar"};
  io.print({sol, exact, error});
  t.stop();

  t.print();

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  if(std::fabs(norm - 0.001115928841940259) > 1.e-15)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
