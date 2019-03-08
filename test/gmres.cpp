#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

using Elem_T = Quad;
uint constexpr order = 1U;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T, order>::RefFE_T,
                          FEType<Elem_T, order>::RecommendedQR>;

static scalarFun_T rhs = [] (Vec3 const& p)
{
  return M_PI*std::sin(M_PI*p(0));
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return std::sin(M_PI*p(0))/M_PI + p(0);
};

int main(int argc, char* argv[])
{
  MilliTimer t;
  uint const numPts = (argc < 2)? 41 : std::stoi(argv[1]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start("mesh build");
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, origin, length, {{numPts, numPts, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  BCList bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [](Vec3 const &){return 0.;});
  t.stop();

  t.start("fe build");
  AssemblyStiffness stiffness{1.0, feSpace};
  AssemblyAnalyticRhs f{rhs, feSpace};
  Builder builder{feSpace.dof.size};
  builder.buildProblem(stiffness, bcs);
  builder.buildProblem(f, bcs);
  builder.closeMatrix();
  t.stop();

  // filelog << "A:\n" << builder.A << std::endl;
  // filelog << "b:\n" << builder.b << std::endl;

  t.start("solver LU");
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  t.start("solver GMRES");
  Var solNew{"uNew"};
  GMRESSolver solverNew;
  solverNew.compute(builder.A);
  solNew.data = solver.solve(builder.b);
  t.stop();

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{feSpace, "output/sol_gmres"};
  io.print({sol, solNew, exact, error});
  t.stop();

  filelog << "the test timings are meaningful only in opt mode with a grid" <<
             " size bigger than 500 x 500" << std::endl;
  t.print(filelog);

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(15) << norm << std::endl;
  if(std::fabs(norm - 1.57263551766678e-07) > 1.e-15)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
