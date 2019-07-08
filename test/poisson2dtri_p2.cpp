#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"

#include <iostream>

using Elem_T = Triangle;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,2>::RefFE_T,
                          FEType<Elem_T,2>::RecommendedQR>;

static scalarFun_T rhs = [] (Vec3 const& p)
{
  return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};

int main(int argc, char* argv[])
{
  uint const numElemsX = (argc < 3)? 10 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3)? 10 : std::stoi(argv[2]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});

  FESpace_T feSpace{*mesh};

  BCList bcs{feSpace};
  bcs.addBC(BCEss{feSpace, side::LEFT, [] (Vec3 const&) {return 0.;}});
  bcs.addBC(BCEss{feSpace, side::BOTTOM, [] (Vec3 const&) {return 0.;}});

  Builder builder{feSpace.dof.size};
  builder.buildLhs(AssemblyStiffness(1.0, feSpace), bcs);
  builder.buildRhs(AssemblyAnalyticRhs(rhs, feSpace), bcs);
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

  IOManager io{feSpace, "output_poisson2dtri_p2/sol"};
  io.print({sol, exact, error});

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  if(std::fabs(norm - 0.004067660888944465) > 1.e-14)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
