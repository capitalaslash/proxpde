#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "geo.hpp"
#include "iomanager.hpp"

int main()
{
  scalarFun_T const rhs = [](Vec3 const & p) { return M_PI * std::sin(M_PI * p(0)); };

  scalarFun_T const exactSol = [](Vec3 const & p)
  { return std::sin(M_PI * p(0)) / M_PI + p(0); };

  // const scalarFun_T rhs = [] (Vec3 const& p) { return p(0); };
  // const scalarFun_T exactSol = [] (Vec3 const& p) { return 0.5*p(0) -
  // p(0)*p(0)*p(0)/6.; };

  // const scalarFun_T rhs = [] (Vec3 const&) { return 8.; };
  // const scalarFun_T exactSol = [] (Vec3 const& p) { return 4.*p(0)*(2.-p(0)); };
  // const scalarFun_T exactSol = [] (Vec3 const& p) { return 4.*p(0)*(1.-p(0)); };

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<Mesh_T, LagrangeFE<Elem_T, 1>::RefFE_T, TrapQR<Elem_T>>;

  uint const numElemsX = 10;
  uint const numElemsY = 2;

  Vec3 const origin(0., 0., 0.);
  Vec3 const length(1., 1., 0.);

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});
  // std::cout << *mesh << std::endl;

  FESpace_T feSpace{*mesh};

  // bc setup
  auto const bcs = std::tuple{};

  AssemblyScalarMass mass{1.0, feSpace};

  AssemblyAnalyticRhs f{exactSol, feSpace};

  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{mass}, bcs);
  builder.buildRhs(std::tuple{f}, bcs);
  builder.closeMatrix();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;

  // std::ofstream fout("output.m");
  // for( int k=0; k<builder.A.outerSize(); k++)
  // {
  //   for (Mat::InnerIterator it(builder.A,k); it; ++it)
  //   {
  //     std::cout << it.row() << " " << it.col() << " " << it.value() << " " <<
  //     it.index() << std::endl; fout << it.row()+1 << " " << it.col()+1 << " " <<
  //     it.value() << std::endl;
  //   }
  //   std::cout << "-----" << std::endl;
  // }
  // std::cout << "=====" << std::endl;
  // fout.close();

  Var sol{"sol"};
  LUSolver solver;
  // Compute the ordering permutation vector from the structural pattern of builder.A
  solver.analyzePattern(builder.A);
  // Compute the numerical factorization
  solver.factorize(builder.A);

  sol.data = solver.solve(builder.b);
  std::cout << "sol:\n" << sol.data << std::endl;

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);

  Var error{"error"};
  error.data = sol.data - exact.data;

  IOManager io{feSpace, "output_trap/sol"};
  io.print({sol, exact, error});

  std::cout << "error: " << error.data.norm() << std::endl;

  return 0;
}
