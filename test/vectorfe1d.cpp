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
  uint constexpr numVars = 2U;

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR,
      numVars>;

  auto const rhs = [](Vec3 const & p)
  {
    FVec<numVars> rhs;
    for (uint d = 0; d < numVars; ++d)
      rhs[d] = M_PI * std::sin(M_PI * p(0));
    return rhs;
  };

  scalarFun_T const exactSol = [](Vec3 const & p)
  { return numVars + std::sin(M_PI * p(0)) / M_PI + p(0); };

  MilliTimer t;
  uint const numElems = (argc < 2) ? 3 : std::stoi(argv[1]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start();
  buildHyperCube(*mesh, origin, length, {{numElems, 0, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  //  // rotation matrix
  //  double theta = M_PI / 3.;
  //  FMat<3,3> R;
  //  R << std::cos(theta), std::sin(theta), 0.0,
  //      -std::sin(theta), std::cos(theta), 0.0,
  //      0.0, 0.0, 1.0;
  //  auto Rt = R.transpose();

  //  // rotate mesh
  //  for (auto & p: mesh->pointList)
  //  {
  //    p.coord = R * p.coord;
  //  }

  t.start();
  FESpace_T feSpace{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto bcValue = [](Vec3 const &)
  {
    FVec<numVars> vars;
    for (uint v = 0; v < numVars; ++v)
    {
      vars[v] = v;
    }
    return vars;
  };
  auto bc = BCEss{feSpace, side::LEFT};
  bc << bcValue;
  auto const bcs = std::tuple{bc};
  std::cout << "bcs: " << t << " ms" << std::endl;

  t.start();
  auto const size = numVars * feSpace.dof.size;
  auto const stiffness = AssemblyStiffness{1.0, feSpace};
  // auto rotatedRhs = [&Rt] (Vec3 const& p) {return rhs(Rt * p);};
  // AssemblyAnalyticRhs f{rotatedRhs, feSpace};
  auto const f = AssemblyAnalyticRhs{rhs, feSpace};
  auto const bcNat = AssemblyBCNatural{bcValue, side::RIGHT, feSpace};
  Builder builder{size};
  builder.buildLhs(std::tuple{stiffness}, bcs);
  builder.buildRhs(
      std::tuple{
          f,
          bcNat,
      },
      bcs);
  builder.closeMatrix();
  std::cout << "fe assembly: " << t << " ms" << std::endl;

  // filelog << "A:\n" << builder.A << std::endl;
  // filelog << "b:\n" << builder.b << std::endl;

  t.start();
  Var sol{"sol"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  std::cout << "solve: " << t << " ms" << std::endl;

  std::cout << "sol:\n" << sol.data << std::endl;

  //  Var exact{"exact"};
  //  // auto rotatedESol = [&Rt] (Vec3 const& p) {return exact_sol(Rt * p);};
  //  // interpolateAnalyticFunction(rotatedESol, feSpace, exact.data);
  //  interpolateAnalyticFunction(exact_sol, feSpace, exact.data);
  //  Var error{"e"};
  //  error.data = sol.data - exact.data;

  t.start();
  IOManager io{feSpace, "output_vectorfe1d/sol"};
  io.print({sol}); // io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  //  double norm = error.data.norm();
  //  std::cout << "the norm of the error is " << norm << std::endl;
  //  if(std::fabs(norm - 2.61664e-11) > 1.e-10)
  //  {
  //    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
  //    return 1;
  //  }

  return 0;
}
