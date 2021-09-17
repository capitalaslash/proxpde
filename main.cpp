#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

enum class SolverType
{
  CHOLESKY,
  BICGSTAB,
  SPARSELU
};

int main()
{
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<Mesh_T, LagrangeFE<Elem_T, 1>::RefFE_T, GaussQR<Elem_T, 4>>;

  auto const solverType = SolverType::SPARSELU;

  auto const rhs = [](Vec3 const & p) { return M_PI * std::sin(M_PI * p(0)); };

  auto const exactSol = [](Vec3 const & p)
  { return std::sin(M_PI * p(0)) / M_PI + p(0); };

  // auto const rhs = [] (Vec3 const& p) { return p(0); };
  // auto const exact_sol = [] (Vec3 const& p) { return 0.5*p(0) -  p(0)*p(0)*p(0)/6.;
  // };

  // auto const rhs = [] (Vec3 const&) { return 8.; };
  // auto const exact_sol = [] (Vec3 const& p) { return 4.*p(0)*(2.-p(0)); };
  // auto const exact_sol = [] (Vec3 const& p) { return 4.*p(0)*(1.-p(0)); };

  uint const numElemsX = 20;
  uint const numElemsY = 1;

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0.02, 0.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});
  // std::cout << *mesh << std::endl;

  FESpace_T feSpace{*mesh};

  auto bc = BCEss{feSpace, side::LEFT};
  bc << [](Vec3 const &) { return 0.; };
  auto const bcs = std::tuple{bc};

  AssemblyStiffness stiffness{1.0, feSpace};
  AssemblyAnalyticRhs f{rhs, feSpace};

  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{stiffness}, bcs);
  builder.buildRhs(std::tuple{f}, bcs);
  builder.closeMatrix();

  // std::cout << "builder.A:\n" << builder.A << std::endl;
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

  Var sol{"u"};
  switch (solverType)
  {
  case SolverType::CHOLESKY:
  {
    Eigen::SimplicialCholesky<Mat<StorageType::ClmMajor>> solver(builder.A);
    sol.data = solver.solve(builder.b);
    break;
  }
  case SolverType::BICGSTAB:
  {
    Eigen::BiCGSTAB<Mat<StorageType::RowMajor>> solver(builder.A);
    sol.data = solver.solve(builder.b);
    break;
  }
  case SolverType::SPARSELU:
  {
    Eigen::SparseLU<Mat<StorageType::ClmMajor>, Eigen::COLAMDOrdering<int>> solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(builder.A);
    // Compute the numerical factorization
    solver.factorize(builder.A);

    sol.data = solver.solve(builder.b);
    break;
  }
  }
  // std::cout<< "sol:\n" << sol << std::endl;

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);

  Var error{"error"};
  error.data = sol.data - exact.data;
  std::cout << "error: " << error.data.norm() << std::endl;

  IOManager io{feSpace, "sol"};
  io.print({sol, exact, error});

  return 0;
}
