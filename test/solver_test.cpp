#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

#include <unsupported/Eigen/src/IterativeSolvers/Scaling.h>
#include <yaml-cpp/yaml.h>

using Elem_T = Hexahedron;
uint constexpr order = 2U;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<
    Mesh_T,
    LagrangeFE<Elem_T, order>::RefFE_T,
    LagrangeFE<Elem_T, order>::RecommendedQR>;

template <StorageType Storage, typename Solver>
void solve(Mat<Storage> const & A, Vec const & b, YAML::Node const & config, Vec & x)
{
  Solver solver;
  solver.setTolerance(config["tolerance"].as<double>());
  solver.setMaxIterations(config["maxNumIters"].as<int>());
  solver.compute(A);
  x = solver.solve(b);
  std::cout << "iters: " << solver.iterations() << std::endl;
  std::cout << "error: " << solver.error() << std::endl;
}

int main(int argc, char * argv[])
{
  const scalarFun_T rhs = [](Vec3 const & p)
  {
    return 2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
  };
  const scalarFun_T exactSol = [](Vec3 const & p)
  { return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1)); };

  MilliTimer t;
  uint const numElems = (argc < 2) ? 10 : std::stoi(argv[1]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 1.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start("mesh build");
  buildHyperCube(*mesh, origin, length, {{numElems, numElems, numElems}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::tuple{bcLeft, bcBottom};
  t.stop();

  t.start("fe build");
  AssemblyStiffness stiffness{1.0, feSpace};
  AssemblyAnalyticRhs f{rhs, feSpace};
  Builder<StorageType::RowMajor> builderR{feSpace.dof.size};
  builderR.buildLhs(std::tuple{stiffness}, bcs);
  builderR.buildRhs(std::tuple{f}, bcs);
  builderR.closeMatrix();
  t.stop();

  // filelog << "A:\n" << builder.A << std::endl;
  // filelog << "b:\n" << builder.b << std::endl;

  YAML::Node config;
  config["tolerance"] = 1.e-8;
  config["maxNumIters"] = 2000;
  using DiagPrec = Eigen::DiagonalPreconditioner<double>;
  // using ILUTPrec = Eigen::IncompleteLUT<double>;

  t.start("BiCGSTAB DiagPrec");
  std::cout << "BiCGSTAB DiagPrec" << std::endl;
  Var solCG{"uBiCGSTAB"};
  solve<StorageType::RowMajor, Eigen::BiCGSTAB<Mat<StorageType::RowMajor>, DiagPrec>>(
      builderR.A, builderR.b, config, solCG.data);
  t.stop();

  t.start("MINRES DiagPrec");
  std::cout << "MINRES DiagPrec" << std::endl;
  Var solMINRES{"uMINRES"};
  solve<
      StorageType::RowMajor,
      Eigen::MINRES<Mat<StorageType::RowMajor>, Eigen::Lower, DiagPrec>>(
      builderR.A, builderR.b, config, solMINRES.data);
  t.stop();

  t.start("GMRES DiagPrec");
  std::cout << "GMRES DiagPrec" << std::endl;
  Var solGMRES{"uGMRES"};
  solve<StorageType::RowMajor, Eigen::GMRES<Mat<StorageType::RowMajor>, DiagPrec>>(
      builderR.A, builderR.b, config, solGMRES.data);
  t.stop();

  // crashes in opt build
  // t.start("DGMRES DiagPrec");
  // std::cout << "DGMRES DiagPrec" << std::endl;
  // {
  //   Eigen::DGMRES<Mat, DiagPrec> solver;
  //   solver.setTolerance(config["tolerance"].as<double>());
  //   solver.setMaxIterations(config["maxNumIters"].as<int>());
  //   solver.set_restart(30); // Set restarting value
  //   solver.setEigenv(1); // Set the number of eigenvalues to deflate
  //   solver.compute(builder.A);
  //   Vec x = solver.solve(builder.b);
  //   std::cout << "iters: " << solver.iterations() << std::endl;
  //   std::cout << "error: " << solver.error() << std::endl;
  // }
  // t.stop();

  t.start("BiCGSTAB scaling");
  std::cout << "BiCGSTAB scaling" << std::endl;
  {
    Vec x;
    Eigen::IterScaling<Mat<StorageType::RowMajor>> scal;

    // Compute the left and right scaling vectors. The matrix is equilibrated at output
    scal.computeRef(builderR.A);

    // Scale the right hand side
    builderR.b = scal.LeftScaling().cwiseProduct(builderR.b);

    // Now, solve the equilibrated linear system with any available solver
    Eigen::BiCGSTAB<Mat<StorageType::RowMajor>, DiagPrec> solver;
    solver.setTolerance(config["tolerance"].as<double>());
    solver.setMaxIterations(config["maxNumIters"].as<int>());
    solver.compute(builderR.A);
    x = solver.solve(builderR.b);
    std::cout << "iters: " << solver.iterations() << std::endl;
    std::cout << "error: " << solver.error() << std::endl;

    // Scale back the computed solution
    x = scal.RightScaling().cwiseProduct(x);
  }
  t.stop();

  // t.start("BiCGSTAB ILUTPrec");
  // std::cout << "BiCGSTAB ILUTPrec" << std::endl;
  // Var solCGILUT{"uBiCGSTABILUT"};
  // solve<StorageType::RowMajor, Eigen::BiCGSTAB<Mat<StorageType::RowMajor>,
  // ILUTPrec>>(builderR.A, builderR.b, config, solCGILUT.data); t.stop();
  //
  // t.start("MINRES ILUTPrec");
  // std::cout << "MINRES ILUTPrec" << std::endl;
  // Var solMINRESILUT{"uMINRESILUT"};
  // solve<StorageType::RowMajor, Eigen::MINRES<Mat<StorageType::RowMajor>,
  // Eigen::Lower, ILUTPrec>>(builderR.A, builderR.b, config, solMINRESILUT.data);
  // t.stop();
  //
  // t.start("GMRES ILUTPrec");
  // std::cout << "GMRES ILUTPrec" << std::endl;
  // Var solGMRESILUT{"uGMRESILUT"};
  // solve<StorageType::RowMajor, Eigen::GMRES<Mat<StorageType::RowMajor>,
  // ILUTPrec>>(builderR.A, builderR.b, config, solGMRESILUT.data); t.stop();

  // t.start("solver LU");
  // Var solLU{"uLU"};
  // {
  //   LUSolver solver;
  //   solver.analyzePattern(builder.A);
  //   solver.factorize(builder.A);
  //   solLU.data = solver.solve(builder.b);
  // }
  // t.stop();

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = solGMRES.data - exact.data;

  // Var diff{"diff"};
  // diff.data = solGMRES.data - solLU.data;

  t.start("output");
  IOManager io{feSpace, "output/sol_gmres"};
  io.print({solGMRES, solMINRES, solCG, /*solLU, diff,*/ exact, error});
  t.stop();

  // filelog << "the test timings are meaningful only in opt mode with a grid" <<
  //            " size bigger than 500 x 500" << std::endl;
  t.print(/*filelog*/);

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(15) << norm
            << std::endl;
  return checkError({norm}, {6.579883739106e-04});
}
