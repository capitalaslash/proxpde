#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timer.hpp"

#include <unsupported/Eigen/src/IterativeSolvers/Scaling.h>

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Hexahedron;
  using Mesh_T = Mesh<Elem_T>;
  uint constexpr order = 2U;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, order>::RefFE_T,
      LagrangeFE<Elem_T, order>::RecommendedQR>;

  ParameterDict config;
  config["mesh"]["origin"] = Vec3{0., 0., 0.};
  config["mesh"]["length"] = Vec3{1., 1., 1.};
  config["mesh"]["n"] = std::array{10U, 10U, 10U};
  config["solver"]["package"] = "eigen";
  config["solver"]["type"] = "iterative";
  config["solver"]["tolerance"] = 1.e-8;
  config["solver"]["max_iters"] = 2000;

  if (argc > 1)
  {
    config.override(argv[1]);
  }

  const scalarFun_T rhs = [](Vec3 const & p)
  {
    return 2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
  };
  const scalarFun_T exactSol = [](Vec3 const & p)
  { return std::sin(0.5 * M_PI * p(0)) * std::sin(1.5 * M_PI * p(1)); };

  MilliTimer t;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start("mesh build");
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [](Vec3 const &) { return 0.; };
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << [](Vec3 const &) { return 0.; };
  auto const bcs = std::vector{bcLeft, bcBottom};
  t.stop();

  t.start("fe build");
  Builder<StorageType::RowMajor> builderR{feSpace.dof.size};
  builderR.buildLhs(std::tuple{AssemblyStiffness{1.0, feSpace}}, bcs);
  builderR.buildRhs(std::tuple{AssemblyAnalyticRhs{rhs, feSpace}}, bcs);
  builderR.closeMatrix();
  fmt::print("the matrix is compressed? {}\n", builderR.A.isCompressed());
  t.stop();

  // filelog << "A:\n" << builder.A << std::endl;
  // filelog << "b:\n" << builder.b << std::endl;

  t.start("iterative solver");
  // generic iterative solver
  auto solver = getSolver(ParameterDict{config["solver"]});
  solver->compute(builderR.A);
  Vec u;
  auto const [iters, err] = solver->solve(builderR.b, u);
  fmt::print("iterative solver: iterations = {}, error = {}\n", iters, err);
  t.stop();

  using DiagPrec = Eigen::DiagonalPreconditioner<double>;
  // using ILUTPrec = Eigen::IncompleteLUT<double>;

  t.start("BiCGSTAB DiagPrec");
  Var solCG{"uBiCGSTAB"};
  config["solver"]["type"] = "bicgstab";
  config["solver"]["preconditioner"] = "diag";
  auto solverCG = getSolver(ParameterDict{config["solver"]});
  solverCG->compute(builderR.A);
  auto const [itersCG, errCG] = solverCG->solve(builderR.b, solCG.data);
  fmt::print("BiCGSTAB DiagPrec\niters: {}\nerror: {}\n", itersCG, errCG);
  t.stop();

  t.start("MINRES DiagPrec");
  Var solMINRES{"uMINRES"};
  config["solver"]["type"] = "minres";
  config["solver"]["preconditioner"] = "diag";
  auto solverMINRES = getSolver(ParameterDict{config["solver"]});
  solverMINRES->compute(builderR.A);
  auto const [itersMINRES, errMINRES] = solverMINRES->solve(builderR.b, solMINRES.data);
  fmt::print("MINRES DiagPrec\niters: {}\nerror: {}\n", itersMINRES, errMINRES);
  t.stop();

  t.start("GMRES DiagPrec");
  Var solGMRES{"uGMRES"};
  config["solver"]["type"] = "gmres";
  config["solver"]["preconditioner"] = "diag";
  auto solverGMRES = getSolver(ParameterDict{config["solver"]});
  solverGMRES->compute(builderR.A);
  auto const [itersGMRES, errGMRES] = solverGMRES->solve(builderR.b, solGMRES.data);
  fmt::print("GMRES DiagPrec\niters: {}\nerror: {}\n", itersGMRES, errGMRES);
  t.stop();

  // crashes in opt build
  // t.start("DGMRES DiagPrec");
  // std::cout << "DGMRES DiagPrec" << std::endl;
  // {
  //   Eigen::DGMRES<Mat, DiagPrec> solver;
  //   solver.setTolerance(config["solver"]["tolerance"].as<double>());
  //   solver.setMaxIterations(config["solver"]["max_iters"].as<int>());
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
    solver.setTolerance(config["solver"]["tolerance"].as<double>());
    solver.setMaxIterations(config["solver"]["max_iters"].as<int>());
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
  // ILUTPrec>>(builderR.A, builderR.b, config["solver"], solCGILUT.data); t.stop();
  //
  // t.start("MINRES ILUTPrec");
  // std::cout << "MINRES ILUTPrec" << std::endl;
  // Var solMINRESILUT{"uMINRESILUT"};
  // solve<StorageType::RowMajor, Eigen::MINRES<Mat<StorageType::RowMajor>,
  // Eigen::Lower, ILUTPrec>>(builderR.A, builderR.b, config["solver"],
  // solMINRESILUT.data); t.stop();
  //
  // t.start("GMRES ILUTPrec");
  // std::cout << "GMRES ILUTPrec" << std::endl;
  // Var solGMRESILUT{"uGMRESILUT"};
  // solve<StorageType::RowMajor, Eigen::GMRES<Mat<StorageType::RowMajor>,
  // ILUTPrec>>(builderR.A, builderR.b, config["solver"], solGMRESILUT.data); t.stop();

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
