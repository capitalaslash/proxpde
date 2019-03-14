#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

#include <yaml-cpp/yaml.h>
#include <unsupported/Eigen/src/IterativeSolvers/Scaling.h>

using Elem_T = Hexahedron;
uint constexpr order = 2U;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T, order>::RefFE_T,
                          FEType<Elem_T, order>::RecommendedQR>;

static scalarFun_T rhs = [] (Vec3 const& p)
{
  return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
};

template <typename Solver>
void solve(Mat const & A, Vec const & b, YAML::Node const & config, Vec & x)
{
  Solver solver;
  solver.setTolerance(config["tolerance"].as<double>());
  solver.setMaxIterations(config["maxNumIters"].as<int>());
  solver.compute(A);
  x = solver.solve(b);
  std::cout << "iters: " << solver.iterations() << std::endl;
  std::cout << "error: " << solver.error() << std::endl;
}

int main(int argc, char* argv[])
{
  MilliTimer t;
  uint const numPts = (argc < 2)? 11 : std::stoi(argv[1]);

  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 1.};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start("mesh build");
  MeshBuilder<Elem_T> meshBuilder;
  meshBuilder.build(*mesh, origin, length, {{numPts, numPts, numPts}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  BCList bcs{feSpace};
  bcs.addEssentialBC(side::LEFT, [](Vec3 const &){return 0.;});
  bcs.addEssentialBC(side::BOTTOM, [](Vec3 const &){return 0.;});
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

  YAML::Node config;
  config["tolerance"] = 1.e-8;
  config["maxNumIters"] = 2000;
  using DiagPrec = Eigen::DiagonalPreconditioner<Mat::Scalar>;
  using ILUTPrec = Eigen::IncompleteLUT<Mat::Scalar>;

  t.start("BiCGSTAB DiagPrec");
  std::cout << "BiCGSTAB DiagPrec" << std::endl;
  Var solCG{"uBiCGSTAB"};
  solve<Eigen::BiCGSTAB<Mat, DiagPrec>>(builder.A, builder.b, config, solCG.data);
  t.stop();

  t.start("MINRES DiagPrec");
  std::cout << "MINRES DiagPrec" << std::endl;
  Var solMINRES{"uMINRES"};
  solve<Eigen::MINRES<Mat, Eigen::Lower, DiagPrec>>(builder.A, builder.b, config, solMINRES.data);
  t.stop();

  t.start("GMRES DiagPrec");
  std::cout << "GMRES DiagPrec" << std::endl;
  Var solGMRES{"uGMRES"};
  solve<Eigen::GMRES<Mat, DiagPrec>>(builder.A, builder.b, config, solGMRES.data);
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
    Eigen::IterScaling<Mat> scal;

    // Compute the left and right scaling vectors. The matrix is equilibrated at output
    scal.computeRef(builder.A);

    // Scale the right hand side
    builder.b = scal.LeftScaling().cwiseProduct(builder.b);

    // Now, solve the equilibrated linear system with any available solver
    Eigen::BiCGSTAB<Mat, DiagPrec> solver;
    solver.setTolerance(config["tolerance"].as<double>());
    solver.setMaxIterations(config["maxNumIters"].as<int>());
    solver.compute(builder.A);
    x = solver.solve(builder.b);
    std::cout << "iters: " << solver.iterations() << std::endl;
    std::cout << "error: " << solver.error() << std::endl;

    // Scale back the computed solution
    x = scal.RightScaling().cwiseProduct(x);
  }
  t.stop();

  // t.start("BiCGSTAB ILUTPrec");
  // std::cout << "BiCGSTAB ILUTPrec" << std::endl;
  // Var solCGILUT{"uBiCGSTABILUT"};
  // solve<Eigen::BiCGSTAB<Mat, ILUTPrec>>(builder.A, builder.b, config, solCGILUT.data);
  // t.stop();
  //
  // t.start("MINRES ILUTPrec");
  // std::cout << "MINRES ILUTPrec" << std::endl;
  // Var solMINRESILUT{"uMINRESILUT"};
  // solve<Eigen::MINRES<Mat, Eigen::Lower, ILUTPrec>>(builder.A, builder.b, config, solMINRESILUT.data);
  // t.stop();
  //
  // t.start("GMRES ILUTPrec");
  // std::cout << "GMRES ILUTPrec" << std::endl;
  // Var solGMRESILUT{"uGMRESILUT"};
  // solve<Eigen::GMRES<Mat, ILUTPrec>>(builder.A, builder.b, config, solGMRESILUT.data);
  // t.stop();

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
  std::cout << "the norm of the error is " << std::setprecision(15) << norm << std::endl;
  if(std::fabs(norm - 0.00065811014197935) > 1.e-12)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
