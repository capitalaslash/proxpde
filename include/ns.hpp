#pragma once

#include "def.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "var.hpp"
#include "iomanager.hpp"

// #include <amgcl/backend/eigen.hpp>
// #include <amgcl/make_solver.hpp>
// #include <amgcl/amg.hpp>
// #include <amgcl/coarsening/smoothed_aggregation.hpp>
// #include <amgcl/relaxation/as_preconditioner.hpp>
// #include <amgcl/relaxation/spai0.hpp>
// #include <amgcl/relaxation/spai1.hpp>
// #include <amgcl/relaxation/ilu0.hpp>
// #include <amgcl/relaxation/ilut.hpp>
// #include <amgcl/relaxation/gauss_seidel.hpp>
// #include <amgcl/relaxation/damped_jacobi.hpp>
// #include <amgcl/solver/bicgstab.hpp>
// #include <amgcl/solver/gmres.hpp>
// #include <amgcl/solver/skyline_lu.hpp>
// #include <amgcl/preconditioner/schur_pressure_correction.hpp>
// #include <amgcl/preconditioner/dummy.hpp>
// #include <amgcl/profiler.hpp>

std::vector<uint> const uComp = {0};
std::vector<uint> const vComp = {1};
std::vector<uint> const wComp = {2};
std::vector<uint> const uvComp = {0, 1};
std::vector<uint> const uwComp = {0, 2};
std::vector<uint> const vwComp = {1, 2};
std::vector<uint> const uvwComp = {0, 1, 2};

struct NSParameters
{
  double const dt;
  double const nu;
  fs::path const outputDir = "output";
  // avoid padding warning
  // long const dummy = 0;
};

template <typename Mesh>
struct NSSolverMonolithic
{
  static int constexpr dim = Mesh::Elem_T::dim;

  using Elem_T = typename Mesh::Elem_T;
  using FESpaceVel_T = FESpace<Mesh,
                               typename FEType<Elem_T, 2>::RefFE_T,
                               typename FEType<Elem_T, 2>::RecommendedQR,
                               dim>;
  using FESpaceP_T = FESpace<Mesh,
                             typename FEType<Elem_T, 1>::RefFE_T,
                             typename FEType<Elem_T, 2>::RecommendedQR>;
  // using Backend = amgcl::backend::eigen<double>;
  // using USolver = amgcl::make_solver<
  //     amgcl::amg<
  //         Backend,
  //         amgcl::coarsening::smoothed_aggregation,
  //         amgcl::relaxation::spai0
  //         // amgcl::relaxation::ilu0
  //         // amgcl::relaxation::gauss_seidel
  //     >,
  //     // amgcl::solver::bicgstab<Backend>
  //     amgcl::solver::gmres<Backend>
  // >;
  //
  // using PSolver = amgcl::make_solver<
  //     amgcl::amg<
  //         Backend,
  //         amgcl::coarsening::smoothed_aggregation,
  //         amgcl::relaxation::gauss_seidel
  //     >,
  //     amgcl::solver::gmres<Backend>
  // >;
  // // using PSolver = amgcl::make_solver<
  // //   amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>,
  // //   // amgcl::preconditioner::dummy<Backend>,
  // //   amgcl::solver::bicgstab<Backend>
  // // >;
  // // using PSolver = USolver;
  //
  // using SchurSolver = amgcl::make_solver<
  //     amgcl::preconditioner::schur_pressure_correction<USolver, PSolver>,
  //     amgcl::solver::gmres<Backend>
  // >;
  // // using SchurSolver = amgcl::make_solver<
  // //   amgcl::amg<
  // //     Backend,
  // //     amgcl::coarsening::smoothed_aggregation,
  // //     amgcl::relaxation::ilut
  // //   >,
  // //   // amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilut>,
  // //   // amgcl::preconditioner::dummy<Backend>,
  // //   amgcl::solver::gmres<Backend>
  // // >;
  //
  // using SolverParams = boost::property_tree::ptree;
  using SchurSolver = IterSolver;

  explicit NSSolverMonolithic(Mesh const & m, NSParameters const & par):
    mesh{m},
    parameters{par},
    feSpaceVel{m},
    feSpaceP{m},
    bcsVel{feSpaceVel},
    bcsP{feSpaceP},
    builder{feSpaceVel.dof.size * dim + feSpaceP.dof.size},
    sol{"vel", feSpaceVel.dof.size * dim + feSpaceP.dof.size},
    p{"p", feSpaceP.dof.size},
    velOld{feSpaceVel.dof.size * dim},
    assemblyRhs{1. / par.dt, velOld, feSpaceVel},
    assemblyAdvection{1.0, velOld, feSpaceVel, feSpaceVel},
    // pMask(feSpaceVel.dof.size * dim + feSpaceP.dof.size, 0),
    ioVel{feSpaceVel, par.outputDir / "vel"},
    ioP{feSpaceP, par.outputDir / "p"}
  {
    // // set up size and position of pressure dofs
    // std::fill(pMask.begin() + feSpaceVel.dof.size * dim, pMask.end(), 1);
    // solverParams.put("precond.pmask_size", pMask.size());
    // solverParams.put("precond.pmask", static_cast<void *>(pMask.data()));
    // // solverParams.put("precond.approx_schur", false);
    // solverParams.put("precond.psolver.precond.direct_coarse", false);
  }

  void init()
  {
    // TODO: assert that this comes after setting up bcs
    auto const dofU = feSpaceVel.dof.size;
    builder.buildProblem(AssemblyMass{1. / parameters.dt, feSpaceVel}, bcsVel);
    builder.buildProblem(AssemblyTensorStiffness{parameters.nu, feSpaceVel}, bcsVel);
    builder.buildProblem(AssemblyGrad{-1.0, feSpaceVel, feSpaceP, allComp<FESpaceVel_T>(), 0, dofU*dim}, bcsVel, bcsP);
    builder.buildProblem(AssemblyDiv{-1.0, feSpaceP, feSpaceVel, allComp<FESpaceVel_T>(), dofU*dim, 0}, bcsP, bcsVel);
    builder.buildProblem(AssemblyMass{0., feSpaceP, {0}, dofU*dim, dofU*dim}, bcsP);
    builder.closeMatrix();
    matFixed = builder.A;
    rhsFixed = builder.b;
  }

  void assemblyStep()
  {
    velOld = sol.data;
    builder.clear();
    builder.buildProblem(assemblyRhs, bcsVel);
    builder.buildProblem(assemblyAdvection, bcsVel);
    builder.closeMatrix();
    builder.A += matFixed;
    builder.b += rhsFixed;
  }

  void solve()
  {
    // // std::cout << builder.A << std::endl;
    // // std::cout << builder.b << std::endl;
    // auto const uSize = dim*feSpaceVel.dof.size;
    // auto const pSize = feSpaceP.dof.size;
    // Eigen::SparseMatrix<double> A  = builder.A.block(    0,     0, uSize, uSize);
    // Eigen::MatrixXd Bt = builder.A.block(    0, uSize, uSize, pSize);
    // Eigen::MatrixXd B  = builder.A.block(uSize,     0, pSize, uSize);
    // Eigen::MatrixXd M  = builder.A.block(uSize, uSize, pSize, pSize);
    // Eigen::VectorXd f  = builder.b.block(    0,     0, uSize,     1);
    // Eigen::VectorXd g  = builder.b.block(uSize,     0, pSize,     1);
    //
    // // step 1
    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solverA(A);
    // solverA.analyzePattern(A);
    // solverA.factorize(A);
    // Eigen::VectorXd u1 = solverA.solve(f);
    // // auto const Ainv = A.inverse();
    // // Eigen::VectorXd u1 = Ainv * f;
    //
    // // step 2
    // // auto const Sfull = M - B * Ainv * Bt;
    // Eigen::MatrixXd S = M - B * A.diagonal().asDiagonal().inverse() * Bt;
    // Eigen::PartialPivLU<Eigen::MatrixXd> solverS(S);
    // Eigen::VectorXd p = solverS.solve(g - B * u1);
    // // Eigen::VectorXd p = S.inverse() * (g - B * u1);
    //
    // // step 3
    // // Eigen::VectorXd u = Ainv * (f - Bt * p);
    // Eigen::VectorXd u = solverA.solve(f - Bt * p);
    //
    // sol.data.block(    0,     0, uSize,     1) = u;
    // sol.data.block(uSize,     0, pSize,     1) = p;
    //
    // // print
    // // std::cout << "det Sfull:  " << Sfull.determinant() << std::endl;
    // // std::cout << "det S: " << S.determinant() << std::endl;
    // // std::cout << "norm u: " << u.norm() << std::endl;
    // // std::cout << "norm p: " << p.norm() << std::endl;
    //

    // Eigen::saveMarket(builder.A, "mat.m");
    // Eigen::saveMarketVector(builder.b, "rhs.m");
    // Eigen::VectorXi pm = Eigen::VectorXi::Zero(uSize + pSize);
    // for (uint k=0; k<pSize; ++k)
    // {
    //   pm[uSize + k] = 1;
    // }
    // Eigen::saveMarketVector(pm, "pm.m");
    // abort();


    // SchurSolver solve(builder.A, solverParams);
    // auto [numIter, res] = solve(builder.b, sol.data);
    // std::cout << "iter: " << numIter << ", res: " << res << std::endl;
    solver.compute(builder.A);
    sol.data = solver.solve(builder.b);
  }

  void ic(std::function<FVec<dim> (Vec3 const &)> const & f)
  {
    interpolateAnalyticFunction(f, feSpaceVel, sol.data);
    velOld = sol.data;
  }

  void print(double const time = 0.0)
  {
    ioVel.print({sol}, time);
    p.data = sol.data.block(feSpaceVel.dof.size*dim, 0, feSpaceP.dof.size, 1);
    ioP.print({p}, time);
  }

  Mesh const & mesh;
  NSParameters parameters;
  FESpaceVel_T feSpaceVel;
  FESpaceP_T feSpaceP;
  BCList<FESpaceVel_T> bcsVel;
  BCList<FESpaceP_T> bcsP;
  Builder<StorageType::RowMajor> builder;
  Var sol;
  Var p;
  Vec velOld;
  AssemblyS2SProjection<FESpaceVel_T> assemblyRhs;
  AssemblyAdvection<FESpaceVel_T, FESpaceVel_T> assemblyAdvection;
  typename Builder<StorageType::RowMajor>::Mat_T matFixed;
  Vec rhsFixed;
  // std::vector<char> pMask;
  // SolverParams solverParams;
  SchurSolver solver;
  IOManager<FESpaceVel_T> ioVel;
  IOManager<FESpaceP_T> ioP;
};

template <typename Mesh>
struct NSSolverSplit2D
{
  static int constexpr dim = Mesh::Elem_T::dim;

  using Elem_T = typename Mesh::Elem_T;
  using FESpaceU_T = FESpace<Mesh,
                             typename FEType<Elem_T, 2>::RefFE_T,
                             typename FEType<Elem_T, 2>::RecommendedQR>;
  using FESpaceVel_T = FESpace<Mesh,
                               typename FEType<Elem_T, 2>::RefFE_T,
                               typename FEType<Elem_T, 2>::RecommendedQR, dim>;
  using FESpaceP_T = FESpace<Mesh,
                             typename FEType<Elem_T, 1>::RefFE_T,
                             typename FEType<Elem_T, 2>::RecommendedQR>;
  using BuilderVel_T = Builder<StorageType::RowMajor>;
  using BuilderP_T = Builder<StorageType::ClmMajor>;
  // using Backend = amgcl::backend::eigen<double>;
  // using Solver = amgcl::make_solver<
  //     amgcl::amg<
  //         Backend,
  //         amgcl::coarsening::smoothed_aggregation,
  //         amgcl::relaxation::spai0
  //         // amgcl::relaxation::ilu0
  //         // amgcl::relaxation::gauss_seidel
  //     >,
  //     amgcl::solver::bicgstab<Backend>
  //     // amgcl::solver::gmres<Backend>
  // >;
  using Solver = IterSolver;

  explicit NSSolverSplit2D(Mesh const & m, NSParameters const & par):
    mesh{m},
    parameters{par},
    feSpaceVel{m},
    feSpaceU{m},
    feSpaceP{m},
    bcsU{feSpaceU},
    bcsV{feSpaceU},
    bcsP{feSpaceP},
    builderUStar{feSpaceU.dof.size},
    builderVStar{feSpaceU.dof.size},
    builderP{feSpaceP.dof.size},
    builderU{feSpaceU.dof.size},
    builderV{feSpaceU.dof.size},
    uStar{"uStar", feSpaceU.dof.size},
    vStar{"vStar", feSpaceU.dof.size},
    p{"p", feSpaceP.dof.size},
    u{"u", feSpaceU.dof.size},
    v{"v", feSpaceU.dof.size},
    velStar{feSpaceU.dof.size * dim},
    pOld{feSpaceP.dof.size},
    dp{feSpaceP.dof.size},
    vel{feSpaceU.dof.size * dim},
    assemblyMassUStar{1. / par.dt, feSpaceU},
    assemblyMassVStar{1. / par.dt, feSpaceU},
    assemblyStiffnessUStar{par.nu, feSpaceU},
    assemblyStiffnessVStar{par.nu, feSpaceU},
    assemblyURhs{1. / par.dt, u.data, feSpaceU},
    assemblyVRhs{1. / par.dt, v.data, feSpaceU},
    assemblyAdvectionU{1.0, vel, feSpaceVel, feSpaceU},
    assemblyAdvectionV{1.0, vel, feSpaceVel, feSpaceU},
    assemblyPOldU{1.0, pOld, feSpaceU, feSpaceP, {0}},
    assemblyPOldV{1.0, pOld, feSpaceU, feSpaceP, {1}},
    assemblyDivVelStar{-1.0, velStar, feSpaceP, feSpaceVel},
    assemblyUStarRhs{1.0, uStar.data, feSpaceU},
    assemblyVStarRhs{1.0, vStar.data, feSpaceU},
    assemblyGradPRhsU{-par.dt, dp, feSpaceU, feSpaceP, {0}},
    assemblyGradPRhsV{-par.dt, dp, feSpaceU, feSpaceP, {1}},
    ioVel{feSpaceU, par.outputDir / "vel"},
    ioP{feSpaceP, par.outputDir / "p"}
  {}

  void init()
  {
    // TODO: assert that this comes after setting up bcs
    builderP.buildProblem(AssemblyStiffness{parameters.dt, feSpaceP}, bcsP);
    builderP.closeMatrix();
    solverP.analyzePattern(builderP.A);
    solverP.factorize(builderP.A);
    rhsFixedP = builderP.b;

    builderU.buildProblem(AssemblyMass{1.0, feSpaceU}, bcsU);
    builderU.closeMatrix();
    solverU.analyzePattern(builderU.A);
    solverU.factorize(builderU.A);
    rhsFixedVel[0] = builderU.b;
    builderV.buildProblem(AssemblyMass{1.0, feSpaceU}, bcsV);
    builderV.closeMatrix();
    solverV.analyzePattern(builderV.A);
    solverV.factorize(builderV.A);
    rhsFixedVel[1] = builderV.b;
  }

  void assemblyStepVelStar()
  {
    setComponent(vel, feSpaceVel, u.data, feSpaceU, 0);
    setComponent(vel, feSpaceVel, v.data, feSpaceU, 1);
    pOld += dp;
    builderUStar.clear();
    builderUStar.buildProblem(assemblyMassUStar, bcsU);
    builderUStar.buildProblem(assemblyAdvectionU, bcsU);
    builderUStar.buildProblem(assemblyStiffnessUStar, bcsU);
    builderUStar.buildProblem(assemblyURhs, bcsU);
    builderUStar.buildProblem(assemblyPOldU, bcsU);
    builderUStar.closeMatrix();
    builderVStar.clear();
    builderVStar.buildProblem(assemblyMassVStar, bcsV);
    builderVStar.buildProblem(assemblyAdvectionV, bcsV);
    builderVStar.buildProblem(assemblyStiffnessVStar, bcsV);
    builderVStar.buildProblem(assemblyVRhs, bcsV);
    builderVStar.buildProblem(assemblyPOldV, bcsV);
    builderVStar.closeMatrix();
  }

  void solveVelStar()
  {
    // Solver solveUStar(builderUStar.A);
    // auto const [numIterUStar, resUStar] = solveUStar(builderUStar.b, uStar.data);
    // std::cout << "ustar - iter: " << numIterUStar << ", res: " << resUStar << std::endl;
    // Solver solveVStar(builderVStar.A);
    // auto const [numIterVStar, resVStar] = solveVStar(builderVStar.b, vStar.data);
    // std::cout << "vstar - iter: " << numIterVStar << ", res: " << resVStar << std::endl;
    solverUStar.compute(builderUStar.A);
    uStar.data = solverUStar.solve(builderUStar.b);
    solverVStar.compute(builderVStar.A);
    vStar.data = solverVStar.solve(builderVStar.b);
  }

  void assemblyStepP()
  {
    setComponent(velStar, feSpaceVel, uStar.data, feSpaceU, 0);
    setComponent(velStar, feSpaceVel, vStar.data, feSpaceU, 1);
    builderP.clearRhs();
    builderP.buildProblem(assemblyDivVelStar, bcsP);
    builderP.b += rhsFixedP;
  }

  void solveP()
  {
    dp = solverP.solve(builderP.b);
    p.data += dp;
  }

  void assemblyStepVel()
  {
    builderU.clearRhs();
    builderU.buildProblem(assemblyUStarRhs, bcsU);
    builderU.buildProblem(assemblyGradPRhsU, bcsU);
    builderU.b += rhsFixedVel[0];
    builderV.clearRhs();
    builderV.buildProblem(assemblyVStarRhs, bcsV);
    builderV.buildProblem(assemblyGradPRhsV, bcsV);
    builderV.b += rhsFixedVel[1];
  }

  void solveVel()
  {
    u.data = solverU.solve(builderU.b);
    v.data = solverV.solve(builderV.b);
  }

  void ic(std::function<FVec<dim> (Vec3 const &)> const & f)
  {
    interpolateAnalyticFunction([&f](Vec3 const & p) { return f(p)[0]; }, feSpaceU, u.data);
    uStar.data = u.data;
    interpolateAnalyticFunction([&f](Vec3 const & p) { return f(p)[1]; }, feSpaceU, v.data);
    vStar.data = v.data;
    dp = Vec::Zero(feSpaceP.dof.size);
    pOld = Vec::Zero(feSpaceP.dof.size);
    p.data = pOld;
  }

  void print(double const time = 0.0)
  {
    ioVel.print({uStar, vStar, u, v}, time);
    Var pPrint{"p"};
    ioP.print({p}, time);
  }

  Mesh const & mesh;
  NSParameters parameters;
  FESpaceVel_T feSpaceVel;
  FESpaceU_T feSpaceU;
  FESpaceP_T feSpaceP;
  BCList<FESpaceU_T> bcsU;
  BCList<FESpaceU_T> bcsV;
  BCList<FESpaceP_T> bcsP;
  Builder<StorageType::RowMajor> builderUStar;
  Builder<StorageType::RowMajor> builderVStar;
  Builder<StorageType::ClmMajor> builderP;
  Builder<StorageType::ClmMajor> builderU;
  Builder<StorageType::ClmMajor> builderV;
  Solver solverUStar;
  Solver solverVStar;
  LUSolver solverP;
  LUSolver solverU;
  LUSolver solverV;
  Var uStar;
  Var vStar;
  Var p;
  Var u;
  Var v;
  Vec velStar;
  Vec pOld;
  Vec dp;
  Vec vel;
  AssemblyMass<FESpaceU_T> assemblyMassUStar;
  AssemblyMass<FESpaceU_T> assemblyMassVStar;
  AssemblyStiffness<FESpaceU_T> assemblyStiffnessUStar;
  AssemblyStiffness<FESpaceU_T> assemblyStiffnessVStar;
  AssemblyS2SProjection<FESpaceU_T> assemblyURhs;
  AssemblyS2SProjection<FESpaceU_T> assemblyVRhs;
  AssemblyAdvection<FESpaceU_T, FESpaceVel_T> assemblyAdvectionU;
  AssemblyAdvection<FESpaceU_T, FESpaceVel_T> assemblyAdvectionV;
  AssemblyGradRhs2<FESpaceU_T, FESpaceP_T> assemblyPOldU;
  AssemblyGradRhs2<FESpaceU_T, FESpaceP_T> assemblyPOldV;
  AssemblyDivRhs<FESpaceP_T, FESpaceVel_T> assemblyDivVelStar;
  AssemblyS2SProjection<FESpaceU_T, FESpaceU_T> assemblyUStarRhs;
  AssemblyS2SProjection<FESpaceU_T, FESpaceU_T> assemblyVStarRhs;
  AssemblyGradRhs<FESpaceU_T, FESpaceP_T> assemblyGradPRhsU;
  AssemblyGradRhs<FESpaceU_T, FESpaceP_T> assemblyGradPRhsV;
  array<typename Builder<StorageType::RowMajor>::Mat_T, dim> matFixedVelStar;
  array<Vec, dim> rhsFixedVelStar;
  Vec rhsFixedP;
  array<Vec, dim> rhsFixedVel;
  IOManager<FESpaceU_T> ioVel;
  IOManager<FESpaceP_T> ioP;
};
