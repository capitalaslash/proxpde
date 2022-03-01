#pragma once

#include "def.hpp"

#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "var.hpp"

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

struct Comp
{
  static const std::vector<short_T> u;
  static const std::vector<short_T> v;
  static const std::vector<short_T> w;
  static const std::vector<short_T> uv;
  static const std::vector<short_T> uw;
  static const std::vector<short_T> vw;
  static const std::vector<short_T> uvw;
};

struct NSParameters
{
  double const dt;
  double const nu;
  fs::path const outputDir = "output";
  // avoid padding warning
  // long const dummy = 0;
};

template <typename FESpaceVel, typename FESpaceP, typename BCSVel, typename BCSP>
struct NSSolverMonolithic
{
  using FESpaceVel_T = FESpaceVel;
  using FESpaceP_T = FESpaceP;
  using Elem_T = typename FESpaceVel_T::Mesh_T::Elem_T;
  static int constexpr dim = Elem_T::dim;

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

  explicit NSSolverMonolithic(
      FESpaceVel_T const & feVel,
      FESpaceP_T const & feP,
      BCSVel const & bcsVV,
      BCSP const & bcsPP,
      NSParameters const & par):
      feSpaceVel{feVel},
      feSpaceP{feP},
      bcsVel{bcsVV},
      bcsP{bcsPP},
      parameters{par},
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
    builder.buildLhs(
        std::tuple{
            AssemblyScalarMass{1. / parameters.dt, feSpaceVel},
            AssemblyTensorStiffness{parameters.nu, feSpaceVel}},
        bcsVel);
    builder.buildCoupling(AssemblyGrad{-1.0, feSpaceVel, feSpaceP}, bcsVel, bcsP);
    builder.buildCoupling(AssemblyDiv{-1.0, feSpaceP, feSpaceVel}, bcsP, bcsVel);
    builder.buildLhs(std::tuple{AssemblyScalarMass{0., feSpaceP}}, bcsP);
    builder.closeMatrix();
    matFixed = builder.A;
    rhsFixed = builder.b;
  }

  void assemblyStep()
  {
    velOld = sol.data;
    builder.clear();
    builder.buildRhs(std::tuple{assemblyRhs}, bcsVel);
    builder.buildLhs(std::tuple{assemblyAdvection}, bcsVel);
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

  void ic(std::function<FVec<dim>(Vec3 const &)> const & f)
  {
    interpolateAnalyticFunction(f, feSpaceVel, sol.data);
    velOld = sol.data;
  }

  void print(double const time = 0.0)
  {
    ioVel.print({sol}, time);
    p.data = sol.data.block(feSpaceVel.dof.size * dim, 0, feSpaceP.dof.size, 1);
    ioP.print({p}, time);
  }

  FESpaceVel_T const & feSpaceVel;
  FESpaceP_T const & feSpaceP;
  BCSVel const & bcsVel;
  BCSP const & bcsP;
  NSParameters parameters;
  Builder<StorageType::RowMajor> builder;
  Var sol;
  Var p;
  Vec velOld;
  AssemblyS2SProjection<FESpaceVel_T> assemblyRhs;
  AssemblyAdvection<FESpaceVel_T, FESpaceVel_T> assemblyAdvection;
  typename Builder<StorageType::RowMajor>::Mat_T matFixed;
  Vec rhsFixed;
  // std::vector<uint8_t> pMask;
  // SolverParams solverParams;
  SchurSolver solver;
  IOManager<FESpaceVel_T> ioVel;
  IOManager<FESpaceP_T> ioP;
};

template <
    typename FESpaceU,
    typename FESpaceP,
    typename BCSU,
    typename BCSV,
    typename BCSP>
struct NSSolverSplit2D
{
  using FESpaceU_T = FESpaceU;
  using FESpaceP_T = FESpaceP;
  static int constexpr dim = FESpaceU_T::Mesh_T::Elem_T::dim;
  using Elem_T = typename FESpaceU_T::Mesh_T::Elem_T;
  using FESpaceVel_T = FESpace<
      typename FESpaceU_T::Mesh_T,
      typename FESpaceU_T::RefFE_T,
      typename FESpaceU_T::QR_T,
      dim>;
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

  explicit NSSolverSplit2D(
      FESpaceU_T const & feU,
      FESpaceP_T const & feP,
      BCSU const & bcsUU,
      BCSV const & bcsVV,
      BCSP const & bcsPP,
      NSParameters const & par):
      feSpaceU{feU},
      feSpaceP{feP},
      feSpaceVel{*feU.mesh},
      bcsU{bcsUU},
      bcsV{bcsVV},
      bcsP{bcsPP},
      parameters{par},
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
      assemblyPOldU{1.0, pOld, feSpaceP, feSpaceU, {0}},
      assemblyPOldV{1.0, pOld, feSpaceP, feSpaceU, {1}},
      assemblyDivVelStar{-1.0, velStar, feSpaceVel, feSpaceP},
      assemblyUStarRhs{1.0, uStar.data, feSpaceU},
      assemblyVStarRhs{1.0, vStar.data, feSpaceU},
      assemblyGradPRhsU{-par.dt, dp, feSpaceP, feSpaceU, {0}},
      assemblyGradPRhsV{-par.dt, dp, feSpaceP, feSpaceU, {1}},
      ioVel{feSpaceU, par.outputDir / "vel"},
      ioP{feSpaceP, par.outputDir / "p"}
  {}

  void init()
  {
    // TODO: assert that this comes after setting up bcs
    builderP.buildLhs(std::tuple{AssemblyStiffness{parameters.dt, feSpaceP}}, bcsP);
    builderP.closeMatrix();
    solverP.analyzePattern(builderP.A);
    solverP.factorize(builderP.A);
    rhsFixedP = builderP.b;

    builderU.buildLhs(std::tuple{AssemblyScalarMass{1.0, feSpaceU}}, std::tuple{});
    builderU.closeMatrix();
    solverU.analyzePattern(builderU.A);
    solverU.factorize(builderU.A);
  }

  void assemblyStepVelStar()
  {
    setComponent(vel, feSpaceVel, u.data, feSpaceU, 0);
    setComponent(vel, feSpaceVel, v.data, feSpaceU, 1);
    pOld += dp;
    builderUStar.clear();
    builderUStar.buildLhs(
        std::tuple{assemblyMassUStar, assemblyAdvectionU, assemblyStiffnessUStar},
        bcsU);
    builderUStar.buildRhs(std::tuple{assemblyURhs, assemblyPOldU}, bcsU);
    builderUStar.closeMatrix();
    builderVStar.clear();
    builderVStar.buildLhs(
        std::tuple{assemblyMassVStar, assemblyAdvectionV, assemblyStiffnessVStar},
        bcsV);
    builderVStar.buildRhs(std::tuple{assemblyVRhs, assemblyPOldV}, bcsV);
    builderVStar.closeMatrix();
  }

  void solveVelStar()
  {
    // Solver solveUStar(builderUStar.A);
    // auto const [numIterUStar, resUStar] = solveUStar(builderUStar.b, uStar.data);
    // std::cout << "ustar - iter: " << numIterUStar << ", res: " << resUStar <<
    // std::endl; Solver solveVStar(builderVStar.A); auto const [numIterVStar, resVStar]
    // = solveVStar(builderVStar.b, vStar.data); std::cout << "vstar - iter: " <<
    // numIterVStar << ", res: " << resVStar << std::endl;
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
    builderP.buildRhs(std::tuple{assemblyDivVelStar}, bcsP);
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
    builderU.buildRhs(std::tuple{assemblyUStarRhs, assemblyGradPRhsU}, std::tuple{});
    builderV.clearRhs();
    builderV.buildRhs(std::tuple{assemblyVStarRhs, assemblyGradPRhsV}, std::tuple{});
  }

  void solveVel()
  {
    u.data = solverU.solve(builderU.b);
    v.data = solverU.solve(builderV.b);
  }

  void ic(std::function<FVec<dim>(Vec3 const &)> const & f)
  {
    interpolateAnalyticFunction(
        [&f](Vec3 const & p) { return f(p)[0]; }, feSpaceU, u.data);
    uStar.data = u.data;
    interpolateAnalyticFunction(
        [&f](Vec3 const & p) { return f(p)[1]; }, feSpaceU, v.data);
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

  FESpaceU_T const & feSpaceU;
  FESpaceP_T const & feSpaceP;
  FESpaceVel_T feSpaceVel;
  BCSU const & bcsU;
  BCSV const & bcsV;
  BCSP const & bcsP;
  NSParameters parameters;
  Builder<StorageType::RowMajor> builderUStar;
  Builder<StorageType::RowMajor> builderVStar;
  Builder<StorageType::ClmMajor> builderP;
  Builder<StorageType::ClmMajor> builderU;
  Builder<StorageType::ClmMajor> builderV;
  Solver solverUStar;
  Solver solverVStar;
  LUSolver solverP;
  LUSolver solverU;
  Var uStar;
  Var vStar;
  Var p;
  Var u;
  Var v;
  Vec velStar;
  Vec pOld;
  Vec dp;
  Vec vel;
  AssemblyScalarMass<FESpaceU_T> assemblyMassUStar;
  AssemblyScalarMass<FESpaceU_T> assemblyMassVStar;
  AssemblyStiffness<FESpaceU_T> assemblyStiffnessUStar;
  AssemblyStiffness<FESpaceU_T> assemblyStiffnessVStar;
  AssemblyS2SProjection<FESpaceU_T> assemblyURhs;
  AssemblyS2SProjection<FESpaceU_T> assemblyVRhs;
  AssemblyAdvection<FESpaceU_T, FESpaceVel_T> assemblyAdvectionU;
  AssemblyAdvection<FESpaceU_T, FESpaceVel_T> assemblyAdvectionV;
  AssemblyGradRhs2<FESpaceU_T, FESpaceP_T> assemblyPOldU;
  AssemblyGradRhs2<FESpaceU_T, FESpaceP_T> assemblyPOldV;
  AssemblyDivRhs<FESpaceP_T, FESpaceVel_T> assemblyDivVelStar;
  AssemblyS2SProjection<FESpaceU_T> assemblyUStarRhs;
  AssemblyS2SProjection<FESpaceU_T> assemblyVStarRhs;
  AssemblyGradRhs<FESpaceU_T, FESpaceP_T> assemblyGradPRhsU;
  AssemblyGradRhs<FESpaceU_T, FESpaceP_T> assemblyGradPRhsV;
  std::array<typename Builder<StorageType::RowMajor>::Mat_T, dim> matFixedVelStar;
  std::array<Vec, dim> rhsFixedVelStar;
  Vec rhsFixedP;
  std::array<Vec, dim> rhsFixedVel;
  IOManager<FESpaceU_T> ioVel;
  IOManager<FESpaceP_T> ioP;
};

template <typename FESpaceWSS, typename FESpaceVel>
void computeElemWSS(
    Vec & wss,
    FESpaceWSS const & feSpaceWSS,
    Vec const & vel,
    FESpaceVel const & feSpaceVel,
    std::unordered_set<marker_T> const & markers,
    double const nu)
{
  int constexpr dim = FESpaceVel::dim;
  wss = Vec::Zero(feSpaceWSS.dof.size);
  FEVar velFE{feSpaceVel};
  velFE.data = vel;

  for (auto & facet: feSpaceWSS.mesh->facetList)
  {
    if (markers.contains(facet.marker))
    {
      // std::cout << "facet " << facet.id << std::endl;
      Vec3 const normal = facet.normal();
      Vec3 const tangent = (facet.pts[1]->coord - facet.pts[0]->coord).normalized();
      auto elem = facet.facingElem[0].ptr;
      auto const id = feSpaceWSS.dof.getId(elem->id);
      velFE.reinit(*elem);
      for (uint q = 0; q < FESpaceVel::QR_T::numPts; ++q)
      {
        FMat<3, 3> tau = FMat<3, 3>::Zero();
        FMat<3, dim> const grad = velFE.evaluateGrad(q);
        tau.template block<dim, 3>(0, 0) = nu * grad.transpose();
        // minus sign due to the outwards normal
        wss[id] -= velFE.feSpace.curFE.JxW[q] * (tau * normal).dot(tangent);
      }
      wss[id] /= elem->volume();
    }
  }
}

template <typename FESpaceWSS, typename FESpaceVel>
void computeFEWSS(
    Vec & wss,
    FESpaceWSS const & feSpaceWSS,
    Vec const & vel,
    FESpaceVel const & feSpaceVel,
    std::unordered_set<marker_T> const & markers,
    double const nu)
{
  uint constexpr dim = FESpaceVel::dim;
  using FacetFEVel_T = typename FESpaceVel::RefFE_T::FacetFE_T;
  using FacetQR_T = SideQR_T<typename FESpaceVel::QR_T>;
  using FacetCurFEVel_T = CurFE<FacetFEVel_T, FacetQR_T>;
  using FacetFESpaceVel_T = FESpace<
      typename FESpaceVel::Mesh_T,
      typename FESpaceVel::RefFE_T,
      SideGaussQR<typename FESpaceVel::Mesh_T::Elem_T, FacetQR_T::numPts>,
      dim>;

  wss = Vec::Zero(feSpaceWSS.dof.size);
  FacetCurFEVel_T facetCurFEVel;

  Grad_T<FESpaceVel> feSpaceGrad{*feSpaceVel.mesh};
  Vec grad;
  computeGradient(grad, feSpaceGrad, vel, feSpaceVel);
  Grad_T<FacetFESpaceVel_T> facetFESpaceGrad{*feSpaceVel.mesh};
  FEVar gradFacet{facetFESpaceGrad};
  gradFacet.data = grad;

  for (auto & facet: feSpaceWSS.mesh->facetList)
  {
    if (markers.contains(facet.marker))
    {
      // std::cout << "facet " << facet.id << std::endl;
      Vec3 const normal = facet.normal();
      Vec3 const tangent = (facet.pts[1]->coord - facet.pts[0]->coord).normalized();
      facetCurFEVel.reinit(facet);
      auto elem = facet.facingElem[0].ptr;
      auto const side = facet.facingElem[0].side;
      auto const id = feSpaceWSS.dof.getId(elem->id);
      gradFacet.reinit(*elem);
      for (uint q = 0; q < FacetQR_T::numPts; ++q)
      {
        FMat<1, dim * dim> tauVec =
            nu * gradFacet.evaluate(side * FacetQR_T::numPts + q);
        FMat<3, 3> tau = FMat<3, 3>::Zero();
        for (uint i = 0; i < dim; ++i)
          for (uint j = 0; j < dim; ++j)
          {
            tau(i, j) = tauVec[j + FESpaceVel::dim * i];
          }
        // minus sign due to the outwards normal
        wss[id] -= facetCurFEVel.JxW[q] * (tau * normal).dot(tangent);
      }
      wss[id] /= facet.volume();
    }
  }
}
