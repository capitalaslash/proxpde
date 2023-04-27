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

template <typename Mesh>
struct NSSolverMonolithic
{
  using Mesh_T = Mesh;
  using Elem_T = typename Mesh_T::Elem_T;
  static int constexpr dim = Elem_T::dim;
  using FESpaceVel_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 2>::RefFE_T,
      typename LagrangeFE<Elem_T, 2>::RecommendedQR,
      dim>;
  using FESpaceP_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 1>::RefFE_T,
      typename LagrangeFE<Elem_T, 2>::RecommendedQR>;

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

  explicit NSSolverMonolithic(Mesh_T const & mesh, ParameterDict const & c):
      feSpaceVel{mesh},
      feSpaceP{mesh, feSpaceVel.dof.size * dim},
      config{c},
      builder{feSpaceVel.dof.size * dim + feSpaceP.dof.size},
      sol{"vel", feSpaceVel.dof.size * dim + feSpaceP.dof.size},
      p{"p", feSpaceP.dof.size},
      velOld{feSpaceVel.dof.size * dim},
      assemblyRhs{1. / c["dt"].as<double>(), velOld, feSpaceVel},
      assemblyAdvection{1.0, velOld, feSpaceVel, feSpaceVel}
  // pMask(feSpaceVel.dof.size * dim + feSpaceP.dof.size, 0)
  {
    // // set up size and position of pressure dofs
    // std::fill(pMask.begin() + feSpaceVel.dof.size * dim, pMask.end(), 1);
    // solverParams.put("precond.pmask_size", pMask.size());
    // solverParams.put("precond.pmask", static_cast<void *>(pMask.data()));
    // // solverParams.put("precond.approx_schur", false);
    // solverParams.put("precond.psolver.precond.direct_coarse", false);

    config.validate({"nu", "dt", "output_dir"});
    auto const outDir = fs::path{config["output_dir"].as<std::string>()};
    ioVel.init(feSpaceVel, outDir / "vel");
    ioP.init(feSpaceP, outDir / "p");
  }

  template <typename BCsVel, typename BCsP>
  void init(BCsVel const & bcsVel, BCsP const & bcsP)
  {
    // TODO: assert that this comes after setting up bcs
    builder.buildLhs(
        std::tuple{
            AssemblyScalarMass{1. / config["dt"].as<double>(), feSpaceVel},
            AssemblyTensorStiffness{config["nu"].as<double>(), feSpaceVel}},
        bcsVel);
    builder.buildCoupling(AssemblyGrad{-1.0, feSpaceVel, feSpaceP}, bcsVel, bcsP);
    builder.buildCoupling(AssemblyDiv{-1.0, feSpaceP, feSpaceVel}, bcsP, bcsVel);
    builder.buildLhs(std::tuple{AssemblyDummy{feSpaceP}}, bcsP);
    builder.closeMatrix();
    matFixed = builder.A;
    rhsFixed = builder.b;
  }

  template <typename BCsVel>
  void assemblyStep(BCsVel const & bcsVel)
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
    // Eigen::VectorXd f  = builder.b.segment(    0, uSize);
    // Eigen::VectorXd g  = builder.b.segment(uSize, pSize);
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
    // sol.data.segent(0, uSize) = u;
    // sol.data.segment(uSize, pSize) = p;
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
    p.data = sol.data.tail(feSpaceP.dof.size);
    ioP.print({p}, time);
  }

  FESpaceVel_T const feSpaceVel;
  FESpaceP_T const feSpaceP;
  ParameterDict const config;
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

enum class VelStarSolverType
{
  MONOLITHIC,
  SEGREGATED,
};

template <typename Mesh>
struct NSSolverSplit
{
  using Mesh_T = Mesh;
  using Elem_T = typename Mesh_T::Elem_T;
  static int constexpr dim = Elem_T::dim;
  using FESpaceVel_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 2>::RefFE_T,
      typename LagrangeFE<Elem_T, 2>::RecommendedQR,
      dim>;
  using FESpaceP_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 1>::RefFE_T,
      typename LagrangeFE<Elem_T, 2>::RecommendedQR>;
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

  explicit NSSolverSplit(Mesh_T const & mesh, ParameterDict const & c):
      feSpaceVel{mesh},
      feSpaceP{mesh},
      config{c},
      builderVelStar{feSpaceVel.dof.size * dim},
      builderP{feSpaceP.dof.size},
      builderVel{feSpaceVel.dof.size * dim},
      velStar{"velStar", feSpaceVel.dof.size * dim},
      velStarOld{feSpaceVel.dof.size * dim},
      p{"p", feSpaceP.dof.size},
      vel{"vel", feSpaceVel.dof.size * dim},
      pOld{feSpaceP.dof.size},
      dp{feSpaceP.dof.size},
      assemblyMassVelStar{1. / c["dt"].as<double>(), feSpaceVel},
      assemblyStiffness{c["nu"].as<double>(), feSpaceVel},
      assemblyVelRhs{1. / c["dt"].as<double>(), vel.data, feSpaceVel},
      assemblyAdvection{1.0, vel.data, feSpaceVel, feSpaceVel},
      assemblyPOld{1.0, pOld, feSpaceP, feSpaceVel},
      assemblyDivVelStar{-1.0, velStar.data, feSpaceVel, feSpaceP},
      assemblyVelStarRhs{1.0, velStar.data, feSpaceVel},
      assemblyGradPRhs{-c["dt"].as<double>(), dp, feSpaceP, feSpaceVel}
  {
    config.validate({"nu", "dt", "output_dir"});
    auto const outDir = fs::path{config["output_dir"].as<std::string>()};
    ioVel.init(feSpaceVel, outDir / "vel");
    ioP.init(feSpaceP, outDir / "p");
  }

  template <typename BCsVel, typename BCsP>
  void init(BCsVel const & bcsVel, BCsP const & bcsP)
  {
    // TODO: assert that this comes after setting up bcs
    builderP.buildLhs(
        std::tuple{AssemblyStiffness{config["dt"].as<double>(), feSpaceP}}, bcsP);
    builderP.closeMatrix();
    solverP.analyzePattern(builderP.A);
    solverP.factorize(builderP.A);
    rhsFixedP = builderP.b;

    builderVel.buildLhs(std::tuple{AssemblyScalarMass{1.0, feSpaceVel}}, bcsVel);
    builderVel.closeMatrix();
    solverVel.analyzePattern(builderVel.A);
    solverVel.factorize(builderVel.A);
  }

  template <typename BCsVel>
  void assemblyStep(BCsVel const & bcsVel)
  {
    assemblyStepVelStar(bcsVel);
  }

  template <VelStarSolverType SolverType, typename BCsP>
  void solve(BCsP const & bcsP)
  {
    solveVelStar<SolverType>();
    assemblyStepP(bcsP);
    solveP();
    assemblyStepVel();
    solveVel();
  }

  template <typename BCsVel>
  void assemblyStepVelStar(BCsVel const & bcsVel)
  {
    velStarOld = velStar.data;
    pOld += dp;
    builderVelStar.clear();
    builderVelStar.buildLhs(
        std::tuple{assemblyMassVelStar, assemblyAdvection, assemblyStiffness}, bcsVel);
    builderVelStar.buildRhs(std::tuple{assemblyVelRhs, assemblyPOld}, bcsVel);
    builderVelStar.closeMatrix();
  }

  template <VelStarSolverType SolverType>
  void solveVelStar()
  {
    if constexpr (SolverType == VelStarSolverType::MONOLITHIC)
    {
      // solve all components of the velocity monolithically
      solverVelStar.compute(builderVelStar.A);
      velStar.data = solverVelStar.solve(builderVelStar.b);
    }
    else if constexpr (SolverType == VelStarSolverType::SEGREGATED)
    {
      // solve each component separately
      auto const dofU = feSpaceVel.dof.size;
      for (short_T d = 0; d < dim; ++d)
      {
        solverVelStar.compute(builderVelStar.A.block(d * dofU, d * dofU, dofU, dofU));
        auto b = builderVelStar.b.segment(d * dofU, dofU);

        // subtract all the other variables considering them fixed
        for (short_T od = 0; od < dim; ++od)
        {
          if (od != d)
          {
            b -= builderVelStar.A.block(d * dofU, od * dofU, dofU, dofU) *
                 velStarOld.segment(od * dofU, dofU);
          }
        }

        velStar.data.segment(d * dofU, dofU) = solverVelStar.solve(b);
      }
    }
  }

  template <typename BCsP>
  void assemblyStepP(BCsP const & bcsP)
  {
    builderP.clearRhs();
    builderP.buildRhs(std::tuple{assemblyDivVelStar}, bcsP);
    builderP.b += rhsFixedP;
  }

  void solveP()
  {
    dp = solverP.solve(builderP.b);
    p.data += dp;
  }

  template <typename BCsVel>
  void assemblyStepVel(BCsVel const & bcsVel)
  {
    builderVel.clearRhs();
    builderVel.buildRhs(std::tuple{assemblyVelStarRhs, assemblyGradPRhs}, bcsVel);
  }

  void solveVel()
  {
    // this can also be split in components, but since it is pre-calculated once it is
    // not so important
    vel.data = solverVel.solve(builderVel.b);
  }

  void ic(std::function<FVec<dim>(Vec3 const &)> const & f)
  {
    interpolateAnalyticFunction(f, feSpaceVel, vel.data);
    velStar.data = vel.data;
    velStarOld = vel.data;
    dp = Vec::Zero(feSpaceP.dof.size);
    pOld = Vec::Zero(feSpaceP.dof.size);
    p.data = pOld;
  }

  void print(double const time = 0.0)
  {
    ioVel.print({velStar, vel}, time);
    Var pPrint{"p"};
    ioP.print({p}, time);
  }

  FESpaceVel_T const feSpaceVel;
  FESpaceP_T const feSpaceP;
  ParameterDict const config;
  Builder<StorageType::RowMajor> builderVelStar;
  Builder<StorageType::ClmMajor> builderP;
  Builder<StorageType::ClmMajor> builderVel;
  Solver solverVelStar;
  LUSolver solverP;
  LUSolver solverVel;
  Var velStar;
  Vec velStarOld;
  Var p;
  Var vel;
  Vec pOld;
  Vec dp;
  AssemblyScalarMass<FESpaceVel_T> assemblyMassVelStar;
  AssemblyStiffness<FESpaceVel_T> assemblyStiffness;
  AssemblyS2SProjection<FESpaceVel_T> assemblyVelRhs;
  AssemblyAdvection<FESpaceVel_T, FESpaceVel_T> assemblyAdvection;
  AssemblyGradRhs2<FESpaceVel_T, FESpaceP_T> assemblyPOld;
  AssemblyDivRhs<FESpaceP_T, FESpaceVel_T> assemblyDivVelStar;
  AssemblyS2SProjection<FESpaceVel_T> assemblyVelStarRhs;
  AssemblyGradRhs<FESpaceVel_T, FESpaceP_T> assemblyGradPRhs;
  std::array<typename Builder<StorageType::RowMajor>::Mat_T, dim> matFixedVelStar;
  std::array<Vec, dim> rhsFixedVelStar;
  Vec rhsFixedP;
  std::array<Vec, dim> rhsFixedVel;
  IOManager<FESpaceVel_T> ioVel;
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
        wss[id] -= velFE.feSpace->curFE.JxW[q] * (tau * normal).dot(tangent);
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
  using FacetCurFEVel_T = typename CurFETraits<FacetFEVel_T, FacetQR_T>::type;
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
