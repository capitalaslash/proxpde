#include "def.hpp"

#include "fespace.hpp"
#include "feutils.hpp"
#include "geo.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "mesh_refine.hpp"
#include "multigrid.hpp"
#include "reffe.hpp"
#include "timer.hpp"
#include "var.hpp"

#include <unsupported/Eigen/SparseExtra>

using namespace proxpde;

template <
    typename Prol,
    typename Rest,
    typename Smoother = IterSolver,
    typename CoarseSolver = LUSolver>
struct MGSolver
{
  using Prol_T = Prol;
  using Rest_T = Rest;
  using Smoother_T = Smoother;
  using CoarseSolver_T = CoarseSolver;

  MGSolver() = delete;

  MGSolver(short_T const n):
      numLevels{n},
      As{static_cast<size_t>(n)},
      prols{static_cast<size_t>(n) - 1},
      rests{static_cast<size_t>(n) - 1},
      smoothers{static_cast<size_t>(n) - 1}
  {}

  MGSolver(
      short_T const n,
      std::vector<Mat<StorageType::RowMajor> const *> mats,
      std::vector<Prol_T const *> ps,
      std::vector<Rest_T const *> rs):
      numLevels{n},
      As{mats},
      prols{ps},
      rests{rs},
      smoothers{static_cast<size_t>(n) - 1}
  {
    init(mats);
  }

  void init(
      short_T const level,
      Mat<StorageType::RowMajor> const & mat,
      Prol_T const & p,
      Rest_T const & r)
  {
    assert(level > 0);
    As[level] = &mat;
    smoothers[level - 1].compute(mat);
    smoothers[level - 1].setMaxIterations(smoothIters);
    prols[level - 1] = &p;
    rests[level - 1] = &r;
  }

  void solve(Vec const & b, Vec & x) { step(b, x, numLevels - 1); }

  void step(Vec const & b, Vec & x, short_T const level)
  {
    // 1) relax with smoother
    Vec xTilda = smoothers[level - 1].solveWithGuess(b, x);
    std::cout << "pre-smoother n. iterations: " << smoothers[level - 1].iterations()
              << std::endl;
    std::cout << "pre-smoother rel. residual: " << smoothers[level - 1].error()
              << std::endl;

    // 2) compute residual
    Vec rFine = b - (*As[level]) * xTilda;

    // 3) restrict residual
    Vec rCoarse = rests[level - 1]->mat * rFine;

    // 4) coarse grid correction
    Vec error = Vec::Zero(As[level - 1]->rows());
    if (level == 1)
    {
      error = coarseSolver.solve(rCoarse);
      // std::cout << "coarse n. iterations: " << coarseSolver.iterations() <<
      // std::endl; std::cout << "coarse rel. residual: " << coarseSolver.error() <<
      // std::endl;
    }
    else
    {
      step(rCoarse, error, level - 1);
    }

    // 5) prolongate coarse grid correction
    xTilda += prols[level - 1]->mat * error;

    // 6) relax with smoother
    x = smoothers[level - 1].solveWithGuess(b, xTilda);
    std::cout << "post-smoother n. iterations: " << smoothers[level - 1].iterations()
              << std::endl;
    std::cout << "post-smoother rel. residual: " << smoothers[level - 1].error()
              << std::endl;
  }

  short_T numLevels;
  std::vector<Mat<StorageType::RowMajor> const *> As;
  std::vector<Prol_T const *> prols;
  std::vector<Rest_T const *> rests;
  std::vector<Smoother_T> smoothers;
  CoarseSolver_T coarseSolver;
  short_T smoothIters = 2;
};

template <typename BC, typename FEDest>
struct BCRecast
{
  // static_assert() compatibility of FESpaces
  using type = BCEss<FEDest>;
};

template <typename BC, typename FEDest>
using BCRecast_T = typename BCRecast<BC, FEDest>::type;

template <typename... BCs, typename FESpaceDest>
auto convertToHomogeneous(
    std::tuple<BCs...> const & bcsOrig, FESpaceDest const & feDest)
{
  // std::tuple<BCRecast_T<BCs, FESpaceDest>...> bcsDest;

  auto const zero = [](Vec3 const &) { return 0.; };
  static_for(
      bcsOrig,
      // bcsDest,
      [&feDest, &zero](uint const /*i*/, auto const & bcOrig)
      {
        auto bcDest = BCEss{feDest, bcOrig.marker};
        bcDest << zero;
      });
  // return bcsDest;
}

template <typename RefFE, typename Function>
int test(Function const & rhs, double const /*expectedNorm*/)
{
  using Elem_T = typename RefFE::GeoElem_T;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T =
      FESpace<Mesh_T, RefFE, typename LagrangeFE<Elem_T, 1>::RecommendedQR>;

  MilliTimer t;

  uint const nCoarse = 4;

  t.start("mesh");
  // std::unique_ptr<Mesh_T> meshCoarse{new Mesh_T};
  // // referenceMesh(*meshCoarse);
  // buildHyperCube(
  //     *meshCoarse,
  //     {0.0, 0.0, 0.0},
  //     {1.0, 1.0, 0.0},
  //     {nCoarse, nCoarse, 0},
  //     MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);
  // // readGMSH(
  // //     *meshCoarse,
  // //     "square_uns.msh",
  // //     MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);

  short_T constexpr numRefs = 3;

  std::array<std::unique_ptr<Mesh_T>, numRefs + 1> meshes;
  for (short_T ref = 0; ref < numRefs + 1; ++ref)
  {
    meshes[ref] = std::make_unique<Mesh_T>(Mesh_T{});
  }
  buildHyperCube(
      *meshes[0],
      {0.0, 0.0, 0.0},
      {1.0, 1.0, 0.0},
      {nCoarse, nCoarse, 0},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);
  t.stop();

  t.start("mesh refine");
  // std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // uniformRefine(*meshCoarse, *mesh);
  for (short_T ref = 0; ref < numRefs; ++ref)
  {
    uniformRefine(*meshes[ref], *meshes[ref + 1]);
  }
  // auto const & meshTop = meshes.back();
  t.stop();

  t.start("fespace");
  // FESpace_T feSpaceCoarse{*meshCoarse};
  // FESpace_T feSpace{*mesh};
  std::array<FESpace_T, numRefs + 1> feSpaces;
  for (short_T level = 0; level < numRefs + 1; ++level)
  {
    feSpaces[level].init(*meshes[level]);
  }
  auto const & feSpaceTop = feSpaces.back();
  t.stop();

  t.start("bcs");
  // auto bcLeft = BCEss{feSpaceTop, side::LEFT};
  // bcLeft << [](Vec3 const &) { return 0.; };
  // auto bcBottom = BCEss{feSpaceTop, side::BOTTOM};
  // bcBottom << [](Vec3 const &) { return 0.; };
  // auto const bcs = std::array{bcLeft, bcBottom};

  std::array<std::vector<BCEss<FESpace_T>>, numRefs + 1> bcLists;
  auto const zero = [](Vec3 const &) { return 0.; };
  for (short_T level = 0; level < numRefs + 1; ++level)
  {
    bcLists[level].push_back(BCEss{feSpaces[level], side::LEFT, zero});
    bcLists[level].push_back(BCEss{feSpaces[level], side::BOTTOM, zero});
  }
  auto const & bcListTop = bcLists.back();
  t.stop();

  t.start("mg operators");
  // Prolongator prol{feSpaceCoarse, feSpace};
  // Restrictor rest{feSpace, feSpaceCoarse};
  // rest.mat = prol.mat.transpose();
  using Prol_T = Prolongator<FESpace_T>;
  using Rest_T = Restrictor<FESpace_T, decltype(bcListTop)>;
  std::array<Prol_T, numRefs> prols;
  std::array<Rest_T, numRefs> rests;
  for (short_T level = 0; level < numRefs; ++level)
  {
    prols[level].init(feSpaces[level], feSpaces[level + 1]);
    rests[level].init(feSpaces[level + 1], feSpaces[level], bcLists[level]);
  }
  t.stop();

  t.start("assembly");
  std::array<Builder<StorageType::RowMajor>, numRefs + 1> builders;
  for (short_T level = 0; level < numRefs + 1; ++level)
  {
    builders[level].init(feSpaces[level].dof.size);
    builders[level].buildLhs(
        std::tuple{AssemblyStiffness{1.0, feSpaces[level]}}, bcLists[level]);
    builders[level].buildRhs(
        std::tuple{AssemblyAnalyticRhs{rhs, feSpaces[level]}}, bcLists[level]);
    builders[level].closeMatrix();
    // Mat<StorageType::RowMajor> Acoarse = rest.mat * builder.A * rest.mat.transpose();
  }
  auto const & builderTop = builders.back();
  t.stop();

  // Eigen::saveMarket(builderTop.A, "mat.mtx");
  // Eigen::saveMarket(builderTop.b, "b.mtx");

  // double const absToll = 1.e-14;
  double const relToll = 1.e-10;
  int const maxIter = 10;
  // int const smoothIter = 1;

  t.start("solve iter");
  Var sol{"u", feSpaceTop.dof.size};
  sol.data = Vec::Zero(feSpaceTop.dof.size);
  IterSolver solverTop;
  solverTop.setMaxIterations(maxIter);
  solverTop.compute(builderTop.A);

  // auto const bNorm = builderTop.b.norm();
  solverTop.setTolerance(relToll);
  sol.data = solverTop.solveWithGuess(builderTop.b, sol.data);
  std::cout << "solver n. iterations: " << solverTop.iterations() << std::endl;
  std::cout << "solver rel. residual: " << solverTop.error() << std::endl;
  t.stop();

  t.start("solve mg");
  MGSolver<Prol_T, Rest_T> mgSolver{numRefs + 1};
  for (short_T l = 1; l < numRefs + 1; ++l)
  {
    mgSolver.init(l, builders[l].A, prols[l - 1], rests[l - 1]);
  }
  mgSolver.As[0] = &builders[0].A;
  mgSolver.coarseSolver.compute(builders[0].A);

  Var solMG{"uMG", feSpaceTop.dof.size};
  solMG.data = Vec::Zero(feSpaceTop.dof.size);

  mgSolver.solve(builderTop.b, solMG.data);

  // std::array<IterSolver, numRefs + 1> smoothers;
  // for (short_T level = 1; level < numRefs + 1; ++level)
  // {
  //   smoothers[level].setMaxIterations(smoothIter);
  //   smoothers[level].compute(builders[level].A);
  // }
  // // bottom level should solve as much as possible, even use direct solver
  // smoothers[0].setMaxIterations(maxIter);
  // smoothers[0].setTolerance(1e-16);
  // smoothers[0].compute(builders[0].A);

  // std::array<Vec, numRefs + 1> residuals;
  // std::array<Vec, numRefs + 1> errors;
  // for (short_T level = 0; level < numRefs; ++level)
  // {
  //   // residuals[level] = Vec::Zero(feSpaces[level].dof.size);
  //   errors[level] = Vec::Zero(feSpaces[level].dof.size);
  // }
  // residuals.back() = builderTop.b - builderTop.A * solMG.data;
  // // set from initial condition
  // errors.back() = solMG.data;

  // for (short_T level = numRefs; level > 0; --level)
  // {
  //   Vec const e = smoothers[level].solveWithGuess(residuals[level], errors[level]);
  //   std::cout << "pre-smoother n. iterations: " << smoothers[level].iterations()
  //             << std::endl;
  //   std::cout << "pre-smoother residual:      " << smoothers[level].error()
  //             << std::endl;
  //   residuals[level - 1] = rests[level - 1].mat * builders[level].A * e;
  // }

  // // bottom
  // errors[0] = smoothers[0].solveWithGuess(residuals[0], errors[0]);
  // std::cout << "coarse n. iterations: " << smoothers[0].iterations() << std::endl;
  // std::cout << "coarse residual:      " << smoothers[0].error() << std::endl;

  // for (short_T level = 0; level < numRefs; ++level)
  // {
  //   errors[level + 1] = prols[level].mat * errors[level];
  //   Vec const res = builders[level + 1].A * errors[level + 1];
  //   errors[level + 1] = smoothers[level + 1].solveWithGuess(res, errors[level + 1]);
  //   std::cout << "post-smoother n. iterations: " << smoothers[level + 1].iterations()
  //             << std::endl;
  //   std::cout << "post-smoother residual:      " << smoothers[level + 1].error()
  //             << std::endl;
  // }

  // solMG.data += errors.back();

  // smoothers[2].setTolerance(relToll);
  // errors[2] = smoothers[2].solveWithGuess(residuals[2], errors[2]);
  // std::cout << "pre-smoother n. iterations: " << smoothers[2].iterations() <<
  // std::endl; std::cout << "pre-smoother rel. residual: " << smoothers[2].error() <<
  // std::endl;

  // residuals[2] = builders[2].A * errors[2];
  // residuals[1] = rests[1].mat * residuals[2];
  // errors[1] = smoothers[1].solveWithGuess(residuals[1], errors[1]);
  // std::cout << "pre-smoother n. iterations: " << smoothers[1].iterations() <<
  // std::endl; std::cout << "pre-smoother residual:      " << smoothers[1].error()
  // << std::endl; residuals[1] = builders[1].A * errors[1];

  // residuals[0] = rests[0].mat * residuals[1];
  // errors[0] = smoothers[0].solveWithGuess(residuals[0], errors[0]);
  // std::cout << "coarse n. iterations: " << smoothers[0].iterations() <<
  // std::endl; std::cout << "coarse residual:      " << smoothers[0].error() <<
  // std::endl;

  // errors[1] = prols[0].mat * errors[0];
  // residuals[1] = builders[1].A * errors[1];
  // errors[1] = smoothers[1].solveWithGuess(residuals[1], errors[1]);
  // std::cout << "post-smoother n. iterations: " << smoothers[1].iterations()
  //           << std::endl;
  // std::cout << "post-smoother residual:      " << smoothers[1].error() <<
  // std::endl;

  // // errors[2] = prols[1].mat * errors[1];

  // errors[2] = smoothers[2].solveWithGuess(builders[2].b, errors[2]);
  // std::cout << "post-smoother n. iterations: " << smoothers[2].iterations()
  //           << std::endl;
  // std::cout << "post-smoother residual:      " << smoothers[2].error() <<
  // std::endl;

  // solMG.data += errors[2];

  // auto & smootherTop = smoothers.back();
  // smootherTop.setTolerance(relToll);
  // Vec u = smootherTop.solveWithGuess(builderTop.b, Vec::Zero(builderTop.b.size()));
  // std::cout << "post-smoother n. iterations: " << smootherTop.iterations() <<
  // std::endl; std::cout << "post-smoother rel. residual: " << smootherTop.error() <<
  // std::endl;

  // Vec const r = builderTop.b - builderTop.A * errors[2];
  // auto const rNorm = r.norm();
  // smootherTop.setTolerance(absToll / rNorm);
  // u = smootherTop.solveWithGuess(builderTop.b, u);
  // std::cout << "post-smoother n. iterations: " << smootherTop.iterations() <<
  // std::endl; std::cout << "post-smoother rel. residual: " << smootherTop.error() <<
  // std::endl;

  // solMG.data = errors[2];

  t.stop();

  std::cout << "difference between the 2 solutions: " << (sol.data - solMG.data).norm()
            << std::endl;

  // Var exact{"exact"};
  // interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  // Var error{"e"};
  // error.data = sol.data - exact.data;

  t.start("output");
  IOManager io{feSpaceTop, "output_multigrid/sol"};
  io.print({sol, solMG});
  t.stop();

  t.print();

  // double norm = error.data.norm();
  // std::cout << "the norm of the error is " << std::setprecision(16) << norm
  //           << std::endl;
  // return checkError({norm}, {0.02049777877937642});
  return 0;
}

int main()
{
  std::bitset<4> tests;

  auto const fScalar = [](Vec3 const & p)
  {
    return 2.5 * M_PI * M_PI * std::sin(0.5 * M_PI * p(0)) *
           std::sin(1.5 * M_PI * p(1));
    ;
  };
  // auto const fVector = [](Vec3 const & p) { return Vec3{p[0], -p[1], 0.0}; };

  // tests[0] = test<RefTriangleP1>(fScalar, 1.762571852935e+00);
  // tests[1] = test<RefTriangleRT0>(fVector, 1.433720877840e+00);
  tests[2] = test<RefQuadQ1>(fScalar, 1.759929209669e+00);
  // tests[3] = test<RefQuadRT0>(fVector, 1.092906420717e+00);

  return tests.any();
}
