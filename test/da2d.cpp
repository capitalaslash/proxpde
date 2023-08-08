#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "fv.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

namespace proxpde
{

template <int N>
FVec<N>
addSUPG(double const h, double const eps, Vec3 const & vel, FMat<N, 3> const & dphi)
{
  double const velNorm = vel.norm();
  double const peclet = 0.5 * velNorm * h / eps;
  double const aOpt = (1. / std::tanh(peclet) - 1. / peclet);
  return 0.5 * aOpt * h * (dphi * vel) / velNorm;
}

template <typename FESpace, typename FESpaceVel>
struct AssemblyMassSUPG: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceVel_T = FESpaceVel;
  using Super_T = Diagonal<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;

  AssemblyMassSUPG(
      double const c,
      FEVar<FESpaceVel_T> & u,
      double const e,
      FESpace_T const & fe,
      bool supg,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>{fe, cmp},
      coef{c},
      vel{&u},
      eps{e},
      enableSUPG{supg}
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    CurFE_T const & curFE = this->feSpace->curFE;
    double const h = curFE.elem->hMin();
    vel->reinit(*curFE.elem);

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      Vec3 const velQPoint = promote<3>(vel->evaluate(q));
      auto const phiiSUPG =
          curFE.phi[q] + enableSUPG * addSUPG(h, eps, velQPoint, curFE.dphi[q]);

      Ke += coef * this->feSpace->curFE.JxW[q] * phiiSUPG *
            this->feSpace->curFE.phi[q].transpose();
    }
  }

  double coef = 1.0;
  FEVar<FESpaceVel_T> * vel;
  double const eps;
  bool const enableSUPG;
};

template <typename FESpace, typename FESpaceVel>
struct AssemblyAdvectionSUPG: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using FESpaceVel_T = FESpaceVel;
  using Super_T = Diagonal<FESpace_T>;
  using LMat_T = typename Super_T::LMat_T;

  AssemblyAdvectionSUPG(
      double const c,
      FEVar<FESpaceVel_T> & u,
      double const e,
      FESpace_T const & fe,
      bool supg,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>{fe, cmp},
      coef{c},
      vel{&u},
      eps{e},
      enableSUPG{supg}
  {}

  void build(LMat_T & Ke) const override
  {
    using CurFE_T = typename FESpace_T::CurFE_T;

    CurFE_T const & curFE = this->feSpace->curFE;
    double const h = curFE.elem->hMin();
    vel->reinit(*curFE.elem);

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      auto const velQPoint = promote<3>(vel->evaluate(q));
      auto const phiiSUPG =
          curFE.phi[q] + enableSUPG * addSUPG(h, eps, velQPoint, curFE.dphi[q]);

      Ke += coef * this->feSpace->curFE.JxW[q] * phiiSUPG *
            (this->feSpace->curFE.dphi[q] * velQPoint).transpose();
    }
  }

  double coef = 1.0;
  FEVar<FESpaceVel_T> * vel;
  double const eps;
  bool const enableSUPG;
};

template <typename ElemType>
class DAEqn
{
public:
  using Elem_T = ElemType;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 1>::RefFE_T,
      typename LagrangeFE<Elem_T, 1>::RecommendedQR>;
  using FESpaceVel_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 1>::RefFE_T,
      typename LagrangeFE<Elem_T, 1>::RecommendedQR,
      Elem_T::dim>;

  DAEqn() = default;
  DAEqn(double const t, bool const supg): theta(t), enableSUPG{supg} {}
  void run();

  // theta = 0 explicit Euler
  // theta = 1 implicit Euler
  // theta = 0.5 Crank-Nicholson
  double const theta = 0.5;
  bool const enableSUPG = true;

  FEVar<FESpace_T> c;

private:
  void setupSystem();
  void solveTimeStep();
  void outputResults();

  MilliTimer t;
  std::unique_ptr<Mesh_T> mesh;
  FESpace_T feSpaceP1;
  std::vector<BCEss<FESpace_T>> bcs;
  Builder<StorageType::ClmMajor> builder;
  Mat<StorageType::RowMajor> massMatrix;
  Mat<StorageType::RowMajor> diffMatrix;
  Mat<StorageType::RowMajor> advMatrix;
  Mat<StorageType::ClmMajor> systemMatrix;
  Vec rhs;
  LUSolver solver;
  IOManager<FESpace_T> io;

  double time = 0.0;
  double const dt = 0.005;
  uint const ntime = 400;
  uint const printStep = 10U;

  Vec cOld;

  scalarFun_T const ic = [](Vec3 const & p)
  { return (p - Vec3(0.5, 0.75, 0.0)).squaredNorm() < 0.15 * 0.15 ? 1.0 : 0.0; };

  // circle
  Fun<Elem_T::dim, 3> const velocityField = [](Vec3 const & p) {
    return Vec2{M_PI * (p(1) - 0.5), -M_PI * (p(0) - 0.5)};
  };

  // // linear
  // Fun<Elem_T::dim, 3> const velocityField = [](Vec3 const & ) {
  //   return Vec2{0.0, -0.25};
  // };

  double const eps = 1.e-3;
};

template <typename ElemType>
void DAEqn<ElemType>::setupSystem()
{
  t.start("fespace");
  feSpaceP1.init(*mesh);
  t.stop();

  c.init("conc", feSpaceP1);
  cOld.resize(feSpaceP1.dof.size);

  t.start("bcs");
  auto const zero = [](Vec3 const &) { return 0.0; };
  bcs.reserve(AllSides<Elem_T::dim>::values.size());
  for (auto const s: AllSides<Elem_T::dim>::values)
  {
    auto bc = BCEss{feSpaceP1, s};
    bc << zero;
    bcs.push_back(bc);
  }
  t.stop();

  t.start("velocity");
  FESpaceVel_T feSpaceVel{*mesh};
  FEVar vel{"velocity", feSpaceVel};
  interpolateAnalyticFunction(velocityField, *vel.feSpace, vel.data);

  double const cfl = computeMaxCFL(feSpaceVel, vel.data, dt);
  fmt::print("cfl max: {}\n", cfl);
  t.stop();

  t.start("io");
  io.init(feSpaceP1, "output_da2d/sol");
  t.stop();

  t.start("build");
  Builder<StorageType::RowMajor> builderMatrix{feSpaceP1.dof.size};

  // AssemblyMass timeDer{1.0 / dt, feSpaceP1};
  AssemblyMassSUPG<FESpace_T, FESpaceVel_T> timeDer{
      1.0 / dt, vel, eps, feSpaceP1, enableSUPG};
  builderMatrix.buildLhs(std::tuple{timeDer}, bcs);
  builderMatrix.closeMatrix();
  massMatrix = builderMatrix.A;

  // with linear fe the second derivative is null, no term in SUPG
  AssemblyStiffness diffusion{eps, feSpaceP1};
  builderMatrix.clearLhs();
  builderMatrix.buildLhs(std::tuple{diffusion}, bcs);
  builderMatrix.closeMatrix();
  diffMatrix = builderMatrix.A;

  builderMatrix.clearLhs();
  // AssemblyAdvection advection{vel, feSpaceP1};
  AssemblyAdvectionSUPG<FESpace_T, FESpaceVel_T> advection{
      1.0, vel, eps, feSpaceP1, enableSUPG};
  builderMatrix.buildLhs(std::tuple{advection}, bcs);
  builderMatrix.closeMatrix();
  advMatrix = builderMatrix.A;

  systemMatrix = massMatrix + theta * (diffMatrix + advMatrix);
  t.stop();

  t.start("solver setup");
  solver.analyzePattern(systemMatrix);
  solver.factorize(systemMatrix);
  t.stop();
}

template <typename ElemType>
void DAEqn<ElemType>::solveTimeStep()
{
  t.start("solve");
  c.data = solver.solve(rhs);
  fmt::print("c max value: {:8.6f}\n", c.data.maxCoeff());
  t.stop();
}

template <typename ElemType>
void DAEqn<ElemType>::outputResults()
{
  t.start("print");
  io.print(std::tuple{c}, time);
  t.stop();
}

template <typename ElemType>
void DAEqn<ElemType>::run()
{
  t.start("mesh");
  mesh.reset(new Mesh_T);
  buildHyperCube(*mesh, {0., 0., 0.}, {1.0, 1.0, 0.0}, {32, 32, 0});
  t.stop();

  setupSystem();

  t.start("init");
  // FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, MiniQR<Elem_T, 10>>
  //     feSpaceIC{*mesh};
  // FESpace<Mesh_T, typename LagrangeFE<Elem_T, 0>::RefFE_T, typename
  // FESpaceP1_T::QR_T>
  //     feSpaceP0{*mesh};
  // Vec cCell(feSpaceIC.dof.size);
  // integrateAnalyticFunction(ic, feSpaceIC, cCell);
  // L2Projector projector{feSpaceP1, feSpaceP0};
  // projector.setRhs(cCell);
  // cOld = projector.apply();

  interpolateAnalyticFunction(ic, feSpaceP1, cOld);

  c.data = cOld;
  t.stop();

  outputResults();

  AssemblyProjection timeDerOld{1.0 / dt, cOld, feSpaceP1};
  for (uint itime = 0; itime < ntime; itime++)
  {
    time += dt;
    fmt::print("solving timestep {:6d}\n", itime + 1);

    t.start("update");
    cOld = c.data;
    // builder.clearRhs();
    // builder.buildRhs(std::tuple{timeDerOld}, bcs);
    // rhs = builder.b;

    // rhs = (massMatrix - (1.0 - theta) * (diffMatrix + advMatrix)) * cOld;
    // multiplying and adding vectors is 10x faster than adding sparse matrices in dbg
    // with n = 1000
    rhs = massMatrix * cOld;
    rhs -= (1.0 - theta) * (diffMatrix * cOld);
    rhs -= (1.0 - theta) * (advMatrix * cOld);
    t.stop();

    solveTimeStep();

    if ((itime + 1) % printStep == 0)
      outputResults();
  }

  t.print();
}

} // namespace proxpde

int main(/*int argc, char * argv[]*/)
{
  using namespace proxpde;

  std::bitset<8> tests;
  uint testCounter = 0U;
  auto const expectedValues = std::vector{{
      7.997857020567e-01,
      7.937748763851e-01,
      9.496437540798e-01,
      9.492604651605e-01,
      7.904235782602e-01,
      7.834533565465e-01,
      9.385255164133e-01,
      9.433922227426e-01,
  }};

  for (auto const theta: {1.0, 0.5})
    for (auto const enableSUPG: {false, true})
    {
      DAEqn<Triangle> daEqn{theta, enableSUPG};
      daEqn.run();
      tests[testCounter] =
          checkError({daEqn.c.data.maxCoeff()}, {expectedValues[testCounter]});
      testCounter++;
    }

  for (auto const theta: {1.0, 0.5})
    for (auto const enableSUPG: {false, true})
    {
      DAEqn<Quad> daEqn{theta, enableSUPG};
      daEqn.run();
      tests[testCounter] =
          checkError({daEqn.c.data.maxCoeff()}, {expectedValues[testCounter]});
      testCounter++;
    }

  fmt::print("tests: {}\n", tests.to_string());
  return tests.any();
}
