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

using namespace proxpde;

template <typename ElemType>
class DAEqn
{
public:
  using Elem_T = ElemType;
  using Mesh_T = Mesh<Elem_T>;
  using FESpaceP1_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 1>::RefFE_T,
      typename LagrangeFE<Elem_T, 1>::RecommendedQR>;
  using FESpaceVel_T = FESpace<
      Mesh_T,
      typename LagrangeFE<Elem_T, 1>::RefFE_T,
      typename LagrangeFE<Elem_T, 1>::RecommendedQR,
      Elem_T::dim>;

  DAEqn() = default;
  DAEqn(double const t): theta(t) {}
  void run();

  // theta = 0 explicit Euler
  // theta = 1 implicit Euler
  // theta = 0.5 Crank-Nicholson
  double const theta = 0.5;

  FEVar<FESpaceP1_T> c;

private:
  void setupSystem();
  void solveTimeStep();
  void outputResults();

  MilliTimer t;
  std::unique_ptr<Mesh_T> mesh;
  FESpaceP1_T feSpaceP1;
  std::vector<BCEss<FESpaceP1_T>> bcs;
  Builder<StorageType::ClmMajor> builder;
  Mat<StorageType::RowMajor> massMatrix;
  Mat<StorageType::RowMajor> diffMatrix;
  Mat<StorageType::RowMajor> advMatrix;
  Mat<StorageType::ClmMajor> systemMatrix;
  Vec rhs;
  LUSolver solver;
  IOManager<FESpaceP1_T> io;

  double time = 0.0;
  double const dt = 0.005;
  uint const ntime = 400;
  uint const printStep = 10U;

  Vec cOld;

  scalarFun_T const ic = [](Vec3 const & p)
  { return (p - Vec3(0.5, 0.75, 0.0)).squaredNorm() < 0.15 * 0.15 ? 1.0 : 0.0; };

  Fun<Elem_T::dim, 3> const velocityField = [](Vec3 const & p) {
    return Vec2{M_PI * (p(1) - 0.5), -M_PI * (p(0) - 0.5)};
  };

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

  AssemblyMass timeDer{1.0 / dt, feSpaceP1};
  builderMatrix.buildLhs(std::tuple{timeDer}, bcs);
  builderMatrix.closeMatrix();
  massMatrix = builderMatrix.A;

  AssemblyStiffness diffusion{eps, feSpaceP1};
  builderMatrix.clearLhs();
  builderMatrix.buildLhs(std::tuple{diffusion}, bcs);
  builderMatrix.closeMatrix();
  diffMatrix = builderMatrix.A;

  builderMatrix.clearLhs();
  AssemblyAdvection advection{vel, feSpaceP1};
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

int main(/*int argc, char * argv[]*/)
{
  std::bitset<4> tests;
  {
    DAEqn<Triangle> daEqn{/*theta = */ 1.0};
    daEqn.run();
    tests[0] = checkError({daEqn.c.data.maxCoeff()}, {7.997857020567e-01});
  }
  {
    DAEqn<Triangle> daEqn{/*theta = */ 0.5};
    daEqn.run();
    tests[1] = checkError({daEqn.c.data.maxCoeff()}, {9.496437540798e-01});
  }
  {
    DAEqn<Quad> daEqn{/*theta = */ 1.0};
    daEqn.run();
    tests[2] = checkError({daEqn.c.data.maxCoeff()}, {7.904235782602e-01});
  }
  {
    DAEqn<Quad> daEqn{/*theta = */ 0.5};
    daEqn.run();
    tests[3] = checkError({daEqn.c.data.maxCoeff()}, {9.385255164133e-01});
  }

  fmt::print("tests: {}\n", tests.to_string());
  return tests.any();
}
