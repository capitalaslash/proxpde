#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(int argc, char * argv[])
{
  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  scalarFun_T const qiii = [](Vec3 const &) { return 0.; };
  const double k = 1.0;
  // double const qiiLeft = 1.;
  double const tLeft = 1.;
  // double const qiiRight = 0.; // equivalent to no bc
  double const tRight = 0.0;
  double const dt = 0.1;
  uint const ntime = 50U;

  scalarFun_T const exactSol = [](Vec3 const & p) { return 1. - p[0]; };

  MilliTimer t;

  uint const numElems = (argc < 2) ? 20U : std::stoi(argv[1]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  t.start("mesh build");
  buildHyperCube(*mesh, origin, length, {{numElems, 0, 0}});
  t.stop();

  t.start("fespace");
  FESpace_T feSpace{*mesh};
  t.stop();

  t.start("bcs");
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [tLeft](Vec3 const &) { return tLeft; };
  auto bcRight = BCEss{feSpace, side::RIGHT};
  bcRight << [tRight](Vec3 const &) { return tRight; };
  auto const bcs = std::tuple{bcLeft, bcRight};
  t.stop();

  t.start("fe build");
  Vec tempOld = Vec::Zero(feSpace.dof.size);
  Vec tempOldInc = Vec::Zero(feSpace.dof.size);
  Vec res = Vec::Zero(feSpace.dof.size);
  AssemblyMass mass{1. / dt, feSpace};
  AssemblyStiffness diffusion{k, feSpace};
  AssemblyAnalyticRhs heatGeneration{qiii, feSpace};
  AssemblyProjection massOld{1. / dt, tempOld, feSpace};
  AssemblyProjection massOldInc{1. / dt, tempOldInc, feSpace};
  AssemblyProjection residual{-1., res, feSpace};
  AssemblyStiffnessRhs diffRhs{-1., tempOldInc, feSpace, feSpace};
  // // using an interpolated rhs makes its quality independent of the chosen qr
  // Vec rhsProj;
  // interpolateAnalyticFunction(rotatedRhs, feSpace, rhsProj);
  // AssemblyProjection f(1.0, rhsProj, feSpace);

  Builder builder{feSpace.dof.size};
  Builder builderInc{feSpace.dof.size};
  t.stop();

  Var temp{"temp"};
  Var dTemp{"dTemp"};
  Var tempInc{"tempInc"};
  // ic
  temp.data = Vec::Zero(feSpace.dof.size);
  dTemp.data = Vec::Zero(feSpace.dof.size);
  tempInc.data = Vec::Zero(feSpace.dof.size);
  tempOldInc = Vec::Zero(feSpace.dof.size);
  // the ic must be compatible with the bcs for the incremental solution!
  tempOldInc[0] = 1.0;

  LUSolver solver;
  IOManager io{feSpace, "output_heat1d/temp"};
  double time = 0.;
  io.print({temp, tempInc}, time);
  for (uint itime = 0; itime < ntime; ++itime)
  {
    time += dt;
    std::cout << "time = " << time << std::endl;

    tempOld = temp.data;
    builder.buildLhs(std::tuple{mass, diffusion}, bcs);
    builder.closeMatrix();
    builder.buildRhs(std::tuple{massOld, heatGeneration}, bcs);

    tempOldInc += dTemp.data;
    builderInc.buildLhs(std::tuple{mass, diffusion}, bcs);
    builderInc.closeMatrix();
    Mat<StorageType::ClmMajor> K{feSpace.dof.size, feSpace.dof.size};
    // K = builderInc.A;
    // Vec const Krhs = builderInc.b;
    // res = K * tempOld;
    builderInc.buildRhs(
        std::tuple{
            heatGeneration,
            // massOld,
            // residual,
            diffRhs,
        },
        bcs);
    // builderInc.b -= Krhs;

    t.start("solve");
    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    temp.data = solver.solve(builder.b);
    t.stop();

    t.start("solve inc");
    solver.analyzePattern(builderInc.A);
    solver.factorize(builderInc.A);
    dTemp.data = solver.solve(builderInc.b);
    std::cout << "dTemp: " << dTemp.data.transpose() << std::endl;
    t.stop();

    t.start("output");
    tempInc.data += dTemp.data;
    io.print({temp, tempInc}, time);
    t.stop();
  }

  // std::cout << "temp:\n" << temp.data << std::endl;

  Vec exact;
  interpolateAnalyticFunction(exactSol, feSpace, exact);
  Vec error;
  error = temp.data - exact;

  t.print();

  double norm = error.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm
            << std::endl;
  return checkError({norm}, {2.87785419773588e-07});
  return 0;
}
