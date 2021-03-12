#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

template <typename FESpace>
struct BCPeriodic
{
  using FESpace_T = FESpace;

  BCPeriodic(FESpace_T const & feSpace,
             marker_T const mDest,
             marker_T const mOrig,
             Fun<3,3> const destToOrigFun):
    bcDest{feSpace, mDest},
    bcOrig{feSpace, mOrig}
  {
    // this works only for translated meshes
    assert(bcOrig._constrainedDofMap.size() == bcDest._constrainedDofMap.size());

    std::map<id_T, id_T> invPtMap;
    id_T counter = 0;
    for (auto const id: feSpace.dof.ptMap)
    {
      invPtMap[id] = counter;
      counter++;
    }

    for (auto const & idDest: bcDest._constrainedDofMap)
    {
      Vec3 const pDest = feSpace.mesh.pointList[invPtMap[idDest.first]].coord;
      // std::cout << "pDest: " << pDest.transpose() << std::endl;

      Vec3 const pConv = destToOrigFun(pDest);
      // std::cout << "pConv: " << pConv.transpose() << std::endl;

      for (auto const & idOrig: bcOrig._constrainedDofMap)
      {
        Vec3 const pOrig = feSpace.mesh.pointList[invPtMap[idOrig.first]].coord;
        if ((pOrig - pConv).norm() < 1.e-15)
        {
          destToOrig[idDest.first] = idOrig.first;
          break;
        }
      }
    }

    for (auto const & [k, v]: destToOrig)
    {
      std::cout << k << " -> " << v << std::endl;
    }
  }

  BCEss<FESpace_T> bcDest;
  BCEss<FESpace_T> bcOrig;
  std::unordered_map<id_T, id_T> destToOrig;
};

int main(int argc, char* argv[])
{
  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  using FESpaceP1_T = FESpace<Mesh_T,
                              LagrangeFE<Elem_T, 1>::RefFE_T,
                              LagrangeFE<Elem_T, 1>::RecommendedQR>;

  MilliTimer t;

  YAML::Node config;
  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["n"] = 20u;
    config["dt"] = 0.05;
    config["nu"] = 1.e-4;
    // the only 1D velocity that is divergence free is a constant one
    config["velocity"] = 0.4;
    config["ntime"] = 100u;
  }

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  uint const numElems = config["n"].as<uint>();
  buildHyperCube(*mesh, origin, length, {numElems, 0, 0}, INTERNAL_FACETS | NORMALS);
  t.stop();

  t.start("fespace");
  FESpaceP1_T feSpace{*mesh};
  t.stop();

  t.start("bc");
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << [] (Vec3 const &) { return 0.; };
  auto bcRight = BCEss{feSpace, side::RIGHT};
  bcRight << [] (Vec3 const &) { return 0.; };
  auto const bcPeriodic = BCPeriodic{feSpace, side::LEFT, side::RIGHT, [](Vec3 const &){return Vec3{1.0, 0., 0.};}};
  auto const bcs = std::make_tuple(/*bcLeft, bcRight*/);
  t.stop();

  auto const dt = config["dt"].as<double>();
  auto const nu = config["nu"].as<double>();
  auto const velocity = config["velocity"].as<double>();
  FEVar vel{feSpace};
  vel << velocity;
  double const hinv = numElems;
  std::cout << "cfl = " << velocity * dt * hinv << std::endl;

  double const a = 20.;
  scalarFun_T ic = [a] (Vec3 const& p)
  {
    return std::exp(-(p(0)-0.5)*(p(0)-0.5) * a);
  };
  // \int_{-\infty}^{+\infty} exp(-a (x - b)^2 dx = \sqrt(\pi / a)
  double const refIntegral = std::sqrt(M_PI / a);

  FEVar u{feSpace, "u"};
  u << ic;

  t.start("assembly lhs");
  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{
                     AssemblyScalarMass{1./dt, feSpace},
                     AssemblyStiffness{nu, feSpace},
                     AssemblyAdvection{1.0, vel.data, feSpace, feSpace}},
                   bcs);
  auto const penalty = 1.e10;
  for (auto const & [dest, orig]: bcPeriodic.destToOrig)
  {
    // diagonal term
    builder._triplets.emplace_back(dest, dest, penalty);
    // away term
    builder._triplets.emplace_back(dest, orig, -penalty);
  }
  builder.closeMatrix();
  Vec uOld{feSpace.dof.size};
  auto assemblyRhs = AssemblyProjection{1./dt, uOld, feSpace};
  t.stop();

  auto const ntime = config["ntime"].as<uint>();
  double time = 0.0;
  IOManager io{feSpace, "output_periodic/sol"};
  io.print(std::tuple{u});
  std::cout << "solution integral: " << u.integrate() << std::endl;

  t.start("solver");
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  t.stop();

  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << ", time = " << time << std::endl;

    t.start("assembly rhs");
    uOld = u.data;
    builder.buildRhs(std::tuple{assemblyRhs}, bcs);
    t.stop();

    t.start("solver");
    u.data = solver.solve(builder.b);
    t.stop();

    // print
    io.print(std::tuple{u}, time);
    std::cout << "solution integral: " << u.integrate() << std::endl;

    builder.clearRhs();
  }

  t.print();

  double const error = std::fabs(u.integrate()[0] - refIntegral);
  std::cout << "the error in the solution is " << std::setprecision(16) << error << std::endl;
  return checkError({error}, {0.0005238803728871422});
}
