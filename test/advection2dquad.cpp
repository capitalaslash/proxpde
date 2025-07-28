#include "def.hpp"

#include <yaml-cpp/yaml.h>

#include "bc.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "fv.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

int main(/*int argc, char * argv[]*/)
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpaceP0_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 0>::RefFE_T,
      LagrangeFE<Elem_T, 0>::RecommendedQR>;
  // using FESpaceP1_T = FESpace<
  //     Mesh_T,
  //     LagrangeFE<Elem_T, 1>::RefFE_T,
  //     LagrangeFE<Elem_T, 1>::RecommendedQR>;
  using FESpaceVel_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR,
      2>;

  scalarFun_T const ic = [](Vec3 const & p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
    if (p(0) < .37)
      return 1.;
    return 0.;
  };

  MilliTimer t;

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // uint const numElemsX = (argc < 3)? 10 : std::stoi(argv[1]);
  // uint const numElemsY = (argc < 3)? 10 : std::stoi(argv[2]);
  // buildHyperCube(
  //       *mesh,
  //       {0., 0., 0.},
  //       {1., 1., 0.},
  //       {numElemsX, numElemsY, 0},
  //       MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS);
  readGMSH(*mesh, "square_q.msh", MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS);
  t.stop();

  t.start("fespace");
  FESpaceP0_T feSpace{*mesh};
  // FESpaceP1_T feSpaceP1{*mesh};
  t.stop();

  t.start("bcs");
  auto const one = [](Vec3 const &) { return 1.; };
  auto bc = BCEss{feSpace, side::LEFT};
  bc << one;
  auto const bcs = std::vector{bc};
  t.stop();

  t.start("velocity");
  FESpaceVel_T feSpaceVel{*mesh};
  FEVar vel{"velocity", feSpaceVel};
  vel << Vec2{0.1, 0.};

  double const dt = 0.1;
  auto const cfl = computeMaxCFL(feSpaceVel, vel.data, dt);
  std::cout << "max cfl = " << cfl << std::endl;
  t.stop();

  t.start("init");
  FEVar c{"conc", feSpace};
  Vec cOld(feSpace.dof.size);
  FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, MiniQR<Elem_T, 10>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, c.data);
  t.stop();

  FVSolver fv{feSpace, bcs, SuperBEELimiter{}};

  uint const ntime = 200;
  double time = 0.0;
  IOManager io{feSpace, "output_advection2dquad/sol"};
  io.print({c});
  for (uint itime = 0; itime < ntime; itime++)
  {
    time += dt;
    fmt::println("solving timestep {}", itime);

    t.start("update");
    cOld = c.data;
    fv.update(c.data);
    fv.computeFluxes(vel);
    fv.advance(c.data, dt);
    t.stop();

    // print
    t.start("print");
    io.print({c}, time);
    t.stop();
  }

  t.print();

  Vec oneField;
  interpolateAnalyticFunction(one, feSpace, oneField);

  double errorNorm = (c.data - oneField).norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << errorNorm
            << std::endl;
  return checkError({errorNorm}, {6.389801046171856e-05});
}
