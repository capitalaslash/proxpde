#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"
#include "fv.hpp"
#include "feutils.hpp"

int main(int argc, char* argv[])
{
  using Elem_T = Triangle;
  using Mesh_T = Mesh<Elem_T>;
  // implicit finite element central
  using FESpaceP1_T = FESpace<Mesh_T,
                              FEType<Elem_T,1>::RefFE_T,
                              FEType<Elem_T,1>::RecommendedQR>;
  // explicit finite volume upwind
  using FESpaceP0_T = FESpace<Mesh_T,
                              FEType<Elem_T,0>::RefFE_T,
                              FEType<Elem_T,0>::RecommendedQR>;
  // velocity field
  // using FESpaceVel_T = FESpace<Mesh_T,
  //                              FEType<Elem_T,1>::RefFE_T,
  //                              FEType<Elem_T,1>::RecommendedQR, Elem_T::dim>;
  using FESpaceVel_T = FESpace<Mesh_T,
                               RefTriangleRT0,
                               GaussQR<Triangle, 3>>;
  using FESpaceVelP0_T = FESpace<Mesh_T,
                                 RefTriangleP0,
                                 GaussQR<Triangle, 3>, 2>;
  // flux feSpace
  using MeshFacet_T = Mesh<Elem_T::Facet_T>;
  using FESpaceFacet_T = FESpace<MeshFacet_T,
                                 FEType<Elem_T::Facet_T,0>::RefFE_T,
                                 FEType<Elem_T::Facet_T,0>::RecommendedQR>;

  ParameterDict config;
  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["nx"] = 10U;
    config["ny"] = 10U;

    config["velocity"] = Vec3{0.2, 0., 0.};
    config["threshold"] = 0.37;

    config["dt"] = 0.1;
    config["ntime"] = 50U;
  }
  config.validate({"nx", "ny", "dt", "velocity", "threshold", "ntime"});

  scalarFun_T const ic = [threshold = config["threshold"].as<double>()] (Vec3 const& p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50-(p(1)-0.7)*(p(1)-0.7)*50);
    if (p(0) < .37) return 1.;
    return 0.;
  };

  MilliTimer t;

  t.start("mesh");
  // we need internal facets
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(
        *mesh,
        {0., 0., 0.},
        {1., 1., 0.},
        {config["nx"].as<uint>(), config["ny"].as<uint>(), 0},
        INTERNAL_FACETS | NORMALS | FACET_PTRS);
  // readGMSH(*mesh, "square_uns.msh", INTERNAL_FACETS | NORMALS);
  // hexagonSquare(*mesh, true);
  // buildNormals(*mesh);


  std::unique_ptr<MeshFacet_T> facetMesh{new MeshFacet_T};
  buildFacetMesh(*facetMesh, *mesh);
  t.stop();

  t.start("fespace");
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceVelP0_T feSpaceVelP0{*mesh};
  FESpaceFacet_T feSpaceFacet{*facetMesh};
  t.stop();

  t.start("bcs");
  // auto const zero = [](Vec3 const & ){return 0.;};
  auto const one = [](Vec3 const & ){return 1.;};
  auto bcLeftP1 = BCEss{feSpaceP1, side::LEFT};
  bcLeftP1 << one;
  // bcLeftP1 << zero;
  // auto bcBottomP1 = BCEss{feSpaceP1, side::BOTTOM};
  // bcbottomP1 << zero;
  // auto bcRightP1 = BCEss{feSpaceP1, side::RIGHT};
  // bcRightP1 << zero;
  // auto bcTopP1 = BCEss{feSpaceP1, side::TOP};
  // bcTopP1 << zero;
  auto const bcsP1 = std::make_tuple(bcLeftP1/*, bcBottomP1, bcRightP1, bcTopP1*/);
  auto bcLeftP0 = BCEss{feSpaceP0, side::LEFT};
  bcLeftP0 << one;
  // bcLeftP0 << zero;
  // auto bcBottomP0 = BCEss{feSpaceP0, side::BOTTOM};
  // bcBottomP0 << zero;
  // auto bcRightP0 = BCEss{feSpaceP0, side::RIGHT};
  // bcRightP0 << zero;
  // auto bcTopP0 = BCEss{feSpaceP0, side::TOP};
  // bcTopP0 << zero;
  auto const bcsP0 = std::make_tuple(bcLeftP0/*, bcBottomP0, bcRightP0, bcTopP0*/);
  t.stop();

  t.start("velocity");
  FEVar vel{feSpaceVel, "velocity"};
  vel << narrow<Elem_T::dim>(config["velocity"].as<Vec3>());

  double const dt = config["dt"].as<double>();
  auto const cfl = computeMaxCFL(feSpaceVel, vel.data, dt);
  std::cout << "max cfl = " << cfl << std::endl;
  t.stop();

  t.start("p1 assembly setup");
  auto const & sizeP1 = feSpaceP1.dof.size;
  Builder builder{sizeP1};
  LUSolver solver;
  AssemblyScalarMass timeDer(1./dt, feSpaceP1);
  AssemblyAdvectionFE advection(1.0, vel, feSpaceP1);
  Vec concP1Old{sizeP1};
  AssemblyProjection timeDerRhs(1./dt, concP1Old, feSpaceP1);
  t.stop();

  t.start("init");
  FEVar concP1{feSpaceP1, "concP1"};
  concP1 << ic;
  FEVar concP0{feSpaceP0, "concP0"};
  // we need to use the highest order available QR to integrate discontinuous functions
  // FESpace<Mesh_T, FEType<Elem_T, 0>::RefFE_T, GaussQR<Triangle, 7>> feSpaceIC{*mesh};
  FESpace<Mesh_T, FEType<Elem_T, 0>::RefFE_T, MiniQR<Elem_T, 10>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, concP0.data);
  Var flux{"flux"};
  flux.data = Vec::Zero(static_cast<uint>(facetMesh->elementList.size()), 1);
  t.stop();

  FVSolver fv{feSpaceP0, bcsP0, SuperBEELimiter{}};

  uint const ntime = config["ntime"].as<uint>();
  double time = 0.0;
  IOManager ioP1{feSpaceP1, "output_advection2dtri/solP1"};
  ioP1.print(std::array{concP1});
  IOManager ioP0{feSpaceP0, "output_advection2dtri/solP0"};
  ioP0.print(std::array{concP0});
  IOManager ioFlux{feSpaceFacet, "output_advection2dtri/flux"};
  ioFlux.print({flux});
  Var velP0{"velP0", feSpaceVelP0.dof.size * FESpaceVelP0_T::dim};
  l2Projection(velP0.data, feSpaceVelP0, vel.data, feSpaceVel);
  IOManager ioVel{feSpaceVelP0, "output_advection2dtri/vel"};
  ioVel.print(std::array{velP0});

  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << std::endl;

    // central implicit
    t.start("p1 assemby");
    concP1Old = concP1.data;
    builder.buildLhs(std::tuple{timeDer, advection}, bcsP1);
    builder.buildRhs(std::tuple{timeDerRhs}, bcsP1);
    builder.closeMatrix();
    t.stop();

    t.start("p1 solve");
    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    concP1.data = solver.solve(builder.b);
    builder.clear();
    t.stop();

    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;
    // std::cout << "sol:\n" << c.data << std::endl;

    // explicit upwind
    t.start("p0 update");
    fv.update(concP0.data);
    fv.computeFluxes(vel);
    fv.advance(concP0.data, dt);
    t.stop();

    // print
    t.start("print");
    ioP1.print(std::array{concP1}, time);
    ioP0.print(std::array{concP0}, time);
    flux.data = fv.fluxes;
    ioFlux.print({flux}, time);
    t.stop();
  }

  t.print();

  Vec oneFieldP1;
  interpolateAnalyticFunction(one, feSpaceP1, oneFieldP1);
  Vec oneFieldP0;
  interpolateAnalyticFunction(one, feSpaceP0, oneFieldP0);

  double errorNormP1 = (concP1.data - oneFieldP1).norm();
  std::cout << "the norm of the P1 error is " << std::setprecision(16) << errorNormP1 << std::endl;
  double errorNormP0 = (concP0.data - oneFieldP0).norm();
  std::cout << "the norm of the P0 error is " << std::setprecision(16) << errorNormP0 << std::endl;
  return checkError(
    {errorNormP1, errorNormP0},
    {0.2572474581492306, 0.03789136711107713});
}
