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
  using FESpaceVel_T = FESpace<Mesh_T,
                               FEType<Elem_T,1>::RefFE_T,
                               FEType<Elem_T,1>::RecommendedQR, Elem_T::dim>;
  // flux feSpace
  using MeshFacet_T = Mesh<Elem_T::Facet_T>;
  using FESpaceFacet_T = FESpace<MeshFacet_T,
                                 FEType<Elem_T::Facet_T,0>::RefFE_T,
                                 FEType<Elem_T::Facet_T,0>::RecommendedQR>;

  scalarFun_T const ic = [] (Vec3 const& p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50-(p(1)-0.7)*(p(1)-0.7)*50);
    if (p(0) < .37) return 1.;
    return 0.;
  };

  Fun<2, 3> const velFun = [] (Vec3 const & /*p*/)
  {
    // double const r = std::sqrt((p(0)-0.5)*(p(0)-0.5) + (p(1)-0.5)*(p(1)-0.5));
    // return Vec2(.25*(p(1)-0.5), -.25*(p(0)-0.5));
    return Vec2{0.2, 0.};
  };

  MilliTimer t;

  t.start("mesh");
  // we need internal facets
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // hexagonSquare(*mesh, true);
  // buildNormals(*mesh);

  uint const numElemsX = (argc < 3)? 10 : std::stoul(argv[1]);
  uint const numElemsY = (argc < 3)? 10 : std::stoul(argv[2]);
  buildHyperCube(*mesh,
            {0., 0., 0.},
            {1., 1., 0.},
            {{numElemsX, numElemsY, 0}},
            INTERNAL_FACETS | NORMALS);

  // readGMSH(*mesh, "square_uns.msh", INTERNAL_FACETS | NORMALS);

  std::unique_ptr<MeshFacet_T> facetMesh{new MeshFacet_T};
  buildFacetMesh(*facetMesh, *mesh);
  t.stop();

  t.start("fespace");
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
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
  FESpaceVel_T feSpaceVel{*mesh};
  Var velFE{"velocity"};
  interpolateAnalyticFunction(velFun, feSpaceVel, velFE.data);
  double const dt = 0.1;
  auto const cfl = computeMaxCFL(feSpaceVel, velFE.data, dt);
  std::cout << "max cfl = " << cfl << std::endl;
  t.stop();

  t.start("p1 assembly setup");
  auto const & sizeP1 = feSpaceP1.dof.size;
  Builder builder{sizeP1};
  LUSolver solver;
  AssemblyMass timeDer(1./dt, feSpaceP1);
  AssemblyAdvection advection(1.0, velFE.data, feSpaceVel, feSpaceP1);
  Vec concP1Old{sizeP1};
  AssemblyProjection timeDerRhs(1./dt, concP1Old, feSpaceP1);
  t.stop();

  t.start("init");
  Var concP1{"concP1"};
  interpolateAnalyticFunction(ic, feSpaceP1, concP1.data);
  Var concP0{"concP0"};
  // we need to use the highest order available QR to integrate discontinuous functions
  // FESpace<Mesh_T, FEType<Elem_T, 0>::RefFE_T, GaussQR<Triangle, 7>> feSpaceIC{*mesh};
  FESpace<Mesh_T, FEType<Elem_T, 0>::RefFE_T, MiniQR<Elem_T, 10>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, concP0.data);
  // interpolateAnalyticFunction(ic, feSpaceP0, concP0.data);
  Var flux{"flux"};
  flux.data = Vec::Zero(static_cast<uint>(facetMesh->elementList.size()), 1);
  t.stop();

  FVSolver fv{feSpaceP0, bcsP0, SuperBEELimiter{}};

  auto const & sizeVel = feSpaceVel.dof.size;
  Table<double, FESpaceVel_T::dim> vel(sizeVel, FESpaceVel_T::dim);
  if (FESpaceVel_T::DOF_T::ordering == DofOrdering::BLOCK)
  {
    for (uint k=0; k<FESpaceVel_T::dim; ++k)
    {
      vel.block(0, k, sizeVel, 1) = velFE.data.block(k*sizeVel, 0, sizeVel, 1);
    }
  }
  else // FESpaceVel_T::DOF_T::ordering == DofOrdering::INTERLEAVED
  {
    for (uint r=0; r<sizeVel; ++r)
    {
      vel.row(r) = velFE.data.block(r*FESpaceVel_T::dim, 0, FESpaceVel_T::dim, 1).transpose();
    }
  }

  uint const ntime = 50;
  double time = 0.0;
  IOManager ioP1{feSpaceP1, "output_advection2dtri/solP1"};
  ioP1.print({concP1});
  IOManager ioP0{feSpaceP0, "output_advection2dtri/solP0"};
  ioP0.print({concP0});
  IOManager ioFlux{feSpaceFacet, "output_advection2dtri/flux"};
  ioFlux.print({flux});
  IOManager ioVel{feSpaceVel, "output_advection2dtri/vel"};
  ioVel.print({velFE});

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
    fv.computeFluxes(vel, feSpaceP1);
    fv.advance(concP0.data, dt);
    t.stop();

    // print
    t.start("print");
    ioP1.print({concP1}, time);
    ioP0.print({concP0}, time);
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
  return checkError({errorNormP1, errorNormP0}, {0.2572474581492306, 0.03789136711107713});
  // if (std::fabs(errorNormP1 - 0.5313772283037004) > 1.e-10 ||
  //     std::fabs(errorNormP0 - 0.01568656361986407) > 1.e-10 ||
  //     std::isnan(errorNormP1) ||
  //     std::isnan(errorNormP0))
}
