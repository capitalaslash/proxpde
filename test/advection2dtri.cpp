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

#include <iostream>

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
using FVSolver_T = FVSolver<FESpaceP0_T, LimiterType::SUPERBEE>;
// flux feSpace
using FacetMesh_T = Mesh<Elem_T::Facet_T>;
using FacetFESpace_T = FESpace<FacetMesh_T,
                               FEType<Elem_T::Facet_T,0>::RefFE_T,
                               FEType<Elem_T::Facet_T,0>::RecommendedQR>;
using VelFESpace_T = FESpace<Mesh_T,
                             FEType<Elem_T,1>::RefFE_T,
                             FEType<Elem_T,1>::RecommendedQR, 2>;


template <typename Elem>
void buildFacetMesh(Mesh<typename Elem::Facet_T> & facetMesh, Mesh<Elem> const & mesh)
{
  facetMesh.pointList = mesh.pointList;
  facetMesh.elementList = mesh.facetList;
  facetMesh.buildConnectivity();
}

static scalarFun_T ic = [] (Vec3 const& p)
{
  // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50-(p(1)-0.7)*(p(1)-0.7)*50);
  if (p(0) < .37) return 1.;
  return 0.;
};

static Fun<2, 3> velFun = [] (Vec3 const & /*p*/)
{
  // double const r = std::sqrt((p(0)-0.5)*(p(0)-0.5) + (p(1)-0.5)*(p(1)-0.5));
  // return Vec2(.25*(p(1)-0.5), -.25*(p(0)-0.5));
  return Vec2{0.2, 0.};
};

int main(int argc, char* argv[])
{
  MilliTimer t;

  t.start("mesh");
  // we need internal facets
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  // hexagonSquare(*mesh, true);
  uint const numElemsX = (argc < 3)? 10 : std::stoul(argv[1]);
  uint const numElemsY = (argc < 3)? 10 : std::stoul(argv[2]);
  buildHyperCube(*mesh,
            {0., 0., 0.},
            {1., 1., 0.},
            {{numElemsX, numElemsY, 0}},
            INTERNAL_FACETS | NORMALS);
  // readGMSH(*mesh, "square_uns.msh", INTERNAL_FACETS | NORMALS);

  std::unique_ptr<FacetMesh_T> facetMesh{new FacetMesh_T};
  buildFacetMesh(*facetMesh, *mesh);
  t.stop();

  t.start("fespace");
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
  FacetFESpace_T facetFESpace{*facetMesh};
  t.stop();

  t.start("bcs");
  // auto const zero = [](Vec3 const & ){return 0.;};
  auto const one = [](Vec3 const & ){return 1.;};
  BCList bcsP1{feSpaceP1};
  bcsP1.addBC(BCEss{feSpaceP1, side::LEFT, one});
  // bcsP1.addBC(BCEss{feSpaceP1, side::LEFT, zero});
  // bcsP1.addBC(BCEss{feSpaceP1, side::BOTTOM, zero});
  // bcsP1.addBC(BCEss{feSpaceP1, side::RIGHT, zero});
  // bcsP1.addBC(BCEss{feSpaceP1, side::TOP, zero});
  BCList bcsP0{feSpaceP0};
  bcsP0.addBC(BCEss{feSpaceP0, side::LEFT, one});
  // bcsP0.addBC(BCEss{feSpaceP0, side::LEFT, zero});
  // bcsP0.addBC(BCEss{feSpaceP0, side::BOTTOM, zero});
  // bcsP0.addBC(BCEss{feSpaceP0, side::RIGHT, zero});
  // bcsP0.addBC(BCEss{feSpaceP0, side::TOP, zero});
  t.stop();

  t.start("velocity");
  auto const & sizeP1 = feSpaceP1.dof.size;
  VelFESpace_T velFESpace{*mesh};
  Var velFE("velocity", velFESpace.dof.size*2);
  for (uint d=0; d<2; d++)
  {
    for (uint i=0; i< sizeP1; i++)
    {
      velFE.data(velFESpace.dof.ptMap[i] + d*sizeP1) =
          velFun(mesh->pointList[i].coord)[d];
    }
  }
  double const dt = 0.1;
  auto const cfl = computeMaxCFL(velFESpace, velFE.data, dt);
  std::cout << "max cfl = " << cfl << std::endl;
  t.stop();

  t.start("p1 assembly setup");
  Builder builder{sizeP1};
  LUSolver solver;
  // Vec velOldStyle = Vec::Zero(2*sizeP1);
  // velOldStyle.block(     0, 0, sizeP1, 1) = Vec::Constant(sizeP1, velocity[0]);
  // velOldStyle.block(sizeP1, 0, sizeP1, 1) = Vec::Constant(sizeP1, velocity[1]);
  AssemblyAdvection advection(1.0, velFE.data, feSpaceP1);
  AssemblyMass timeder(1./dt, feSpaceP1);
  Vec concP1Old(sizeP1);
  AssemblyProjection timeder_rhs(1./dt, concP1Old, feSpaceP1);
  t.stop();

  t.start("init");
  Var concP1{"concP1"};
  interpolateAnalyticFunction(ic, feSpaceP1, concP1.data);
  Var concP0{"concP0"};
  // we need to use the highest order available QR to integrate discontinuous functions
  // FESpace<Mesh_T, RefTriangleP0, GaussQR<Triangle, 7>> feSpaceIC{*mesh};
  FESpace<Mesh_T, RefTriangleP0, MiniQR<Triangle, 10>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, concP0.data);
  // interpolateAnalyticFunction(ic, feSpaceP0, concP0.data);
  Var flux{"flux"};
  flux.data = Vec::Zero(static_cast<uint>(facetMesh->elementList.size()), 1);
  t.stop();

  FVSolver_T fv{feSpaceP0, bcsP0};
  Table<double, 2> vel(sizeP1, 2);
  //vel.block(0, 0, sizeP1, 1) = Vec::Constant(sizeP1, 1, velocity[0]);
  //vel.block(0, 1, sizeP1, 1) = Vec::Constant(sizeP1, 1, velocity[1]);
  vel.block(0, 0, sizeP1, 1) = velFE.data.block(0, 0, sizeP1, 1);
  vel.block(0, 1, sizeP1, 1) = velFE.data.block(sizeP1, 0, sizeP1, 1);

  uint const ntime = 50;
  double time = 0.0;
  IOManager ioP1{feSpaceP1, "output_advection2dtri/solP1"};
  ioP1.print({concP1});
  IOManager ioP0{feSpaceP0, "output_advection2dtri/solP0"};
  ioP0.print({concP0});
  IOManager ioFlux{facetFESpace, "output_advection2dtri/flux"};
  ioFlux.print({flux});
  IOManager ioVel{velFESpace, "output_advection2dtri/vel"};
  ioVel.print({velFE});

  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << std::endl;

    // central implicit
    t.start("p1 assemby");
    concP1Old = concP1.data;
    builder.buildProblem(timeder, bcsP1);
    builder.buildProblem(timeder_rhs, bcsP1);
    builder.buildProblem(advection, bcsP1);
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
    ioP1.time = time;
    ioP1.iter += 1;
    ioP1.print({concP1});
    ioP0.time = time;
    ioP0.iter += 1;
    ioP0.print({concP0});
    flux.data = fv.fluxes;
    ioFlux.time = time;
    ioFlux.iter += 1;
    ioFlux.print({flux});
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
  if (std::fabs(errorNormP1 - 0.2572474581492306) > 1.e-10 ||
      std::fabs(errorNormP0 - 0.03789136711107713) > 1.e-10)
  // if (std::fabs(errorNormP1 - 0.5313772283037004) > 1.e-10 ||
  //     std::fabs(errorNormP0 - 0.01568656361986407) > 1.e-10)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
