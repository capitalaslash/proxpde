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
// flux feSpace
using FacetMesh_T = Mesh<Elem_T::Facet_T>;
using FacetFESpace_T = FESpace<FacetMesh_T,
                               FEType<Elem_T::Facet_T,0>::RefFE_T,
                               FEType<Elem_T::Facet_T,0>::RecommendedQR>;

template <typename Elem>
void buildFacetMesh(Mesh<typename Elem::Facet_T> & facetMesh, Mesh<Elem> const & mesh)
{
  facetMesh.pointList = mesh.pointList;
  facetMesh.elementList = mesh.facetList;
  facetMesh.buildConnectivity();
}

static scalarFun_T ic = [] (Vec3 const& p)
{
  // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
  if (p(0) < .4) return 1.;
  return 0.;
};

int main(int argc, char* argv[])
{
  MilliTimer t;

  t.start();
  // we need internal facets
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  hexagonSquare(*mesh, true);
  // MeshBuilder<Elem_T> meshBuilder;
  // uint const numPtsX = (argc < 3)? 17 : std::stoul(argv[1]);
  // uint const numPtsY = (argc < 3)? 17 : std::stoul(argv[2]);
  // meshBuilder.build(
  //       *mesh,
  //       {0., 0., 0.},
  //       {1., 1., 0.},
  //       {{numPtsX, numPtsY, 0}},
  //       true);
  // readGMSH(*mesh, "square_uns.msh");
  // buildFacets(*mesh, true);
  buildNormals(*mesh);

  std::unique_ptr<FacetMesh_T> facetMesh{new FacetMesh_T};
  buildFacetMesh(*facetMesh, *mesh);
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
  FacetFESpace_T facetFESpace{*facetMesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  auto const leftBC = [](Vec3 const &){return 1.;};
  BCList bcsP1{feSpaceP1};
  bcsP1.addEssentialBC(side::LEFT, leftBC);
  BCList bcsP0{feSpaceP0};
  bcsP0.addEssentialBC(side::LEFT, leftBC);
  std::cout << "bcs: " << t << " ms" << std::endl;

  auto const velocity = Vec2(0.2, 0.0);
  double const dt = 0.1;
  auto const cfl = computeMaxCFL(*mesh, velocity, dt);
  std::cout << "max cfl = " << cfl << std::endl;

  auto const & sizeP1 = feSpaceP1.dof.size;
  Builder builder{sizeP1};
  LUSolver solver;
  Vec velOldStyle = Vec::Zero(2*sizeP1);
  velOldStyle.block(   0, 0, sizeP1, 1) = Vec::Constant(sizeP1, velocity[0]);
  velOldStyle.block(sizeP1, 0, sizeP1, 1) = Vec::Constant(sizeP1, velocity[1]);
  AssemblyAdvection advection(1.0, velOldStyle, feSpaceP1);
  AssemblyMass timeder(1./dt, feSpaceP1);
  Vec concP1Old(sizeP1);
  AssemblyProjection timeder_rhs(1./dt, concP1Old, feSpaceP1);

  Var concP1{"concP1"};
  interpolateAnalyticFunction(ic, feSpaceP1, concP1.data);
  Var concP0{"concP0"};
  interpolateAnalyticFunction(ic, feSpaceP0, concP0.data);
  Var flux{"flux"};
  flux.data = Vec::Zero(static_cast<uint>(facetMesh->elementList.size()), 1);

  FVSolver fv{feSpaceP0, bcsP0};
  Table<double, 2> vel(sizeP1, 2);
  vel.block(0, 0, sizeP1, 1) = Vec::Constant(sizeP1, 1, velocity[0]);
  vel.block(0, 1, sizeP1, 1) = Vec::Constant(sizeP1, 1, velocity[1]);

  uint const ntime = 50;
  double time = 0.0;
  IOManager ioP1{feSpaceP1, "output_advection2dtri/solP1"};
  ioP1.print({concP1});
  IOManager ioP0{feSpaceP0, "output_advection2dtri/solP0"};
  ioP0.print({concP0});
  IOManager ioFlux{facetFESpace, "output_advection2dtri/flux"};
  ioFlux.print({flux});

  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << std::endl;

    // central implicit
    concP1Old = concP1.data;

    builder.buildProblem(timeder, bcsP1);
    builder.buildProblem(timeder_rhs, bcsP1);
    builder.buildProblem(advection, bcsP1);
    builder.closeMatrix();

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    concP1.data = solver.solve(builder.b);
    builder.clear();

    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;
    // std::cout << "sol:\n" << c.data << std::endl;

    // explicit upwind
    fv.update(concP0.data);
    fv.computeFluxes(vel, feSpaceP1);
    fv.advance(concP0.data, dt);

    // print
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
  }

  double normP1 = concP1.data.norm();
  std::cout << "the norm of the P1 solution is " << std::setprecision(12) << normP1 << std::endl;
  double normP0 = concP0.data.norm();
  std::cout << "the norm of the P0 solution is " << std::setprecision(12) << normP0 << std::endl;
  if(std::fabs(normP1 - 10.8594759676) > 1.e-10 || std::fabs(normP0 - 13.7780000857) > 1.e-10)
  {
    std::cerr << "the norm of the solution is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
