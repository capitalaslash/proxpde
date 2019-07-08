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
using FVSolver_T = FVSolver<FESpaceP0_T, LimiterType::UPWIND>;
using FESpaceVel_T = FESpace<Mesh_T,
                             FEType<Elem_T,1>::RefFE_T,
                             FEType<Elem_T,1>::RecommendedQR, 2>;
static scalarFun_T ic = [] (Vec3 const& p)
{
  // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
  if (p(0) < .4) return 1.;
  return 0.;
};

double areaTriangle(Vec3 p0, Vec3 p1, Vec3 p2)
{
  return 0.5 * ((p1 - p0).cross(p2 - p0)).norm();
}

struct DualCell //: public GeoElem
{
  using elemNeighborList_T = std::vector<Triangle const *>;
  using edgeNeighborList_T = std::vector<Line const *>;

  void build()
  {
    assert(edgeNeighbors.size() > 0);
    volumes.resize(edgeNeighbors.size());
    for (uint n=0; n<edgeNeighbors.size(); ++n)
    {
      auto const & edge = *edgeNeighbors[n];
      // points.row(n) = edge.midpoint();
      volumes[n] = areaTriangle(center->coord, edge.midpoint(), edge.facingElem[0].ptr->midpoint());
      if (edge.facingElem[1].ptr != nullptr)
      {
        volumes[n] += areaTriangle(center->coord, edge.midpoint(), edge.facingElem[1].ptr->midpoint());
      }
    }
  }

  Point const * center;
  elemNeighborList_T elemNeighbors;
  edgeNeighborList_T edgeNeighbors;
  std::vector<double> volumes;
};

static Fun<2, 3> velFun = [] (Vec3 const &)
{
  return Vec2(0.2, 0.0);
};

int main(int argc, char* argv[])
{
  MilliTimer t;

  t.start("mesh build");
  // we need internal facets
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  hexagonSquare(*mesh, true);
  addElemFacetList(*mesh);

  // uint const numElemsX = (argc < 3)? 16 : std::stoul(argv[1]);
  // uint const numElemsY = (argc < 3)? 16 : std::stoul(argv[2]);
  // buildHyperCube(
  //       *mesh,
  //       {0., 0., 0.},
  //       {1., 1., 0.},
  //       {{numElemsX, numElemsY, 0}},
  //       INTERNAL_FACETS | NORMALS);

  // readGMSH(*mesh, "square_uns.msh");
  // buildFacets(*mesh, true);
  // buildNormals(*mesh);
  t.stop();

  t.start("fe spaces");
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
  t.stop();

  t.start("dual");
  std::vector<DualCell> dualCellList(mesh->pointList.size());
  for (auto const & p: mesh->pointList)
  {
    dualCellList[p.id].center = &p;
  }
  for (auto const & e: mesh->elementList)
  {
    for (auto const & p: e.pointList)
    {
      dualCellList[p->id].elemNeighbors.push_back(&e);
    }
  }
  // TODO: should be using edges (both in 2d and 3d)
  for (auto const & f: mesh->facetList)
  {
    dualCellList[f.pointList[0]->id].edgeNeighbors.push_back(&f);
    dualCellList[f.pointList[1]->id].edgeNeighbors.push_back(&f);
  }

  Var dualVolumes("dual_volumes", mesh->pointList.size());
  for (auto & c: dualCellList)
  {
    c.build();
    dualVolumes.data[feSpaceP1.dof.ptMap[c.center->id]] =
        std::accumulate(c.volumes.begin(), c.volumes.end(), 0.);
  }
  t.stop();

  t.start("bcs");
  auto const leftBC = [](Vec3 const &){return 1.;};
  BCList bcsP1{feSpaceP1};
  bcsP1.addBC(BCEss{feSpaceP1, side::LEFT, leftBC});
  BCList bcsP0{feSpaceP0};
  bcsP0.addBC(BCEss{feSpaceP0, side::LEFT, leftBC});
  t.stop();

  t.start("fe builder");
  FESpaceVel_T feSpaceVel{*mesh};
  Var velFE{"velocity"};
  interpolateAnalyticFunction(velFun, feSpaceVel, velFE.data);
  double const dt = 0.1;
  auto const cfl = computeMaxCFL(feSpaceVel, velFE.data, dt);
  std::cout << "max cfl = " << cfl << std::endl;

  auto const & sizeP1 = feSpaceP1.dof.size;
  Builder builder{sizeP1};
  LUSolver solver;
  AssemblyAdvection advection(1.0, velFE.data, feSpaceVel, feSpaceP1);
  AssemblyMass timeder(1./dt, feSpaceP1);
  Vec concP1Old(sizeP1);
  AssemblyProjection timeder_rhs(1./dt, concP1Old, feSpaceP1);

  Var concP1{"concP1"};
  interpolateAnalyticFunction(ic, feSpaceP1, concP1.data);
  Var concP0{"concP0"};
  interpolateAnalyticFunction(ic, feSpaceP0, concP0.data);

  FVSolver_T fv{feSpaceP0, bcsP0};
  Table<double, 2> vel(sizeP1, 2);
  vel.block(0, 0, sizeP1, 1) = velFE.data.block(0, 0, sizeP1, 1);
  vel.block(0, 1, sizeP1, 1) = velFE.data.block(sizeP1, 0, sizeP1, 1);
  t.stop();

  uint const ntime = 50;
  double time = 0.0;
  IOManager ioP1{feSpaceP1, "output_dual/solP1"};
  ioP1.print({concP1, dualVolumes});
  IOManager ioP0{feSpaceP0, "output_dual/solP0"};
  ioP0.print({concP0});

  t.start("time iter");
  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << std::endl;

    // central implicit
    concP1Old = concP1.data;

    builder.buildLhs(timeder, bcsP1);
    builder.buildRhs(timeder_rhs, bcsP1);
    builder.buildLhs(advection, bcsP1);
    builder.closeMatrix();

    solver.analyzePattern(builder.A);
    solver.factorize(builder.A);
    concP1.data = solver.solve(builder.b);
    builder.clear();

    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;
    // std::cout << "sol:\n" << c.data << std::endl;

    // explicit upwind
    fv.uOld = concP0.data;
    fv.computeFluxes(vel, feSpaceP1);
    fv.advance(concP0.data, dt);

    // print
    ioP1.print({concP1}, time);
    ioP0.print({concP0}, time);
  }
  t.stop();

  t.print();

  double normP1 = concP1.data.norm();
  std::cout << "the norm of the P1 solution is " << std::setprecision(16) << normP1 << std::endl;
  double normP0 = concP0.data.norm();
  std::cout << "the norm of the P0 solution is " << std::setprecision(16) << normP0 << std::endl;
  // if(std::fabs(normP1 - 10.8594759676) > 1.e-10 || std::fabs(normP0 - 13.7780000857) > 1.e-10)
  // {
  //   std::cerr << "the norm of the solution is not the prescribed value" << std::endl;
  //   return 1;
  // }

  return 0;
}
