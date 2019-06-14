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

int main(/*int argc, char* argv[]*/)
{
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpaceP0_T = FESpace<Mesh_T,
                              FEType<Elem_T, 0>::RefFE_T,
                              FEType<Elem_T, 0>::RecommendedQR>;
  using FESpaceP1_T = FESpace<Mesh_T,
                              FEType<Elem_T, 1>::RefFE_T,
                              FEType<Elem_T, 1>::RecommendedQR>;
  using FESpaceVel_T = FESpace<Mesh_T,
                               FEType<Elem_T, 1>::RefFE_T,
                               FEType<Elem_T, 1>::RecommendedQR, 2>;
  using FVSolver_T = FVSolver<FESpaceP0_T, LimiterType::SUPERBEE>;

  scalarFun_T const ic = [] (Vec3 const& p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
    if(p(0) < .37) return 1.;
    return 0.;
  };

  Fun<2, 3> const velFun = [] (Vec3 const & /*p*/)
  {
    return Vec2{0.1, 0.};
    // return Vec2(.25*(p(1)-0.5), -.25*(p(0)-0.5));
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
  //       INTERNAL_FACETS | NORMALS);
  readGMSH(*mesh, "square_q.msh", INTERNAL_FACETS | NORMALS);
  t.stop();

  t.start("fespace");
  FESpaceP0_T feSpace{*mesh};
  FESpaceP1_T feSpaceP1{*mesh};
  t.stop();

  t.start("bcs");
  auto const one = [](Vec3 const & ){return 1.;};
  BCList bcs{feSpace};
  bcs.addBC(BCEss{feSpace, side::LEFT, one});
  t.stop();

  t.start("velocity");
  FESpaceVel_T feSpaceVel{*mesh};
  Var velFE{"velocity"};
  interpolateAnalyticFunction(velFun, feSpaceVel, velFE.data);
  double const dt = 0.1;
  auto const cfl = computeMaxCFL(feSpaceVel, velFE.data, dt);
  std::cout << "max cfl = " << cfl << std::endl;
  t.stop();

  t.start("init");
  Var c{"conc"};
  Vec cOld(feSpace.dof.size);
  FESpace<Mesh_T, FEType<Elem_T, 0>::RefFE_T, MiniQR<Elem_T, 10>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, c.data);
  t.stop();

  FVSolver_T fv{feSpace, bcs};

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

  uint const ntime = 200;
  double time = 0.0;
  IOManager io{feSpace, "output_advection2dquad/sol"};
  io.print({c});
  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << std::endl;

    t.start("update");
    cOld = c.data;
    fv.update(c.data);
    fv.computeFluxes(vel, feSpaceP1);
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
  std::cout << "the norm of the error is " << std::setprecision(16) << errorNorm << std::endl;
  if (std::fabs(errorNorm - 6.389801046171856e-05) > 1.e-10)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
