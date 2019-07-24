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

template <typename FESpace, typename Vel>
struct AssemblyLhsSUPG: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;

  explicit AssemblyLhsSUPG(double const c,
                           Vel & u,
                           FESpace_T & fe,
                           AssemblyBase::CompList const & cmp = allComp<FESpace>()):
    Diagonal<FESpace>(fe, cmp),
    idt(c),
    vel(u)
  {
    // this works only if the same quad rule is defined on all fe spaces
    static_assert(
          std::is_same_v<
          typename FESpace_T::QR_T,
          typename Vel::FESpace_T::QR_T>,
          "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    vel.reinit(elem);
  }

  void build(LMat_T & Ke) const override
  {
    using RefFE_T = typename FESpace_T::RefFE_T;
    using CurFE_T = typename FESpace_T::CurFE_T;

    uint constexpr size = CurFE_T::numDOFs;
    CurFE_T const & curFE = this->feSpace.curFE;

    for(uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      auto const velQPoint = promote<3>(vel.evaluate(q));

      FVec<RefFE_T::numFuns> phiSUPG =
          curFE.phi[q] +
          0.5 * curFE.elem->hMin() * (curFE.dphi[q] * velQPoint);

      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Ke.template block<size, size>(d * size, d * size) +=
            curFE.JxW[q] *
            phiSUPG *
            (idt * curFE.phi[q] +
             curFE.dphi[q] * velQPoint).transpose();
      }
    }
  }

  double const idt;
  Vel & vel;
};

template <typename FESpace,
          typename Rhs,
          typename Vel>
struct AssemblyRhsSUPG: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyRhsSUPG(double const c,
                  Rhs & u,
                  Vel & v,
                  FESpace_T & fe,
                  AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
    AssemblyVector<FESpace_T>(fe, cmp),
    idt(c),
    uOld(u),
    vel(v)
  {
    // this works only if the same quad rule is defined on all fe spaces
    static_assert(
          std::is_same_v<
          typename FESpace_T::QR_T,
          typename Rhs::FESpace_T::QR_T>,
          "the two quad rule are not the same");
    static_assert(
          std::is_same_v<
          typename FESpace_T::QR_T,
          typename Vel::FESpace_T::QR_T>,
          "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace.curFE.reinit(elem);
    uOld.reinit(elem);
    vel.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using RefFE_T = typename FESpace_T::RefFE_T;
    using CurFE_T = typename FESpace_T::CurFE_T;

    uint constexpr size = CurFE_T::numDOFs;
    CurFE_T const & curFE = this->feSpace.curFE;

    for (uint q=0; q<CurFE_T::QR_T::numPts; ++q)
    {
      auto const uOldQPoint = uOld.evaluate(q);
      auto const velQPoint = promote<3>(vel.evaluate(q));

      FVec<RefFE_T::numFuns> phiSUPG =
          curFE.phi[q] +
          0.5 * curFE.elem->hMin() * (curFE.dphi[q] * velQPoint);

      for (uint d=0; d<FESpace_T::dim; ++d)
      {
        Fe.template block<size, 1>(d* size, 0) +=
            curFE.JxW[q] *
            phiSUPG *
            (idt * uOldQPoint[d]);
      }
    }
  }

  double idt;
  Rhs & uOld;
  Vel & vel;
};

int main(int argc, char* argv[])
{
  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  // implicit finite element central
  using FESpaceP1_T = FESpace<Mesh_T,
                              FEType<Elem_T, 1>::RefFE_T,
                              FEType<Elem_T, 1>::RecommendedQR>;
  // explicit finite volume upwind
  using FESpaceP0_T = FESpace<Mesh_T,
                              FEType<Elem_T, 0>::RefFE_T,
                              FEType<Elem_T, 0>::RecommendedQR>;
  // velocity field
  using FESpaceVel_T = FESpaceP1_T;

  MilliTimer t;

  YAML::Node config;
  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["n"] = 10U;
    config["dt"] = 0.1;
    config["final_time"] = 3.0;
    config["velocity"] = 0.5;
    config["threshold"] = 0.37;
  }

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  uint const numElems = config["n"].as<uint>();
  buildHyperCube(*mesh, origin, length, {numElems, 0, 0}, INTERNAL_FACETS | NORMALS);
  // auto const r = 1. / 1.1;
  // auto const starting = (1. - r) / (1. - pow(r, numElems));
  // auto counter = 0;
  // for (auto & p: mesh->pointList)
  // {
  //   p.coord[0] = starting * (1. - pow(r, counter)) / (1. - r);
  //   counter++;
  // }
  t.stop();

  t.start("fespace");
  FESpaceP1_T feSpaceP1{*mesh};
  FESpaceP0_T feSpaceP0{*mesh};
  t.stop();

  t.start("bcs");
  auto const one = [](Vec3 const & ){return 1.;};
  auto bcLeftP1 = BCEss{feSpaceP1, side::LEFT};
  bcLeftP1 << one;
  auto const bcsP1 = std::make_tuple(bcLeftP1);
  auto bcLeftP0 = BCEss{feSpaceP0, side::LEFT};
  bcLeftP0 << one;
  auto const bcsP0 = std::make_tuple(bcLeftP0);
  t.stop();

  auto const dt = config["dt"].as<double>();
  // the only 1D velocity that is divergence free is a constant one
  auto const velocity = config["velocity"].as<double>();
  FESpaceVel_T feSpaceVel{*mesh};
  FEVar vel{feSpaceVel};
  vel << velocity;
  double const hinv = numElems;
  std::cout << "cfl = " << velocity * dt * hinv << std::endl;

  Builder builder{feSpaceP1.dof.size};
  LUSolver solver;
  AssemblyLhsSUPG advection(1./dt, vel, feSpaceP1);
  // AssemblyStiffness artificialDiff(0.5 * velocity / numElems, feSpaceP1);
  FEVar concP1Old{feSpaceP1};
  AssemblyRhsSUPG timeDerRhs(1./dt, concP1Old, vel, feSpaceP1);

  // FEVar c{feSpace, "conc"};
  // FEList feList{feSpace, feSpace};
  // auto fe1 = std::get<0>(feList);
  // BlockFEVar tmp{"tmp", feList};
  auto const threshold = config["threshold"].as<double>();
  scalarFun_T ic = [threshold] (Vec3 const& p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
    if (p(0) < threshold) return 1.;
    return 0.;
  };
  Var concP1{"conc"};
  interpolateAnalyticFunction(ic, feSpaceP1, concP1.data);

  t.start("p0 ic");
  Var concP0{"concP0"};
  // we need to use the highest order available QR to integrate discontinuous functions
  // FESpace<Mesh_T, RefLineP0, GaussQR<Line, 4>> feSpaceIC{*mesh};
  FESpace<Mesh_T, RefLineP0, MiniQR<Line, 20>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, concP0.data);
  t.stop();

  FVSolver fv{feSpaceP0, bcsP0, MinModLimiter{}};

  IOManager ioP1{feSpaceP1, "output_supg1d/solP1"};
  ioP1.print({concP1});
  IOManager ioP0{feSpaceP0, "output_supg1d/solP0"};
  ioP0.print({concP0});

  auto const lhs = std::tuple{advection};
  auto const rhs = std::tuple{timeDerRhs};

  auto const ntime = static_cast<uint>(std::nearbyint(config["final_time"].as<double>() / dt));
  double time = 0.0;
  for(uint itime=0; itime<ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << ", time = " << time << std::endl;

    // central implicit
    t.start("p1 assemby");
    concP1Old.data = concP1.data;
    builder.buildLhs(lhs, bcsP1);
    builder.buildRhs(rhs, bcsP1);
    builder.closeMatrix();
    t.stop();

    t.start("p1 solve");
    solver.compute(builder.A);
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
    ioP1.print({concP1}, time);
    ioP0.print({concP0}, time);
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
  return checkError({errorNormP1, errorNormP0}, {0.0005419506968368649, 2.68467866957268e-06});
}