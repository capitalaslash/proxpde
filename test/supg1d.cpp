#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "fv.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

namespace proxpde
{

// int constexpr sigma = +1; // MS
int constexpr sigma = 0; // SUPG
// int constexpr sigma = -1; // GLS

int constexpr orderCG = 1;

template <typename FESpace, typename Vel>
struct AssemblyLhsSUPG: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;

  explicit AssemblyLhsSUPG(
      double const c,
      Vel & u,
      FESpace_T & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      idt(c),
      vel(u)
  {
    // this works only if the same quad rule is defined on all fe spaces
    static_assert(
        std::is_same_v<typename FESpace_T::QR_T, typename Vel::FESpace_T::QR_T>,
        "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace->curFE.reinit(elem);
    vel.reinit(elem);
  }

  void build(LMat_T & Ke) const override
  {
    using RefFE_T = typename FESpace_T::RefFE_T;
    using CurFE_T = typename FESpace_T::CurFE_T;

    uint constexpr size = CurFE_T::numDOFs;
    CurFE_T const & curFE = this->feSpace->curFE;

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      auto const velQPoint = promote<3>(vel.evaluate(q));
      double const iVelNorm = 1. / velQPoint.norm();

      FVec<RefFE_T::numDOFs> phiSUPG = curFE.phi[q] + 0.5 * curFE.elem->hMin() *
                                                          iVelNorm *
                                                          (curFE.dphi[q] * velQPoint);

      for (uint d = 0; d < FESpace_T::dim; ++d)
      {
        Ke.template block<size, size>(d * size, d * size) +=
            curFE.JxW[q] * phiSUPG *
            (idt * curFE.phi[q] + curFE.dphi[q] * velQPoint).transpose();
      }
    }
  }

  double const idt;
  Vel & vel;
};

template <typename FESpace, typename Rhs, typename Vel>
struct AssemblyRhsSUPG: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyRhsSUPG(
      double const c,
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
        std::is_same_v<typename FESpace_T::QR_T, typename Rhs::FESpace_T::QR_T>,
        "the two quad rule are not the same");
    static_assert(
        std::is_same_v<typename FESpace_T::QR_T, typename Vel::FESpace_T::QR_T>,
        "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace->curFE.reinit(elem);
    uOld.reinit(elem);
    vel.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using RefFE_T = typename FESpace_T::RefFE_T;
    using CurFE_T = typename FESpace_T::CurFE_T;

    uint constexpr size = CurFE_T::numDOFs;
    CurFE_T const & curFE = this->feSpace->curFE;

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      auto const uOldQPoint = uOld.evaluate(q);
      auto const velQPoint = promote<3>(vel.evaluate(q));
      double const iVelNorm = 1. / velQPoint.norm();

      FVec<RefFE_T::numDOFs> phiSUPG = curFE.phi[q] + 0.5 * curFE.elem->hMin() *
                                                          iVelNorm *
                                                          (curFE.dphi[q] * velQPoint);

      for (uint d = 0; d < FESpace_T::dim; ++d)
      {
        Fe.template block<size, 1>(d * size, 0) +=
            curFE.JxW[q] * phiSUPG * (idt * uOldQPoint[d]);
      }
    }
  }

  double idt;
  Rhs & uOld;
  Vel & vel;
};

template <typename FESpace, typename Vel>
struct AssemblyLhsSUPGRes: public Diagonal<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = Diagonal<FESpace>;
  using LMat_T = typename Super_T::LMat_T;

  explicit AssemblyLhsSUPGRes(
      double const c,
      Vel & u,
      double const e,
      FESpace_T & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace>()):
      Diagonal<FESpace>(fe, cmp),
      idt(c),
      vel(u),
      eps(e)
  {
    // this works only if the same quad rule is defined on all fe spaces
    static_assert(
        std::is_same_v<typename FESpace_T::QR_T, typename Vel::FESpace_T::QR_T>,
        "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace->curFE.reinit(elem);
    vel.reinit(elem);
  }

  void build(LMat_T & Ke) const override
  {
    using RefFE_T = typename FESpace_T::RefFE_T;
    using CurFE_T = typename FESpace_T::CurFE_T;

    uint constexpr size = CurFE_T::numDOFs;
    CurFE_T const & curFE = this->feSpace->curFE;

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      auto const velQPoint = promote<3>(vel.evaluate(q));
      double const iVelNorm = 1. / velQPoint.norm();

      FVec<RefFE_T::numDOFs> phiSUPG =
          .5 * curFE.elem->hMin() * iVelNorm * (curFE.dphi[q] * velQPoint)
#ifdef PROXPDE_ENABLE_SECONDDERIV
          + sigma * eps *
                (curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(0, 0).col(0) +
                 curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(1, 0).col(1) +
                 curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(2, 0).col(2))
#endif
          ;

      auto const res =
          idt * curFE.phi[q] + curFE.dphi[q] * velQPoint
#ifdef PROXPDE_ENABLE_SECONDDERIV
          + eps * (curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(0, 0).col(0) +
                   curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(1, 0).col(1) +
                   curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(2, 0).col(2))
#endif
          ;
      for (uint d = 0; d < FESpace_T::dim; ++d)
      {
        Ke.template block<size, size>(d * size, d * size) +=
            curFE.JxW[q] * phiSUPG * res.transpose();
      }
    }
  }

  double const idt;
  Vel & vel;
  double const eps;
};

template <typename FESpace, typename Rhs, typename Vel>
struct AssemblyRhsSUPGRes: public AssemblyVector<FESpace>
{
  using FESpace_T = FESpace;
  using Super_T = AssemblyVector<FESpace>;
  using LVec_T = typename Super_T::LVec_T;

  AssemblyRhsSUPGRes(
      double const c,
      Rhs & u,
      Vel & v,
      double const e,
      FESpace_T & fe,
      AssemblyBase::CompList const & cmp = allComp<FESpace_T>()):
      AssemblyVector<FESpace_T>(fe, cmp),
      idt(c),
      uOld(u),
      vel(v),
      eps(e)
  {
    // this works only if the same quad rule is defined on all fe spaces
    static_assert(
        std::is_same_v<typename FESpace_T::QR_T, typename Rhs::FESpace_T::QR_T>,
        "the two quad rule are not the same");
    static_assert(
        std::is_same_v<typename FESpace_T::QR_T, typename Vel::FESpace_T::QR_T>,
        "the two quad rule are not the same");
  }

  void reinit(GeoElem const & elem) const override
  {
    this->feSpace->curFE.reinit(elem);
    uOld.reinit(elem);
    vel.reinit(elem);
  }

  void build(LVec_T & Fe) const override
  {
    using RefFE_T = typename FESpace_T::RefFE_T;
    using CurFE_T = typename FESpace_T::CurFE_T;

    uint constexpr size = CurFE_T::numDOFs;
    CurFE_T const & curFE = this->feSpace->curFE;

    for (uint q = 0; q < CurFE_T::QR_T::numPts; ++q)
    {
      auto const uOldQPoint = uOld.evaluate(q);
      auto const velQPoint = promote<3>(vel.evaluate(q));
      double const iVelNorm = 1. / velQPoint.norm();

      FVec<RefFE_T::numDOFs> phiSUPG =
          .5 * curFE.elem->hMin() * iVelNorm * (curFE.dphi[q] * velQPoint)
#ifdef PROXPDE_ENABLE_SECONDDERIV
          + sigma * eps *
                (curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(0, 0).col(0) +
                 curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(1, 0).col(1) +
                 curFE.d2phi[q].template block<RefFE_T::numDOFs, 3>(2, 0).col(2))
#endif
          ;
      for (uint d = 0; d < FESpace_T::dim; ++d)
      {
        Fe.template block<size, 1>(d * size, 0) +=
            curFE.JxW[q] * phiSUPG * (idt * uOldQPoint[d]);
      }
    }
  }

  double idt;
  Rhs & uOld;
  Vel & vel;
  double eps;
};

} // namespace proxpde

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Line;
  using Mesh_T = Mesh<Elem_T>;
  // implicit finite element central
  using FESpaceCG_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, orderCG>::RefFE_T,
      LagrangeFE<Elem_T, orderCG>::RecommendedQR>;
  // explicit finite volume upwind
  using FESpaceFV_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 0>::RefFE_T,
      LagrangeFE<Elem_T, 0>::RecommendedQR>;
  // velocity field
  using FESpaceVel_T = FESpaceCG_T;

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
    config["final_time"] = 3.;
    config["velocity"] = 0.5;
    config["threshold"] = 0.37;
  }

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 0., 0.};
  uint const numElems = config["n"].as<uint>();
  buildHyperCube(
      *mesh,
      origin,
      length,
      {numElems, 0, 0},
      MeshFlags::INTERNAL_FACETS | MeshFlags::NORMALS);
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
  FESpaceCG_T feSpaceCG{*mesh};
  FESpaceFV_T feSpaceFV{*mesh};
  t.stop();

  t.start("bcs");
  auto const one = [](Vec3 const &) { return 1.; };
  // auto const zero = [] (Vec3 const & ) { return 0.; };
  auto const bcsCG = std::vector{BCEss{feSpaceCG, side::LEFT, one}};
  auto const bcsFV = std::vector{BCEss{feSpaceFV, side::LEFT, one}};
  t.stop();

  auto const dt = config["dt"].as<double>();
  // the only 1D velocity that is divergence free is a constant one
  auto const velocity = config["velocity"].as<double>();
  FESpaceVel_T feSpaceVel{*mesh};
  FEVar vel{feSpaceVel};
  vel << velocity;
  double const cfl = velocity * dt * numElems;
  std::cout << "cfl = " << cfl << std::endl;
  assert(cfl < 1. - 1.e-8);

  Builder builder{feSpaceCG.dof.size};
  LUSolver solver;
  AssemblyLhsSUPG lhsSUPG{1. / dt, vel, feSpaceCG};
  // AssemblyStiffness artificialDiff(0.5 * velocity / numElems, feSpaceCG);
  FEVar uCGOld{feSpaceCG};
  AssemblyRhsSUPG rhsSUPG{1. / dt, uCGOld, vel, feSpaceCG};

  Builder builderRes{feSpaceCG.dof.size};
  LUSolver solverRes;
  AssemblyMass timeDer{1. / dt, feSpaceCG};
  AssemblyAdvection advection{1.0, vel.data, *vel.feSpace, feSpaceCG};
  double const eps = 0.001;
  AssemblyStiffness diffusion{eps, feSpaceCG};
  AssemblyLhsSUPGRes lhsSUPGRes{1. / dt, vel, eps, feSpaceCG};
  FEVar uCGOldRes{feSpaceCG};
  AssemblyProjection timeDerRhs{1. / dt, uCGOldRes.data, *uCGOldRes.feSpace, feSpaceCG};
  AssemblyRhsSUPGRes rhsSUPGRes{1. / dt, uCGOldRes, vel, eps, feSpaceCG};

  // FEVar c{"conc", feSpace};
  // FEList feList{feSpace, feSpace};
  // auto fe1 = std::get<0>(feList);
  // BlockFEVar tmp{"tmp", feList};
  auto const threshold = config["threshold"].as<double>();
  scalarFun_T ic = [threshold](Vec3 const & p)
  {
    // return std::exp(-(p(0)-0.5)*(p(0)-0.5)*50);
    if (p(0) < threshold)
      return 1.;
    return 0.;
  };
  // scalarFun_T ic = [] (Vec3 const& p)
  // {
  //   if (p(0) > 0.25 && p(0) < 0.5) return 1.;
  //   return 0.;
  // };
  FEVar uCG{"uCG", feSpaceCG};
  uCG << ic;
  // uCG.data[8] = 0.8;
  // integrateAnalyticFunction(ic, feSpaceCG, uCG.data);
  FEVar uCGRes{"uCGres", feSpaceCG};
  uCGRes << ic;
  // uCGRes.data[8] = 0.8;
  // integrateAnalyticFunction(ic, feSpaceCG, uCGRes.data);

  t.start("p0 ic");
  Var uFV{"uFV"};
  // we need to use the highest order available QR to integrate discontinuous functions
  // FESpace<Mesh_T, RefLineP0, GaussQR<Line, 4>> feSpaceIC{*mesh};
  FESpace<Mesh_T, RefLineP0, MiniQR<Line, 20>> feSpaceIC{*mesh};
  integrateAnalyticFunction(ic, feSpaceIC, uFV.data);
  t.stop();

  FVSolver fv{feSpaceFV, bcsFV, MinModLimiter{}};

  IOManager ioCG{feSpaceCG, "output_supg1d/cg"};
  ioCG.print(std::tuple{uCG, uCGRes});
  IOManager ioFV{feSpaceFV, "output_supg1d/fv"};
  ioFV.print({uFV});

  auto const lhs = std::tuple{lhsSUPG};
  auto const rhs = std::tuple{rhsSUPG};

  auto const lhsRes = std::tuple{timeDer, advection, diffusion, lhsSUPGRes};
  auto const rhsRes = std::tuple{timeDerRhs, rhsSUPGRes};

  double maxFV = 0.;
  double maxCG = 0.;
  double maxCGRes = 0.;

  auto const ntime =
      static_cast<uint>(std::nearbyint(config["final_time"].as<double>() / dt));
  double time = 0.0;
  for (uint itime = 0; itime < ntime; itime++)
  {
    time += dt;
    std::cout << "solving timestep " << itime << ", time = " << time << std::endl;

    // central implicit
    t.start("p1 assemby");
    uCGOld.data = uCG.data;
    builder.buildLhs(lhs, bcsCG);
    builder.buildRhs(rhs, bcsCG);
    builder.closeMatrix();
    t.stop();

    uCGOldRes.data = uCGRes.data;
    builderRes.buildLhs(lhsRes, bcsCG);
    builderRes.buildRhs(rhsRes, bcsCG);
    builderRes.closeMatrix();

    t.start("p1 solve");
    solver.compute(builder.A);
    uCG.data = solver.solve(builder.b);
    builder.clear();
    t.stop();
    maxCG = std::max(maxCG, uCG.data.maxCoeff());

    // std::cout << "A:\n" << builder.A << std::endl;
    // std::cout << "b:\n" << builder.b << std::endl;
    // std::cout << "sol:\n" << c.data << std::endl;

    solverRes.compute(builderRes.A);
    uCGRes.data = solverRes.solve(builderRes.b);
    builderRes.clear();
    maxCGRes = std::max(maxCGRes, uCGRes.data.maxCoeff());

    // explicit upwind
    t.start("p0 update");
    fv.update(uFV.data);
    fv.computeFluxes(vel);
    fv.advance(uFV.data, dt);
    t.stop();
    maxFV = std::max(maxFV, uFV.data.maxCoeff());

    // print
    t.start("print");
    ioCG.print(std::tuple{uCG, uCGRes}, time);
    ioFV.print({uFV}, time);
    t.stop();
  }

  t.print();

  std::cout << "maxFV:    " << maxFV << std::endl;
  std::cout << "maxCG:    " << maxCG << std::endl;
  std::cout << "maxCGRes: " << maxCGRes << std::endl;

  Vec oneFieldCG;
  interpolateAnalyticFunction(one, feSpaceCG, oneFieldCG);
  Vec oneFieldFV;
  interpolateAnalyticFunction(one, feSpaceFV, oneFieldFV);

  double errorNormCG = (uCG.data - oneFieldCG).norm();
  std::cout << "the norm of the CG error is " << std::setprecision(16) << errorNormCG
            << std::endl;
  double errorNormFV = (uFV.data - oneFieldFV).norm();
  std::cout << "the norm of the FV error is " << std::setprecision(16) << errorNormFV
            << std::endl;
  return checkError(
      {errorNormCG, errorNormFV}, {0.0007185655616791794, 2.68467866957268e-06});
}
