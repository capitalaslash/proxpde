#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "eqn.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "feutils.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "timer.hpp"

template <Component Comp, typename FESpaceOut, typename FESpaceIn>
void interpolateOnBoundary(
    Vec & out,
    FESpaceOut const & feSpaceOut,
    Vec const & in,
    FESpaceIn const & feSpaceIn,
    std::unordered_set<marker_T> const & markers)
{
  uint constexpr dim = FESpaceIn::dim;
  using FEFacet_T = typename FESpaceIn::RefFE_T::FacetFE_T;
  using QRFacet_T = SideQR_T<typename FESpaceIn::QR_T>;
  using CurFEFacet_T = CurFE<FEFacet_T, QRFacet_T>;
  using FESpaceFacet_T = FESpace<
      typename FESpaceIn::Mesh_T,
      typename FESpaceIn::RefFE_T,
      SideGaussQR<typename FESpaceIn::Mesh_T::Elem_T, QRFacet_T::numPts>,
      dim>;

  out = Vec::Zero(feSpaceOut.dof.size);
  CurFEFacet_T curFEFacet;

  FESpaceFacet_T feSpaceFacet{feSpaceIn.mesh};
  FEVar inFacet{feSpaceFacet};
  inFacet.data = in;

  for (auto const & facet: feSpaceOut.mesh.facetList)
  {
    if (markers.contains(facet.marker))
    {
      curFEFacet.reinit(facet);

      auto const & elem = *(facet.facingElem[0].ptr);
      auto const side = facet.facingElem[0].side;
      auto const id = feSpaceOut.dof.getId(elem.id);
      inFacet.reinit(elem);

      // TODO: add SCALAR with Vec1 dotted
      Vec3 dotted;
      if constexpr (Comp == Component::NORMAL)
      {
        dotted = facet.normal();
      }
      else if constexpr (Comp == Component::TANGENTIAL)
      {
        // TODO: make it a class method
        dotted = (facet.pointList[1]->coord - facet.pointList[0]->coord).normalized();
      }

      for (uint q = 0; q < QRFacet_T::numPts; ++q)
      {
        Vec3 const inLocal = inFacet.evaluate(side * QRFacet_T::numPts + q);
        out[id] += curFEFacet.JxW[q] * inLocal.dot(dotted);
      }
    }
  }
}

int main(int argc, char * argv[])
{
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;

  using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
  // using LinearQR = LagrangeFE<Elem_T,1>::RecommendedQR;

  // monolithic
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, Elem_T::dim>;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;

  // split
  using FESpaceU_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR>;

  // split RT
  using FESpaceUStarRT_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;
  using FESpaceVelStarRT_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR, Elem_T::dim>;
  using FESpaceRT0_T =
      FESpace<Mesh_T, RaviartThomasFE<Elem_T, 0>::RefFE_T, QuadraticQR>;
  using FESpaceLambda_T = FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, QuadraticQR>;

  MilliTimer t;

  ParameterDict config;

  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["mesh"]["origin"] = Vec3{0.0, 0.0, 0.0};
    config["mesh"]["length"] = Vec3{1.0, 2.0, 0.0};
    config["mesh"]["n"] = std::array{5U, 5U, 0U};
    config["mesh"]["flags"] =
        MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS;
    config["dt"] = 0.2;
    config["ntime"] = 50U;
    config["nu"] = 0.01;
    config["printStep"] = 1U;
    config["monolithic"] = true;
    config["split"] = true;
    config["split_rt"] = true;
  }
  config.validate(
      {"mesh", "dt", "ntime", "nu", "printStep", "monolithic", "split", "split_rt"});

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildHyperCube(*mesh, ParameterDict{config["mesh"]});
  t.stop();

  t.start("fespace");
  // monolithic
  FESpaceVel_T feSpaceVel{*mesh};
  FESpaceP_T feSpaceP{*mesh, feSpaceVel.dof.size * FESpaceVel_T::dim};

  // split
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpacePSplit{*mesh};

  // split RT
  FESpaceUStarRT_T feSpaceUStarRT{*mesh};
  FESpaceVelStarRT_T feSpaceVelStarRT{*mesh};
  FESpaceRT0_T feSpaceVelRT{*mesh};
  FESpaceLambda_T feSpaceLambda{*mesh};
  FESpaceLambda_T feSpaceLambdaCoupling{*mesh, feSpaceVelRT.dof.size};
  t.stop();

  t.start("bc");
  // monolithic
  // auto const inlet = [] (Vec3 const &) { return Vec2(0.0, 1.0); };
  auto const zeroV = [](Vec3 const &) { return Vec2{0., 0.}; };
  auto const zeroS = [](Vec3 const &) { return 0.; };
  // auto const inlet = [](Vec3 const & p) { return Vec2{0.0, 1.5 * (1. - p[0] * p[0])};
  // };
  auto const inlet = [](Vec3 const &) { return Vec2{0.0, 0.0}; };
  auto const inlet0 = [&inlet](Vec3 const & p) { return inlet(p)[0]; };
  auto const inlet1 = [&inlet](Vec3 const & p) { return inlet(p)[1]; };

  // last essential bc wins on corners
  auto bcsVel = std::tuple{
      BCEss{feSpaceVel, side::BOTTOM},
      BCEss{feSpaceVel, side::RIGHT},
      BCEss{feSpaceVel, side::LEFT, {0}}};
  std::get<0>(bcsVel) << inlet;
  std::get<1>(bcsVel) << zeroV;
  std::get<2>(bcsVel) << zeroV;
  // TODO: pin a single DOF
  auto bcTopP = BCEss{feSpaceP, side::TOP};
  bcTopP << zeroS;
  auto const bcsP = std::tuple{bcTopP};
  // DofSet_T pinSet = {1};
  // bcsP.addBC(BCEss{feSpaceP, pinSet, zeroS});

  // split
  auto bcsUStar = std::tuple{
      BCEss{feSpaceU, side::BOTTOM},
      BCEss{feSpaceU, side::RIGHT},
      BCEss{feSpaceU, side::LEFT}};
  std::get<0>(bcsUStar) << inlet0;
  std::get<1>(bcsUStar) << zeroS;
  std::get<2>(bcsUStar) << zeroS;
  auto bcsVStar =
      std::tuple{BCEss{feSpaceU, side::BOTTOM}, BCEss{feSpaceU, side::RIGHT}};
  std::get<0>(bcsVStar) << inlet1;
  std::get<1>(bcsVStar) << zeroS;
  auto bcTopPSplit = BCEss{feSpacePSplit, side::TOP};
  bcTopPSplit << zeroS;
  auto const bcsPSplit = std::tuple{bcTopPSplit};
  // eqnP.bcList.addBC(BCEss{eqnP.feSpace, side::BOTTOM, pIn});

  // split RT
  auto bcsUStarRT = std::tuple{
      BCEss{feSpaceUStarRT, side::BOTTOM},
      BCEss{feSpaceUStarRT, side::RIGHT},
      BCEss{feSpaceUStarRT, side::LEFT}};
  std::get<0>(bcsUStarRT) << inlet0;
  std::get<1>(bcsUStarRT) << zeroS;
  std::get<2>(bcsUStarRT) << zeroS;
  auto bcsVStarRT = std::tuple{
      BCEss{feSpaceUStarRT, side::BOTTOM}, BCEss{feSpaceUStarRT, side::RIGHT}};
  std::get<0>(bcsVStarRT) << inlet1;
  std::get<1>(bcsVStarRT) << zeroS;
  auto bcsVelRT = std::tuple{
      BCEss{feSpaceVelRT, side::BOTTOM},
      BCEss{feSpaceVelRT, side::RIGHT},
      BCEss{feSpaceVelRT, side::LEFT}};
  // TODO: do integration on facet
  // std::get<0>(bcsVelRT) << [] (Vec3 const & p) {
  //   if (p[0] < 0.5)
  //     return -11./8.;
  //   return -5./8.;
  // };
  std::get<0>(bcsVelRT) << [&inlet](Vec3 const & p)
  {
    return promote<3>(inlet(p));
    // if (p[0] < 0.75)
    //   return -4./ 3.;
    // return 0.;
    // return -1.;
  };
  auto const zero3d = [](Vec3 const &) { return Vec3{0.0, 0.0, 0.0}; };
  std::get<1>(bcsVelRT) << zero3d;
  std::get<2>(bcsVelRT) << zero3d;
  // auto const bcsVelRT = std::tuple{};
  auto const bcsLambda = std::tuple{};
  t.stop();

  auto const dt = config["dt"].as<double>();
  auto const nu = config["nu"].as<double>();
  auto const g = [](Vec3 const &) { return Vec2{0.0, -1.0}; };
  auto const g0 = [&g](Vec3 const & p) { return g(p)[0]; };
  auto const g1 = [&g](Vec3 const & p) { return g(p)[1]; };

  t.start("assembly monolithic");
  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;

  // vel / dt + (velOld \cdot \nabla) vel - \nabla \cdot (nu \nabla vel) + \nabla p =
  // velOld / dt
  // - \nabla \cdot vel = 0
  Vec velOldMonolithic{dofU * FESpaceVel_T::dim};
  AssemblyScalarMass timeder(1. / dt, feSpaceVel);
  AssemblyAdvection advection(1.0, velOldMonolithic, feSpaceVel, feSpaceVel);
  AssemblyTensorStiffness stiffness(nu, feSpaceVel);
  AssemblyGrad grad(-1.0, feSpaceVel, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceVel);
  AssemblyProjection timederRhs(1. / dt, velOldMonolithic, feSpaceVel);
  AssemblyAnalyticRhs gravity(g, feSpaceVel);

  Builder<StorageType::RowMajor> builderM{dofU * FESpaceVel_T::dim + dofP};
  // builderM.buildLhs(dummy, bcsP);
  // builderM.buildLhs(std::tuple{timeder, stiffness}, bcsVel);
  // builderM.buildCoupling(grad, bcsVel, bcsP);
  // builderM.buildCoupling(div, bcsP, bcsVel);
  // builderM.closeMatrix();
  // Mat<StorageType::RowMajor> matFixed = builderM.A;
  // Vec rhsFixed = builderM.b;
  // builderM.clear();
  IterSolver solverM;
  t.stop();

  t.start("assembly split");
  Var uStar{"uStar", feSpaceU.dof.size};
  Var vStar{"vStar", feSpaceU.dof.size};
  Var p{"p", feSpacePSplit.dof.size};
  Var u{"u", feSpaceU.dof.size};
  Var v{"v", feSpaceU.dof.size};

  Vec velStar{dofU * FESpaceVel_T::dim};
  Vec pOld = Vec::Zero(dofP);
  Vec vel{dofU * FESpaceVel_T::dim};

  // uStar / dt + (vel \cdot \nabla) uStar - \nabla \cdot (nu \nabla uStar) =
  // u / dt - d pOld / dx
  auto const uStarLhs = std::tuple{
      AssemblyScalarMass{1. / dt, feSpaceU},
      AssemblyAdvection{1.0, vel, feSpaceVel, feSpaceU},
      AssemblyStiffness{nu, feSpaceU},
  };
  auto const uStarRhs = std::tuple{
      AssemblyProjection{1. / dt, u.data, feSpaceU},
      AssemblyGradRhs2{1.0, pOld, feSpacePSplit, feSpaceU, {0}},
      AssemblyAnalyticRhs{g0, feSpaceU},
  };

  // vStar / dt + (vel \cdot \nabla) vStar - \nabla \cdot (nu \nabla vStar) =
  // v / dt - d pOld / dy
  auto const vStarLhs = std::tuple{
      AssemblyScalarMass{1. / dt, feSpaceU},
      AssemblyAdvection{1.0, vel, feSpaceVel, feSpaceU},
      AssemblyStiffness{nu, feSpaceU},
  };
  auto const vStarRhs = std::tuple{
      AssemblyProjection{1. / dt, v.data, feSpaceU},
      AssemblyGradRhs2{1.0, pOld, feSpacePSplit, feSpaceU, {1}},
      AssemblyAnalyticRhs{g1, feSpaceU},
  };
  // AssemblyGradRhs{-1.0, pOld, feSpaceVel, feSpaceP});

  // dt \nabla^2 \delta p = \nabla \cdot velStar
  auto const pLhs = std::tuple{
      AssemblyStiffness{dt, feSpacePSplit},
  };
  auto const pRhs = std::tuple{
      AssemblyDivRhs{-1.0, velStar, feSpaceVel, feSpacePSplit},
  };
  // AssemblyStiffnessRhs{-dt, pOld, feSpaceP});

  // pOld += \delta p
  // u = uStar - dt d \delta p / dx
  auto const uLhs = std::tuple{
      AssemblyScalarMass{1.0, feSpaceU},
  };
  auto const uRhs = std::tuple{
      AssemblyProjection{1.0, uStar.data, feSpaceU},
      AssemblyGradRhs{-dt, p.data, feSpacePSplit, feSpaceU, {0}},
  };

  // v = vStar - dt d \delta p / dy
  auto const vLhs = std::tuple{
      AssemblyScalarMass{1.0, feSpaceU},
  };
  auto const vRhs = std::tuple{
      AssemblyProjection{1.0, vStar.data, feSpaceU},
      AssemblyGradRhs{-dt, p.data, feSpacePSplit, feSpaceU, {1}},
  };
  t.stop();

  t.start("eqn");
  Eqn eqnUstar{uStar, uStarLhs, uStarRhs, bcsUStar, EqnSolver<StorageType::ClmMajor>()};
  Eqn eqnVstar{vStar, vStarLhs, vStarRhs, bcsVStar, EqnSolver<StorageType::ClmMajor>()};

  Eqn eqnP{p, pLhs, pRhs, bcsPSplit};
  // eqnP lhs does not change in time, we can pre-compute and factorize it
  eqnP.buildLhs();
  eqnP.compute();

  Eqn eqnU{u, uLhs, uRhs, std::tuple{}};
  // eqnU lhs does not change in time, we can pre-compute and factorize it
  eqnU.buildLhs();
  eqnU.compute();
  Eqn eqnV{v, vLhs, vRhs, std::tuple{}};
  // eqnV lhs does not change in time, we can pre-compute and factorize it
  eqnV.buildLhs();
  eqnV.compute();
  t.stop();

  t.start("assembly RT");
  Var uStarRT{"uStarRT", feSpaceUStarRT.dof.size};
  Var vStarRT{"vStarRT", feSpaceUStarRT.dof.size};
  Vec velStarRT = Vec::Zero(feSpaceUStarRT.dof.size * Elem_T::dim);
  FEVar velRT{"velRT", feSpaceVelRT};
  Var lambda{"lambda", feSpaceLambda.dof.size};
  Vec lambdaOld = Vec::Zero(feSpaceLambda.dof.size);
  Vec dLambda = Vec::Zero(feSpaceLambda.dof.size);

  // 2 * uStar / dt + (velRT \cdot \nabla) uStar - \nabla \cdot (nu \nabla uStar) =
  // 2 * uRT / dt
  auto const uStarRTLhs = std::tuple{
      AssemblyScalarMass{2. / dt, feSpaceUStarRT},
      AssemblyAdvection{1.0, velRT.data, feSpaceVelRT, feSpaceUStarRT},
      AssemblyStiffness{nu, feSpaceUStarRT},
  };
  auto const uStarRTRhs = std::tuple{
      AssemblyProjection{2. / dt, velRT.data, feSpaceVelRT, feSpaceUStarRT, {0}},
      // AssemblyProjection{-1. / dt, uStarRT.data, feSpaceUStarRT},
      AssemblyGradRhs2{-1.0, lambdaOld, feSpaceLambda, feSpaceUStarRT, {0}},
      AssemblyAnalyticRhs{g0, feSpaceUStarRT},
  };

  // 2 * vStar / dt + (velRT \cdot \nabla) vStar - \nabla \cdot (nu \nabla vStar) =
  // 2 * vRT / dt
  auto const vStarRTLhs = std::tuple{
      AssemblyScalarMass{2. / dt, feSpaceUStarRT},
      AssemblyAdvection{1.0, velRT.data, feSpaceVelRT, feSpaceUStarRT},
      AssemblyStiffness{nu, feSpaceUStarRT},
  };
  auto const vStarRTRhs = std::tuple{
      AssemblyProjection{2. / dt, velRT.data, feSpaceVelRT, feSpaceUStarRT, {1}},
      // AssemblyProjection{-1. / dt, vStarRT.data, feSpaceUStarRT},
      AssemblyGradRhs2{-1.0, lambdaOld, feSpaceLambda, feSpaceUStarRT, {1}},
      AssemblyAnalyticRhs{g1, feSpaceUStarRT},
  };
  t.stop();

  t.start("eqn RT");
  Eqn eqnUstarRT{
      uStarRT, uStarRTLhs, uStarRTRhs, bcsUStarRT, EqnSolver<StorageType::RowMajor>()};
  Eqn eqnVstarRT{
      vStarRT, vStarRTLhs, vStarRTRhs, bcsVStarRT, EqnSolver<StorageType::RowMajor>()};

  // \vec{uRT} - dt \nabla \lambda = \vec{uStar}
  // \nabla \cdot \vec{uRT} = 0
  Builder<StorageType::RowMajor> builderRT{
      feSpaceVelRT.dof.size + feSpaceLambda.dof.size};
  builderRT.buildLhs(std::tuple{AssemblyVectorMass{1. / dt, feSpaceVelRT}}, bcsVelRT);
  builderRT.buildCoupling(
      AssemblyVectorGrad(1.0, feSpaceVelRT, feSpaceLambdaCoupling),
      bcsVelRT,
      bcsLambda);
  builderRT.buildCoupling(
      AssemblyVectorDiv(1.0, feSpaceLambdaCoupling, feSpaceVelRT), bcsLambda, bcsVelRT);
  auto const velRTRhs = std::tuple{
      AssemblyProjection{1. / dt, velStarRT, feSpaceVelStarRT, feSpaceVelRT}};
  builderRT.closeMatrix();
  auto const velRTRhsFixed = builderRT.b;
  builderRT.clearRhs();
  // IterSolver solverRT;
  LUSolver solverRT;
  solverRT.compute(builderRT.A);
  t.stop();

  t.start("ic");
  auto const ic = [](Vec3 const &) { return Vec2(0., 0.); };
  // auto const ic = [](Vec3 const & p) { return Vec2{0., 1.5 * (1. - p(0)*p(0))}; };
  // auto const ic = zero;

  // monolithic
  Var velM{"velM"};
  velM.data = Vec::Zero(dofU * FESpaceVel_T::dim + dofP);
  interpolateAnalyticFunction(ic, feSpaceVel, velM.data);

  // split
  auto const ic0 = [&ic](Vec3 const & p) { return ic(p)[0]; };
  auto const ic1 = [&ic](Vec3 const & p) { return ic(p)[1]; };
  interpolateAnalyticFunction(ic0, feSpaceU, u.data);
  interpolateAnalyticFunction(ic1, feSpaceU, v.data);
  uStar.data = u.data;
  vStar.data = v.data;

  // split RT
  interpolateAnalyticFunction(ic0, feSpaceUStarRT, uStarRT.data);
  interpolateAnalyticFunction(ic1, feSpaceUStarRT, vStarRT.data);
  auto const ic3d = [&ic](Vec3 const & p) { return promote<3>(ic(p)); };
  interpolateAnalyticFunction(ic3d, feSpaceVelRT, velRT.data);
  t.stop();

  t.start("print");
  IOManager ioVelM{feSpaceVel, "output_ipcs_rt/velm"};
  IOManager ioPM{feSpaceP, "output_ipcs_rt/pm"};
  Var pM{"pM", velM.data.block(dofU * FESpaceVel_T::dim, 0, dofP, 1)};
  if (config["monolithic"].as<bool>())
  {
    ioVelM.print({velM});
    ioPM.print({pM});
  }

  IOManager ioUsplit{feSpaceU, "output_ipcs_rt/usplit"};
  IOManager ioPsplit{feSpaceP, "output_ipcs_rt/psplit"};
  Var pSplit{"p", feSpacePSplit.dof.size};
  if (config["split"].as<bool>())
  {
    ioUsplit.print({eqnUstar.sol, eqnVstar.sol, eqnU.sol, eqnV.sol});
    ioPsplit.print({pSplit});
  }

  IOManager ioUsplitRT{feSpaceUStarRT, "output_ipcs_rt/usplit_rt"};
  IOManagerP0 ioVelRT{feSpaceVelRT, "output_ipcs_rt/velrt"};
  IOManagerFacet ioFacet{feSpaceVelRT, "output_ipcs_rt/facet"};
  IOManager ioLambda{feSpaceLambda, "output_ipcs_rt/lambda"};
  if (config["split_rt"].as<bool>())
  {
    ioUsplitRT.print({eqnUstarRT.sol, eqnVstarRT.sol});
    ioVelRT.print(std::tuple{velRT});
    ioFacet.print(std::tuple{velRT});
    ioLambda.print({lambda});
  }
  t.stop();

  MilliTimer timerStep;
  double time = 0.0;
  auto const ntime = config["ntime"].as<uint>();
  auto const printStep = config["printStep"].as<uint>();
  for (uint itime = 0; itime < ntime; itime++)
  {
    timerStep.start();
    time += dt;
    std::cout << "\n"
              << Utils::separator << "solving timestep " << itime + 1
              << ", time = " << time << std::endl;
    // filelog << "\n" << Utils::separator;

    if (config["monolithic"].as<bool>())
    {
      t.start("build monolithic");
      velOldMonolithic = velM.data;

      builderM.buildLhs(std::tuple{timeder, advection, stiffness}, bcsVel);
      builderM.buildCoupling(grad, bcsVel, bcsP);
      builderM.buildCoupling(div, bcsP, bcsVel);
      builderM.buildRhs(std::tuple{timederRhs}, bcsVel);
      builderM.closeMatrix();
      t.stop();

      t.start("solve monolithic");
      solverM.compute(builderM.A);
      velM.data = solverM.solve(builderM.b);
      auto const res = builderM.A * velM.data - builderM.b;
      std::cout << "residual norm: " << res.norm() << std::endl;
      t.stop();

      t.start("clear monolithic");
      builderM.clear();
      t.stop();
    }

    if (config["split"].as<bool>())
    {
      t.start("build ustar");
      setComponent(vel, feSpaceVel, u.data, feSpaceU, 0);
      setComponent(vel, feSpaceVel, v.data, feSpaceU, 1);
      pOld += eqnP.sol.data;
      eqnUstar.build();
      eqnVstar.build();
      t.stop();

      t.start("solve ustar");
      eqnUstar.compute();
      eqnUstar.solve();
      std::cout << "eqnUstar residual norm: " << eqnUstar.residualNorm() << std::endl;
      eqnVstar.compute();
      eqnVstar.solve();
      std::cout << "eqnVstar residual norm: " << eqnVstar.residualNorm() << std::endl;
      t.stop();

      t.start("build p");
      setComponent(velStar, feSpaceVel, uStar.data, feSpaceU, 0);
      setComponent(velStar, feSpaceVel, vStar.data, feSpaceU, 1);
      eqnP.buildRhs();
      t.stop();

      t.start("solve p");
      eqnP.solve();
      std::cout << "eqnP residual norm: " << eqnP.residualNorm() << std::endl;
      t.stop();

      t.start("build u");
      eqnU.buildRhs();
      eqnV.buildRhs();
      t.stop();

      t.start("solve u");
      eqnU.solve();
      std::cout << "eqnU residual norm: " << eqnU.residualNorm() << std::endl;
      eqnV.solve();
      std::cout << "eqnV residual norm: " << eqnV.residualNorm() << std::endl;
      t.stop();

      t.start("clear split");
      eqnUstar.builder.clear();
      eqnVstar.builder.clear();
      eqnP.builder.clear();
      eqnU.builder.clear();
      eqnV.builder.clear();
      t.stop();
    }

    if (config["split_rt"].as<bool>())
    {
      t.start("build ustarRT");
      lambdaOld += dLambda;
      eqnUstarRT.build();
      eqnVstarRT.build();
      t.stop();

      t.start("solve ustarRT");
      eqnUstarRT.compute();
      eqnUstarRT.solve();
      std::cout << "eqnUstarRT residual norm: " << eqnUstar.residualNorm() << std::endl;
      eqnVstarRT.compute();
      eqnVstarRT.solve();
      std::cout << "eqnVstarRT residual norm: " << eqnVstar.residualNorm() << std::endl;
      t.stop();

      t.start("build proj RT");
      setComponent(velStarRT, feSpaceVelStarRT, uStarRT.data, feSpaceUStarRT, 0);
      setComponent(velStarRT, feSpaceVelStarRT, vStarRT.data, feSpaceUStarRT, 1);
      builderRT.buildRhs(velRTRhs, bcsVelRT);
      builderRT.b += velRTRhsFixed;
      std::cout << "RT rhs norm: " << builderRT.b.norm() << std::endl;
      t.stop();

      t.start("solve RT");
      // std::cout << "ART:\n" << builderRT.A << std::endl;
      auto const solRT = solverRT.solve(builderRT.b);
      velRT.data = solRT.block(0, 0, feSpaceVelRT.dof.size, 1);
      dLambda = solRT.block(feSpaceVelRT.dof.size, 0, feSpaceLambda.dof.size, 1);
      std::cout << "RT residual norm: " << (builderRT.A * solRT - builderRT.b).norm()
                << std::endl;
      t.stop();

      t.start("clear RT");
      eqnUstarRT.builder.clear();
      eqnVstarRT.builder.clear();
      builderRT.clearRhs();
      t.stop();
    }

    t.start("print");
    if (config["split"].as<bool>())
    {
      pSplit.data += eqnP.sol.data;
    }
    if (config["split"].as<bool>())
    {
      lambda.data += dLambda;
    }
    if ((itime + 1) % printStep == 0)
    {
      std::cout << "printing" << std::endl;
      if (config["monolithic"].as<bool>())
      {
        ioVelM.print({velM}, time);
        pM.data = velM.data.block(dofU * FESpaceVel_T::dim, 0, dofP, 1);
        ioPM.print({pM}, time);
      }

      if (config["split"].as<bool>())
      {
        ioUsplit.print({eqnUstar.sol, eqnVstar.sol, eqnU.sol, eqnV.sol}, time);
        ioPsplit.print({pSplit}, time);
      }

      if (config["split_rt"].as<bool>())
      {
        ioUsplitRT.print({eqnUstarRT.sol, eqnVstarRT.sol}, time);
        ioVelRT.print(std::tuple{velRT}, time);
        ioFacet.print(std::tuple{velRT}, time);
        ioLambda.print({lambda}, time);
      }
    }
    t.stop();

    std::cout << "time required: " << timerStep << " ms" << std::endl;
  }

  t.start("boundary postpro");
  using MeshFacet_T = Mesh<typename Elem_T::Facet_T>;
  using FESpaceFacet_T = FESpace<
      MeshFacet_T,
      typename LagrangeFE<typename Elem_T::Facet_T, 0>::RefFE_T,
      typename LagrangeFE<typename Elem_T::Facet_T, 0>::RecommendedQR>;
  std::unique_ptr<MeshFacet_T> facetMesh{new MeshFacet_T};
  buildFacetMesh(*facetMesh, *mesh);

  FESpaceFacet_T feSpaceBd{*facetMesh};
  FEVar bd{"bd", feSpaceBd};
  interpolateOnFacets<Component::TANGENTIAL>(
      bd.data, feSpaceBd, velRT.data, feSpaceVelRT);
  if (bd.data.size() < 100)
    std::cout << bd.data.transpose() << std::endl;
  t.stop();

  t.print();

  return 1;
}
