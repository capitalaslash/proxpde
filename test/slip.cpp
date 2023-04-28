#include "def.hpp"

#include "assembler.hpp"
#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"
#include "ns.hpp"
#include "timer.hpp"
#include "var.hpp"

namespace proxpde
{

template <typename FESpace>
struct BCEssNormal: public BCEss<FESpace>
{
  using FESpace_T = FESpace;
  using RefFE_T = typename FESpace_T::RefFE_T;
  static constexpr uint dim = FESpace_T::dim;

  BCEssNormal(FESpace_T const & fe, marker_T const m): BCEss<FESpace_T>{fe, m} {}

  BCEssNormal<FESpace_T> operator<<(std::function<double(FVec<dim - 1> const &)> f)
  {
    auto const refPt = Vec3{0., 0., 0.};
    for (auto const & facet: this->feSpace->mesh->facetList)
    {
      if (facet.marker == this->marker)
      {
        auto const normal = narrow<dim>(facet.normal());
        auto const & [elem, side] = facet.facingElem[0];
        this->feSpace->curFE.reinit(*elem);
        if constexpr (
            order_v < RefFE_T >> 0 || family_v<RefFE_T> != FamilyType::LAGRANGE)
        {
          for (auto const dofFacet: RefFE_T::dofOnFacet[side])
          {
            auto const pt = this->feSpace->curFE.dofPts[dofFacet];
            // TODO: f(s) where s is the curvilinear variable(s) on the facet
            auto const value =
                -f(FVec<dim - 1>::Constant((pt - refPt).norm())) * normal;
            for (auto const c: this->comp)
            {
              DOFid_T const dofId = this->feSpace->dof.getId(elem->id, dofFacet, c);
              this->data[this->_constrainedDofMap.at(dofId)] = value[c];
              if (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
              {
                // value gives the entrant flux, it must be scaled to the facet size
                // ???
                this->data[this->_constrainedDofMap.at(dofId)] *= facet.volume();
              }
            }
          }
        }
        else // order_v<RefFE_T> == 0 && family_v<RefFE_T> == FamilyType::LAGRANGE
        {
          auto const value = -f(FVec<dim - 1>::Constant(0.)) * normal;
          for (auto const c: this->comp)
          {
            DOFid_T const dofId = this->feSpace.dof.getId(elem->id, 0, c);
            this->data[this->_constrainedDofMap.at(dofId)] = value[c];
          }
        }
      }
    }
    return *this;
  }

  BCEssNormal<FESpace_T> operator<<(std::function<double(double const &)> f)
  {
    static_assert(dim == 2, "scalar function on facet only available in 2D.");
    return this->operator<<([f](const Vec1 s) { return f(s[0]); });
  }
};

} // namespace proxpde

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using QuadraticRefFE = LagrangeFE<Elem_T, 2>::RefFE_T;
  using LinearRefFE = LagrangeFE<Elem_T, 1>::RefFE_T;
  using QuadraticQR = LagrangeFE<Elem_T, 2>::RecommendedQR;
  using FESpaceP_T = FESpace<Mesh_T, LinearRefFE, QuadraticQR>;
  using FESpaceVel_T = FESpace<Mesh_T, QuadraticRefFE, QuadraticQR, 2>;
  // using FESpaceComponent_T = FESpace<Mesh_T,QuadraticRefFE,QuadraticQR>;

  MilliTimer t;

  t.start("mesh");
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  uint const numElemsX = (argc < 3) ? 10 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3) ? 10 : std::stoi(argv[2]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0});
  // readGMSH(*mesh, "square_q.msh");
  t.stop();

  // rotated mesh
  // double theta = M_PI / 4.;
  // FMat<3,3> R;
  // R << std::cos(theta), std::sin(theta), 0.0,
  //     -std::sin(theta), std::cos(theta), 0.0,
  //     0.0, 0.0, 1.0;
  // auto Rt = R.transpose();
  // for (auto & p: mesh->pointList)
  // {
  //   p.coord = R * p.coord;
  // }

  // bend
  Vec3 const center = {2., 0., 0.};
  double const totalAngle = M_PI / 3.;
  auto const bend = [center](Vec3 const & p, double const angle)
  {
    double const r = center[0] - p[0];
    return Vec3{
        center[0] + r * std::cos(M_PI - angle),
        center[1] + r * std::sin(M_PI - angle),
        0.};
  };
  for (uint j = 0; j < numElemsY + 1; ++j)
  {
    double const angle = j * totalAngle / numElemsY;
    for (uint i = 0; i < numElemsX + 1; ++i)
    {
      auto & p = mesh->pointList[i + (numElemsX + 1) * j];
      p.coord = bend(p.coord, angle);
      // double const r = center[0] - p[0];
      // p.coord[0] = center[0] + r * std::cos(M_PI - angle);
      // p.coord[1] = center[1] + r * std::sin(M_PI - angle);
    }
  }

  t.start("fespace");
  FESpaceVel_T feSpaceVel{*mesh};
  // FESpaceComponent_T feSpaceComponent{*mesh};
  FESpaceP_T feSpaceP{*mesh, FESpaceVel_T::dim * feSpaceVel.dof.size};
  t.stop();

  t.start("bc");
  auto const zero = [](Vec3 const &) { return Vec2::Constant(0.); };
  // auto const one = [] (Vec3 const &) { return -1.; };
  // auto const oneY = [] (Point const &) { return Vec2{0.0, 1.0}; };
  // auto const inlet = [&R] (Vec3 const & p)
  // {
  //   return narrow<2>(R * Vec3{0., 5 * p[0] * (1. - p[0]), 0.});
  // };
  // auto const inlet = [&R] (Vec3 const & ) { return narrow<2>(R * Vec3{0., 1., 0.});
  // }; auto const inlet = [&R] (Vec3 const & p) {
  //   return narrow<2>(R * Vec3{0., (2. - p[0]), 0.});
  // };
  auto bcBottom = BCEssNormal{feSpaceVel, side::BOTTOM};
  // bcBottom << inlet;
  bcBottom << [](double const & /*s*/) { return 1. /*6. * s * (1. - s)*/; };
  auto const bcsVel = std::tuple{bcBottom};

  // auto const pinPt = bend(Vec3{0.5, 1., 0.}, totalAngle);
  // auto const toll = 1.e-12;
  // auto const close = [toll] (Vec3 const pt1, Vec3 const pt2)
  // {
  //   return (pt1 - pt2).norm() < toll;
  // };
  // DOFCoordSet pinSet{
  //   feSpaceP,
  //       [close, pinPt] (Vec3 const & p) { return close(p, pinPt); }
  // };
  // auto const bcPin = BCEss{feSpaceP, pinSet.ids};
  // auto const bcsP = std::tuple{bcPin};

  // auto bcPTop = BCEss(feSpaceP, side::TOP);
  // bcPTop << [] (Vec3 const &) { return 0.; };
  // auto const bcsP = std::tuple{bcPTop};

  auto const bcsP = std::tuple{};
  t.stop();

  auto const dofU = feSpaceVel.dof.size;
  auto const dofP = feSpaceP.dof.size;

  t.start("build");
  double const nu = 0.1;
  Vec velOld{2 * dofU};
  auto const stiffness = AssemblyTensorStiffness{nu, feSpaceVel};
  auto const advection = AssemblyAdvection{0.0, velOld, feSpaceVel, feSpaceVel};
  auto const grad = AssemblyGrad{-1., feSpaceVel, feSpaceP};
  auto const div = AssemblyDiv{-1., feSpaceP, feSpaceVel};
  auto const dummy = AssemblyDummy{feSpaceP};
  auto const bcNat = AssemblyBCNatural{zero, side::TOP, feSpaceVel};

  Builder builder{dofU * FESpaceVel_T::dim + dofP};
  builder.buildLhs(std::tuple{stiffness, advection}, bcsVel);
  builder.buildCoupling(grad, bcsVel, bcsP);
  builder.buildCoupling(div, bcsP, bcsVel);
  builder.buildLhs(std::tuple{dummy}, bcsP);
  builder.buildRhs(std::tuple{bcNat}, bcsVel);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol{"vel"};
  sol.data = Vec::Zero(2 * dofU + dofP);
  LUSolver solver(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol.data << std::endl;

  // Var exact{"exact", dofU*FESpaceVel_T::dim + dofP};
  // interpolateAnalyticFunction(inlet, feSpaceVel, exact.data);
  // interpolateAnalyticFunction(
  //       [] (Vec3 const & p) { return (1. - p[1]); },
  //       feSpaceP,
  //       exact.data);
  //
  // Var u{"u"};
  // Var v{"v"};
  // getComponent(u.data, feSpaceComponent, sol.data, feSpaceVel, 0);
  // getComponent(v.data, feSpaceComponent, sol.data, feSpaceVel, 1);
  Var p{"p", sol.data, 2 * dofU, dofP};
  //
  // Var ue{"ue"};
  // Var ve{"ve"};
  // getComponent(ue.data, feSpaceComponent, exact.data, feSpaceVel, 0);
  // getComponent(ve.data, feSpaceComponent, exact.data, feSpaceVel, 1);
  // Var pe{"pe", exact.data, 2*dofU, dofP};

  // wall shear stress
  t.start("wssCell");
  FESpace<Mesh_T, LagrangeFE<Elem_T, 0>::RefFE_T, LagrangeFE<Elem_T, 0>::RecommendedQR>
      feSpaceP0{*mesh};
  // mean cell value
  Var wssCell{"wssCell"};
  computeElemWSS(
      wssCell.data, feSpaceP0, sol.data, feSpaceVel, {side::RIGHT, side::LEFT}, nu);
  t.stop();
  std::cout << "wssCell min: " << std::setprecision(16) << wssCell.data.minCoeff()
            << std::endl;
  std::cout << "wssCell max: " << std::setprecision(16) << wssCell.data.maxCoeff()
            << std::endl;

  t.start("wssFacet");
  Var wssFacet{"wssFacet"};
  computeFEWSS(
      wssFacet.data, feSpaceP0, sol.data, feSpaceVel, {side::RIGHT, side::LEFT}, nu);
  t.stop();
  std::cout << "wssFacet min: " << std::setprecision(16) << wssFacet.data.minCoeff()
            << std::endl;
  std::cout << "wssFacet max: " << std::setprecision(16) << wssFacet.data.maxCoeff()
            << std::endl;

  t.start("print");
  IOManager ioVel{feSpaceVel, "output_slip/vel"};
  ioVel.print({sol, /*exact*/});
  IOManager ioP{feSpaceP, "output_slip/p"};
  ioP.print({p, /*pe*/});
  IOManager ioWSS{feSpaceP0, "output_slip/wss"};
  ioWSS.print({wssCell, wssFacet});
  t.stop();

  auto bcLeft = BCEss{feSpaceVel, side::LEFT};
  auto bcRight = BCEss{feSpaceVel, side::RIGHT};
  double normOld = 0.;
  for (uint k = 0; k < 20; ++k)
  {
    std::cout << Utils::separator << "iteration " << k << std::endl;

    t.start("tangent");
    bcLeft << sol.data;
    bcLeft.makeTangent();
    // std::cout << "bcLeft: " << bcLeft.data.transpose() << std::endl;
    bcRight << sol.data;
    bcRight.makeTangent();
    auto const bcsVelNew = std::tuple_cat(bcsVel, std::tuple{bcLeft, bcRight});
    t.stop();

    t.start("build");
    builder.clear();
    velOld = sol.data.head(2 * dofU);
    builder.buildLhs(std::tuple{stiffness, advection}, bcsVelNew);
    builder.buildCoupling(grad, bcsVelNew, bcsP);
    builder.buildCoupling(div, bcsP, bcsVelNew);
    builder.buildLhs(std::tuple{dummy}, bcsP);
    builder.buildRhs(std::tuple{bcNat}, bcsVel);
    builder.closeMatrix();
    t.stop();

    t.start("solve");
    solver.compute(builder.A);
    sol.data = solver.solve(builder.b);
    auto const res = builder.A * sol.data - builder.b;
    std::cout << "residual norm: " << res.norm() << std::endl;
    p.data = sol.data.tail(dofP);
    t.stop();

    t.start("integral");
    double pIntegral = integrateOnBoundary(p.data, feSpaceP, side::TOP);
    t.stop();
    std::cout << "integral of pressure on top face: " << std::setprecision(16)
              << pIntegral << std::endl;

    t.start("wss");
    computeFEWSS(
        wssFacet.data, feSpaceP0, sol.data, feSpaceVel, {side::RIGHT, side::LEFT}, nu);
    t.stop();
    std::cout << "wssFacet min: " << std::setprecision(16) << wssFacet.data.minCoeff()
              << std::endl;
    std::cout << "wssFacet max: " << std::setprecision(16) << wssFacet.data.maxCoeff()
              << std::endl;

    t.start("print");
    ioVel.print({sol}, k + 1);
    ioP.print({p}, k + 1);
    ioWSS.print({wssFacet}, k + 1);
    t.stop();

    double const norm = sol.data.norm();
    std::cout << "solution norm: " << norm << std::endl;
    if (std::fabs(norm - normOld) < 1.e-5)
    {
      break;
    }
    normOld = norm;
  }

  t.print();

  // auto uError = (u.data - ue.data).norm();
  // auto vError = (v.data - ve.data).norm();
  // auto pError = (p.data - pe.data).norm();
  //
  // std::cout << "u error norm: " << uError << std::endl;
  // std::cout << "v error norm: " << vError << std::endl;
  // std::cout << "p error norm: " << pError << std::endl;
  //
  // if (std::fabs(uError) > 1.e-14 ||
  //     std::fabs(vError) > 2.e-14 ||
  //     std::fabs(pError) > 3.e-14 ||
  //     std::fabs(wssCell.data.maxCoeff() - (nu * (1. - .5 / numElemsY))) > 1.e-12 ||
  //     std::fabs(wssFacet.data.maxCoeff() - nu) > 1.e-12)
  // {
  //   std::cerr << "one of the norms of the error is not the prescribed value" <<
  //   std::endl; return 1;
  // }

  return 0;
}
