#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "var.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"
#include "ns.hpp"

int main(int argc, char* argv[])
{
  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpaceU_T = FESpace<Mesh_T,
                             FEType<Elem_T,2>::RefFE_T,
                             FEType<Elem_T,2>::RecommendedQR, 2>;
  using FESpaceP_T = FESpace<Mesh_T,
                             FEType<Elem_T,1>::RefFE_T,
                             FEType<Elem_T,2>::RecommendedQR>;

  MilliTimer t;

  t.start("mesh");
  uint const numElemsX = (argc < 4)? 6 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 4)? 6 : std::stoi(argv[2]);
  uint const numElemsR = (argc < 4)? 6 : std::stoi(argv[3]);

  Vec3 const origin{0., 0., 0.};
  double const radius = 1.;

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  buildCircleMesh(*mesh, origin, radius, {{numElemsX, numElemsY, numElemsR}});
  t.stop();
  // std::cout << *mesh << std::endl;

  t.start("feSpace");
  FESpaceU_T feSpaceU{*mesh};
  FESpaceP_T feSpaceP{*mesh, 2*feSpaceU.dof.size};
  t.stop();

  t.start("bcs");
  auto const zeroFun = [] (Vec3 const&) {return Vec2{0., 0.};};
  auto const bcsU = std::make_tuple(
        BCEss{feSpaceU, side::CIRCLE, zeroFun},
        BCEss{feSpaceU, side::TOP,    zeroFun, {0}},
        BCEss{feSpaceU, side::LEFT,   zeroFun, {1}});
  DOFCoordSet pinSet
  {
    feSpaceP,
    [](Vec3 const & p)
    {
      return std::fabs(p[0] - 1.0) < 1.e-6 &&
          std::fabs(p[1] - 0.0) < 1.e-6;
    }
  };
  auto const bcsP = std::make_tuple(
        BCEss{feSpaceP, pinSet.ids, [] (Vec3 const &) {return 0.;}});
  t.stop();

  t.start("assembly");
  double const nu = 0.1;
  AssemblyStiffness stiffness(nu, feSpaceU);
  AssemblyGrad grad(-1.0, feSpaceU, feSpaceP);
  AssemblyDiv div(-1.0, feSpaceP, feSpaceU);
  // this is required to properly apply the pinning on the pressure
  AssemblyMass dummy{0.0, feSpaceP};

  uint const numDOFs = 2*feSpaceU.dof.size + feSpaceP.dof.size;
  Builder builder{numDOFs};
  builder.buildLhs(stiffness, bcsU);
  builder.buildLhs(dummy, bcsP);
  builder.buildCoupling(grad, bcsU, bcsP);
  builder.buildCoupling(div, bcsP, bcsU);
  builder.buildRhs(AssemblyBCNatural{[](Vec3 const &){return Vec2{1., 0.};}, side::LEFT, feSpaceU}, bcsU);
  builder.closeMatrix();
  t.stop();

  t.start("solve");
  Var sol{"vel", numDOFs};
  LUSolver solver(builder.A);
  sol.data = solver.solve(builder.b);
  t.stop();

  // std::cout << "A:\n" << builder.A << std::endl;
  // std::cout << "b:\n" << builder.b << std::endl;
  // std::cout << "sol:\n" << sol.data << std::endl;

  t.start("postpro");
  Var p{"p", sol.data , 2*feSpaceU.dof.size, feSpaceP.dof.size};
  t.stop();

  // integral on boundary
  t.start("boundary int");
  double pIntegral = 0.;
  {
    using FacetFEP_T = FESpaceP_T::RefFE_T::FacetFE_T;
    using FacetQR_T = SideQR_T<FESpaceP_T::QR_T>;
    using facetCurFEP_T = CurFE<FacetFEP_T, FacetQR_T>;
    using FacetFESpaceP_T =
      FESpace<Mesh_T,
              FESpaceP_T::RefFE_T,
              SideGaussQR<Elem_T, FacetQR_T::numPts>>;

    facetCurFEP_T facetCurFEP;
    FacetFESpaceP_T facetFESpaceP{*mesh};
    FEVar pFacet{facetFESpaceP};
    pFacet.data = p.data;

    for (auto & facet: mesh->facetList)
    {
      if (facet.marker == side::LEFT)
      {
        // std::cout << "facet " << facet.id << std::endl;
        facetCurFEP.reinit(facet);
        auto elem = facet.facingElem[0].ptr;
        auto const side = facet.facingElem[0].side;
        pFacet.reinit(*elem);
        for (uint q=0; q<FacetQR_T::numPts; ++q)
        {
          double const pValue = pFacet.evaluate(side * FacetQR_T::numPts + q);
          pIntegral += facetCurFEP.JxW[q] * pValue;
        }
      }
    }
    std::cout << "integral of pressure on left face: " << std::setprecision(16) << pIntegral << std::endl;
  }
  t.stop();

  // wall shear stress on circle
  t.start("wss");
  FESpace<Mesh_T, FEType<Elem_T, 0>::RefFE_T, FEType<Elem_T, 0>::RecommendedQR> feSpaceWSS{*mesh};
  // mean cell value
  Var wss{"wss"};
  computeElemWSS(wss.data, feSpaceWSS, sol.data, feSpaceU, {side::CIRCLE}, nu);
  t.stop();
  std::cout << "wss min: " << std::setprecision(16) << wss.data.minCoeff() << std::endl;
  std::cout << "wss max: " << std::setprecision(16) << wss.data.maxCoeff() << std::endl;

  // Var exact{"exact"};
  // interpolateAnalyticFunction(exact_sol, feSpaceU, exact.data);
  // Var error{"e"};
  // error.data = sol.data /*- exact.data*/;

  t.start("print");
  IOManager ioU{feSpaceU, "output_stokes2dcircle/vel"};
  ioU.print({sol});
  IOManager ioP{feSpaceP, "output_stokes2dcircle/p"};
  ioP.print({p});
  IOManager ioWSS{feSpaceWSS, "output_stokes2dcircle/wss"};
  ioWSS.print({wss});
  t.stop();

  t.print();

  // double norm = error.data.norm();
  // std::cout << "the norm of the error is " << norm << std::endl;

  if (std::fabs(pIntegral - 0.9800943987494872) > 1.e-12 ||
      std::fabs(wss.data.maxCoeff() - 0.3068317649009756) > 1.e-12)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
