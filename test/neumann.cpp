#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T, 1>::RefFE_T,
                          FEType<Elem_T, 1>::RecommendedQR>;

static scalarFun_T rhs = [] (Vec3 const & p)
{
  return 6. * p(0);
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return 4. * p(0) - pow(p(0), 3);
};

int main(int argc, char* argv[])
{
  MilliTimer t;

  t.start();
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  uint const numElems = (argc < 3)? 10 : std::stoi(argv[1]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(*mesh, origin, length, {{numElems, numElems, 0}});
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace(*mesh);
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList bcs{feSpace};
  bcs.addBC(BCEss{feSpace, side::LEFT, [] (Vec3 const &) { return 0.; }});
  // bcs.addBC(BCNat<FESpace_T>{side::RIGHT, [] (Vec3 const &) { return 1.; }});
  std::cout << "bcs: " << t << " ms" << std::endl;

  t.start();
  auto const size = feSpace.dof.size;
  AssemblyStiffness stiffness(1.0, feSpace);
  Builder builder{size};
  builder.buildProblem(AssemblyStiffness{1.0, feSpace}, bcs);
  builder.buildProblem(AssemblyAnalyticRhs{rhs, feSpace}, bcs);
  builder.buildProblem(AssemblyBCNatural{[] (Vec3 const &) { return 1.; }, side::RIGHT, feSpace}, bcs);
  builder.closeMatrix();
  std::cout << "fe build: " << t << " ms" << std::endl;

  t.start();
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  std::cout << "solve: " << t << " ms" << std::endl;

  Var flux{"flux"};
  using RecFESpace_T = FESpace<Mesh_T,
                               FEType<Elem_T, 1>::RefFE_T,
                               FEType<Elem_T, 1>::ReconstructionQR, 2>;
  RecFESpace_T feSpaceRec{*mesh};
  reconstructGradient(sol.data, feSpaceRec, flux.data, {0, 1});
  DOFCoordSet fluxSet{
    feSpace,
    [](Vec const & p) { return std::fabs(p[0] - 1.0) < 1e-12; }
  };
  for (auto const id: fluxSet.ids)
  {
    auto const boundaryFlux = flux.data[id];
    if (std::fabs(boundaryFlux - 1.29) > 1.e-8)
    {
      std::cerr << "the flux on the boundary is not the expected value: "<< boundaryFlux << std::endl;
      return 2;
    }
  }

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  t.start();
  IOManager io{feSpace, "output_neumann/sol"};
  io.print({sol, exact, error});
  IOManager ioFlux{feSpaceRec, "output_neumann/flux"};
  ioFlux.print({flux});
  std::cout << "output: " << t << " ms" << std::endl;

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm << std::endl;
  if(std::fabs(norm - 6.800165872394435e-14) > 1.e-15)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }
  return 0;
}
