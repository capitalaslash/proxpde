#include "def.hpp"

#include "assembly.hpp"
#include "bc.hpp"
#include "builder.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "iomanager.hpp"
#include "mesh.hpp"

int main(int argc, char * argv[])
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  scalarFun_T rhs = [](Vec3 const & p)
  { return (1. + 2 * M_PI * M_PI) * std::sin(M_PI * p(0)) * std::sin(M_PI * p(1)); };
  scalarFun_T exactSol = [](Vec3 const & p)
  { return std::sin(M_PI * p(0)) * std::sin(M_PI * p(1)); };

  // vectorFun_T mesh_mod = [] (Vec3 const &p)
  // {
  //   double const x = 1. * (1.-std::cos(p(0) * M_PI / (1.*2.)));
  //   return Vec3(x, p(1), std::sqrt(2.*x-x*x));
  // };
  // vectorFun_T mesh_mod_inv = [] (Vec3 const &p)
  // {
  //   double const x = 1.*2.*std::acos(1 - p(0)/1.) / M_PI;
  //   return Vec3(x, p(1), 1.-std::sqrt(1.-x*x));
  // };

  vectorFun_T mesh_mod = [](Vec3 const & p) { return p; };
  vectorFun_T mesh_mod_inv = [](Vec3 const & p) { return p; };
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  uint const numElemsX = (argc < 3) ? 10 : std::stoi(argv[1]);
  uint const numElemsY = (argc < 3) ? 10 : std::stoi(argv[2]);
  Vec3 const origin{0., 0., 0.};
  Vec3 const length{1., 1., 0.};
  buildHyperCube(*mesh, origin, length, {{numElemsX, numElemsY, 0}});

  // // rotation matrix
  // double thetay = M_PI / 4.;
  // double thetaz = M_PI / 3.;
  // FMat<3,3> Ry, Rz;
  // Ry << std::cos(thetay), 0.0, std::sin(thetay),
  //       0.0, 1.0, 0.0,
  //      -std::sin(thetay), 0.0, std::cos(thetay);
  // Rz << std::cos(thetaz), std::sin(thetaz), 0.0,
  //      -std::sin(thetaz), std::cos(thetaz), 0.0,
  //     0.0, 0.0, 1.0;
  // auto R = Ry * Rz;
  // auto Rt = R.transpose();
  //
  // // rotate mesh
  // for (auto & p: mesh->pointList)
  // {
  //   p.coord = R * p.coord;
  // }

  vectorFun_T meshMod = [](Vec3 const & p) { return p; };

  vectorFun_T meshModInv = [](Vec3 const & p) { return p; };

  // mesh modifier
  for (auto & p: mesh->pointList)
  {
    p.coord = meshMod(p.coord);
  }

  FESpace_T feSpace{*mesh};

  auto const zeroFun = [](Vec3 const &) { return 0.; };
  auto bcRight = BCEss{feSpace, side::RIGHT};
  bcRight << zeroFun;
  auto bcLeft = BCEss{feSpace, side::LEFT};
  bcLeft << zeroFun;
  auto bcBottom = BCEss{feSpace, side::BOTTOM};
  bcBottom << zeroFun;
  auto bcTop = BCEss{feSpace, side::TOP};
  bcTop << zeroFun;
  auto const bcs = std::vector{bcRight, bcLeft, bcBottom, bcTop};

  AssemblyStiffness stiffness{1.0, feSpace};
  AssemblyScalarMass mass{1.0, feSpace};
  // auto rotatedRhs = [&Rt] (Vec3 const& p) { return rhs(Rt * p); };
  auto modifiedRhs = [&rhs, &mesh_mod_inv](Vec3 const & p)
  { return rhs(mesh_mod_inv(p)); };
  AssemblyRhsAnalytic f{modifiedRhs, feSpace};
  // AssemblyRhsAnalytic f{rhs, feSpace};

  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{stiffness, mass}, bcs);
  builder.buildRhs(std::tuple{f}, bcs);
  builder.closeMatrix();

  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);

  Var exact{"exact"};
  // auto rotatedESol = [&Rt] (Vec3 const& p) { return exact_sol(Rt * p); };
  auto modifiedESol = [&exactSol, &mesh_mod_inv](Vec3 const & p)
  { return exactSol(mesh_mod_inv(p)); };
  interpolateAnalyticFunction(modifiedESol, feSpace, exact.data);
  // interpolateAnalyticFunction(exact_sol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  IOManager io{feSpace, "output_dr2dquad/sol"};
  io.print({sol, exact, error});

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(16) << norm
            << std::endl;
  return checkError({norm}, {0.04331597477422448});
}
