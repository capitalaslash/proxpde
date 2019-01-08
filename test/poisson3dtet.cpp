#include "def.hpp"
#include "mesh.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "bc.hpp"
#include "assembly.hpp"
#include "builder.hpp"
#include "iomanager.hpp"
#include "timer.hpp"

#include <yaml-cpp/yaml.h>

template <typename Mesh>
void readMesh(Mesh & mesh, YAML::Node const & config)
{
  auto const mesh_type = config["mesh_type"].as<std::string>();
  if (mesh_type == "structured")
  {
    array<uint,3> numPts =
    {
      config["nx"].as<uint>() + 1,
      config["ny"].as<uint>() + 1,
      config["nz"].as<uint>() + 1
    };
    Vec3 const origin{0., 0., 0.};
    Vec3 const length{1., 1., 1.};
    MeshBuilder<typename Mesh::Elem_T> meshBuilder;
    meshBuilder.build(mesh, origin, length, numPts);
  }
  else if (mesh_type == "msh")
  {
    readGMSH(mesh, config["mesh_file"].as<std::string>());
  }
}

using Elem_T = Tetrahedron;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<Mesh_T,
                          FEType<Elem_T,1>::RefFE_T,
                          FEType<Elem_T,1>::RecommendedQR>;

static scalarFun_T rhs = [] (Vec3 const& p)
{
  return 2.5*M_PI*M_PI*std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  // return 0.;
};
static scalarFun_T exactSol = [] (Vec3 const& p)
{
  return std::sin(0.5*M_PI*p(0))*std::sin(1.5*M_PI*p(1));
  // return 1.;
};

int main(int argc, char* argv[])
{
  MilliTimer t;

  auto const configFile = (argc > 1) ? argv[1] : "poisson3dtet.yaml";
  auto const config = YAML::LoadFile(configFile);

  t.start();
  std::unique_ptr<Mesh_T> mesh{new Mesh_T};
  readMesh(*mesh, config);
  std::cout << "mesh build: " << t << " ms" << std::endl;

  t.start();
  FESpace_T feSpace{*mesh};
  std::cout << "fespace: " << t << " ms" << std::endl;

  t.start();
  BCList bcs{feSpace};
  // face refs with z-axis that exits from the plane, x-axis towards the right
  bcs.addEssentialBC(side::LEFT, [] (Vec3 const&) {return 0.;});
  bcs.addEssentialBC(side::BOTTOM, [] (Vec3 const&) {return 0.;});
  std::cout << "bcs: " << t << " ms" << std::endl;

  t.start();
  MilliTimer tBuild;
  tBuild.start();
  Builder builder{feSpace.dof.size};
  std::cout << tBuild << std::endl;
  tBuild.start();
  builder.buildProblem(AssemblyStiffness(1.0, feSpace), bcs);
  std::cout << tBuild << std::endl;
  tBuild.start();
  // builder.buildProblem(AssemblyAnalyticRhs(rhs, feSpace), bcs);
  // using an interpolated rhs makes its quality independent of the chosen qr
  std::cout << tBuild << std::endl;
  tBuild.start();
  Vec rhsProj;
  interpolateAnalyticFunction(rhs, feSpace, rhsProj);
  std::cout << tBuild << std::endl;
  tBuild.start();
  builder.buildProblem(AssemblyProjection(1.0, rhsProj, feSpace), bcs);
  std::cout << tBuild << std::endl;
  tBuild.start();
  builder.closeMatrix();
  std::cout << tBuild << std::endl;
  tBuild.start();
  std::cout << "fe build: " << t << " ms" << std::endl;

  t.start();
  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);
  std::cout << "solve: " << t << " ms" << std::endl;

  t.start();
  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);
  Var error{"e"};
  error.data = sol.data - exact.data;

  IOManager io{feSpace, "output_poisson3dtet/sol"};
  io.print({sol, exact, error});
  std::cout << "output: " << t << " ms" << std::endl;

  double norm = error.data.norm();
  std::cout << "the norm of the error is " << std::setprecision(12) << norm << std::endl;
  if(std::fabs(norm - config["expected_error"].as<double>()) > 1.e-10)
  {
    std::cerr << "the norm of the error is not the prescribed value" << std::endl;
    return 1;
  }

  return 0;
}
