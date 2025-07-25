#include "def.hpp"

// stl
#include <bitset>

// local
#include "mesh.hpp"

using namespace proxpde;

int vectorTest(Vec3 const & v, Vec3 const & e)
{
  auto const norm = (v - e).norm();
  if (norm > 1e-8)
  {
    std::cerr << "the normal " << v.transpose() << " does not correspond to "
              << v.transpose() << " with an error of " << norm << std::endl;
    return 1;
  }
  return 0;
}

int main()
{
  std::shared_ptr<Mesh<Triangle>> triangleMesh{new Mesh<Triangle>()};
  refTriangleMesh(*triangleMesh);
  buildNormals(*triangleMesh, MeshFlags::NONE);

  std::bitset<3> tests;
  tests[0] = vectorTest(triangleMesh->facetList[0]._normal, Vec3(0.0, -1.0, 0.0));
  tests[1] = vectorTest(triangleMesh->facetList[1]._normal, Vec3(-1.0, 0.0, 0.0));
  tests[2] = vectorTest(
      triangleMesh->facetList[2]._normal,
      Vec3(0.5 * std::sqrt(2), 0.5 * std::sqrt(2), 0.0));

  return tests.any();
}
