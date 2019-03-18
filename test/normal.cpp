#include "def.hpp"
#include "mesh.hpp"

int main()
{
  std::shared_ptr<Mesh<Triangle>> triangleMesh{new Mesh<Triangle>()};
  refTriangleMesh(*triangleMesh);
  buildFacets(*triangleMesh);
  buildNormals(*triangleMesh);
  auto normal0 = triangleMesh->facetList[0]._normal;
  assert((normal0 - Vec3(0.0, -1.0, 0.0)).norm() < 1e-12);
  auto normal1 = triangleMesh->facetList[1]._normal;
  assert((normal1 - Vec3(-1.0, 0.0, 0.0)).norm() < 1e-12);
  auto normal2 = triangleMesh->facetList[2]._normal;
  assert((normal2 - Vec3(.5 * std::sqrt(2.), .5 * std::sqrt(2.), 0.0)).norm() < 1e-12);
  return 0;
}
