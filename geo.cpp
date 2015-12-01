#include "geo.hpp"

void buildMesh1D(std::shared_ptr<Mesh<Line>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 uint const numPts)
{
  Vec3 const h = length / (numPts-1);
  meshPtr->pointList.reserve(numPts);
  for(uint p=0; p<numPts; ++p)
  {
    meshPtr->pointList.emplace_back(origin + p * h, p);
  }
  meshPtr->pointList[0].marker = side::LEFT;
  meshPtr->pointList[numPts-1].marker = side::RIGHT;

  uint const numElems = numPts-1;
  meshPtr->elementList.reserve(numElems);
  for(uint e=0; e<numElems; ++e)
  {
    Line line({&meshPtr->pointList[e], &meshPtr->pointList[e+1]}, e);
    meshPtr->elementList.push_back(std::move(line));
  }

  meshPtr->buildConnectivity();
}

void buildMesh2D(std::shared_ptr<Mesh<Triangle>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 std::array<uint, 2> const numPts)
{
  Vec3 const h = {length(0) / (numPts[0]-1.), length(1) / (numPts[1]-1.), 1.};
  meshPtr->pointList.reserve(numPts[0]*numPts[1]);
  for(uint j=0; j<numPts[1]; ++j)
  {
    for(uint i=0; i<numPts[0]; ++i)
    {
      Vec3 p(i * h(0), j * h(1), 0.0);
      meshPtr->pointList.emplace_back(origin + p, i + numPts[0]*j);
    }
  }

  for(uint s=0; s<numPts[0]; s++)
  {
    meshPtr->pointList[s].marker = side::BOTTOM;
    meshPtr->pointList[(numPts[1]-1)*numPts[0]+s].marker = side::TOP;
  }
  for(uint s=0; s<numPts[1]; s++)
  {
    meshPtr->pointList[0+numPts[0]*s].marker = side::LEFT;
    meshPtr->pointList[numPts[0]-1+numPts[0]*s].marker = side::RIGHT;
  }

  uint const numElems = (numPts[0]-1)*(numPts[1]-1);
  meshPtr->elementList.reserve(numElems);
  id_T counter = 0;
  for(uint j=0; j<numPts[1]-1; ++j)
    for(uint i=0; i<numPts[0]-1; ++i)
    {
      uint const base = i + j*numPts[0];
      Triangle tri_b({&meshPtr->pointList[base],
                      &meshPtr->pointList[base+1],
                      &meshPtr->pointList[base+numPts[0]+1]},
                     counter++);
      Triangle tri_t({&meshPtr->pointList[base],
                      &meshPtr->pointList[base+numPts[0]+1],
                      &meshPtr->pointList[base+numPts[0]]},
                     counter++);
      meshPtr->elementList.push_back(std::move(tri_b));
      meshPtr->elementList.push_back(std::move(tri_t));
  }

  meshPtr->buildConnectivity();

}
