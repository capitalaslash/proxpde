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
