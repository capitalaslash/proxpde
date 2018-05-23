#include "mesh.hpp"

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

  uint const numElems = numPts-1;
  meshPtr->elementList.reserve(numElems);
  for(uint e=0; e<numElems; ++e)
  {
    Line line({&meshPtr->pointList[e], &meshPtr->pointList[e+1]}, e);
    meshPtr->elementList.push_back(std::move(line));
  }

  meshPtr->buildConnectivity();
  buildFacets(meshPtr);
  meshPtr->facetList[0].marker = side::LEFT;
  meshPtr->facetList[1].marker = side::RIGHT;
}

template <typename Mesh>
void markFacetsCube(std::shared_ptr<Mesh> meshPtr,
                    Vec3 const& origin,
                    Vec3 const& length)
{
  for(auto & f: meshPtr->facetList)
  {
    auto const [min, max] = f.bbox();
    if(std::fabs(max[1]-origin[1]) < 1e-6*length[1])
    {
      f.marker = side::BOTTOM;
    }
    else if(std::fabs(min[0]-origin[0]-length[0]) < 1e-6*length[0])
    {
      f.marker = side::RIGHT;
    }
    else if(std::fabs(min[1]-origin[1]-length[1]) < 1e-6*length[1])
    {
      f.marker = side::TOP;
    }
    else if(std::fabs(max[0]-origin[0]) < 1e-6*length[0])
    {
      f.marker = side::LEFT;
    }
    else if(std::fabs(max[2]-origin[2]) < 1e-6*length[2])
    {
      f.marker = side::BACK;
    }
    else if(std::fabs(min[2]-origin[2]-length[2]) < 1e-6*length[2])
    {
      f.marker = side::FRONT;
    }
  }
}

void buildMesh2D(std::shared_ptr<Mesh<Triangle>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 2> const numPts)
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
      id_T const base = i + j*numPts[0];
      array<id_T,3> triplet_b, triplet_t;
      // TODO: make the traingle pattern configurable via a templated flag
      // if ((i-0.5*(numPts[0]-2))*(j-0.5*(numPts[1]-2)) < 0)
      if (base % 2)
      {
        triplet_b[0] = base;
        triplet_b[1] = base+1;
        triplet_b[2] = base+numPts[0];
        triplet_t[0] = base+1;
        triplet_t[1] = base+numPts[0]+1;
        triplet_t[2] = base+numPts[0];
      }
      else
      {
        triplet_b[0] = base;
        triplet_b[1] = base+1;
        triplet_b[2] = base+numPts[0]+1;
        triplet_t[0] = base;
        triplet_t[1] = base+numPts[0]+1;
        triplet_t[2] = base+numPts[0];
      }
      meshPtr->elementList.emplace_back(
        Triangle{{&meshPtr->pointList[triplet_b[0]],
                  &meshPtr->pointList[triplet_b[1]],
                  &meshPtr->pointList[triplet_b[2]]},
                 counter++});
      meshPtr->elementList.emplace_back(
        Triangle{{&meshPtr->pointList[triplet_t[0]],
                  &meshPtr->pointList[triplet_t[1]],
                  &meshPtr->pointList[triplet_t[2]]},
                 counter++});
    }
  meshPtr->buildConnectivity();
  buildFacets(meshPtr);
  markFacetsCube(meshPtr, origin, length);
}

void buildMesh2D(std::shared_ptr<Mesh<Quad>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 2> const numPts)
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
      id_T const base = i + j*numPts[0];
      meshPtr->elementList.emplace_back(
        Quad{{&meshPtr->pointList[base],
              &meshPtr->pointList[base+1],
              &meshPtr->pointList[base+numPts[0]+1],
              &meshPtr->pointList[base+numPts[0]]},
             counter++});
    }
  meshPtr->buildConnectivity();
  buildFacets(meshPtr);
  markFacetsCube(meshPtr, origin, length);
}

void buildCircleMesh(std::shared_ptr<Mesh<Quad>> meshPtr,
                     Vec3 const& origin,
                     double const& radius,
                     array<uint, 3> const numPts)
{
  meshPtr->pointList.resize(
    numPts[0]*numPts[1] +
    numPts[0]*(numPts[2]-1) +
    (numPts[1]-1)*(numPts[2]-1)
  );
  // uint ptCounter = 0;
  Vec3 r0(0.6*radius, 0.0, 0.0);
  Vec3 ra(sqrt(2.)*0.3*radius, -sqrt(2.)*0.3*radius, 0.0);
  Vec3 l0(0.0, 0.0, 0.0);
  Vec3 la(0.0, -0.6*radius, 0.0);
  for(uint j=0; j<numPts[1]; ++j)
  {
    double sj = j/(numPts[1]-1.);
    Vec3 b = (1.-sj)*l0 + sj*la;
    Vec3 e = (1.-sj)*r0 + sj*ra;
    for(uint i=0; i<numPts[0]; ++i)
    {
      double si = i/(numPts[0]-1.);
      Vec3 p = (1.-si)*b + si*e;
      const id_T id = i + numPts[0]*j;
      meshPtr->pointList[id] = Point(origin + p, id);
    }
  }

  uint offset = numPts[0]*numPts[1];
  for(uint i=0; i<numPts[0]; ++i)
  {
    double const alpha = 0.25*M_PI*i/(numPts[0]-1.);
    double si = i/(numPts[0]-1.);
    Vec3 b = la*(1.-si) + ra*si;
    Vec3 e(std::sin(alpha) * radius, -std::cos(alpha) * radius, 0.0);
    for(uint j=1; j<numPts[2]; ++j)
    {
      double sj = j/(numPts[2]-1.);
      Vec3 const p = b*(1-sj) + e*sj;
      const id_T id = offset + i + numPts[0]*(j-1);
      meshPtr->pointList[id] = Point(origin + p, id);
    }
  }

  offset += numPts[0]*(numPts[2]-1);
  for(uint j=0; j<numPts[1]-1; ++j)
  {
    double const alpha = 0.25*M_PI*j/(numPts[1]-1);
    double const sj = j/(numPts[1]-1.);
    Vec3 b = r0*(1.-sj)+ra*sj;
    Vec3 e(std::cos(alpha), -std::sin(alpha), 0.0);
    for(uint i=1; i<numPts[2]; ++i)
    {
      double const si = i/(numPts[2]-1.);
      Vec3 p = b*(1.-si) + e*si;
      const id_T id = offset + i-1 + (numPts[2]-1)*j;
      meshPtr->pointList[id] = Point(origin + p, id);
    }
  }

  // for(uint s=0; s<numPts[0]; s++)
  // {
  //   meshPtr->pointList[s].marker = side::TOP;
  // }
  // for(uint s=0; s<numPts[1]; s++)
  // {
  //   meshPtr->pointList[0+numPts[0]*s].marker = side::LEFT;
  //   meshPtr->pointList[numPts[0]-1+numPts[0]*s].marker = side::RIGHT;
  // }

  uint const numElems =
    (numPts[0]-1)*(numPts[1]-1) +
    (numPts[0]-1)*(numPts[2]-1) +
    (numPts[1]-1)*(numPts[2]-1);

  meshPtr->elementList.reserve(numElems);
  id_T counter = 0;
  for(uint j=0; j<numPts[1]-1 + numPts[2]-1; ++j)
    for(uint i=0; i<numPts[0]-1; ++i)
    {
      id_T const base = i + j*numPts[0];
      meshPtr->elementList.emplace_back(
        Quad{{&meshPtr->pointList[base],
              &meshPtr->pointList[base+1],
              &meshPtr->pointList[base+numPts[0]+1],
              &meshPtr->pointList[base+numPts[0]]},
             counter++});
    }

  uint blockOffset = numPts[0]*numPts[1]+numPts[0]*(numPts[2]-1);
  for(uint j=0; j<numPts[1]-2; ++j)
  {
    id_T const base = numPts[0]-1 + j*numPts[0];
    meshPtr->elementList.emplace_back(
      Quad{{&meshPtr->pointList[base],
            &meshPtr->pointList[blockOffset+j*(numPts[2]-1)],
            &meshPtr->pointList[blockOffset+(j+1)*(numPts[2]-1)],
            &meshPtr->pointList[base+numPts[0]]},
            counter++});
  }
  meshPtr->elementList.emplace_back(
    Quad{{&meshPtr->pointList[numPts[0]-1 + (numPts[1]-2)*numPts[0]],
          &meshPtr->pointList[blockOffset+(numPts[1]-2)*(numPts[2]-1)],
          &meshPtr->pointList[numPts[0]-1 + (numPts[1]-2)*numPts[0]+numPts[0] + numPts[0]],
          &meshPtr->pointList[numPts[0]-1 + (numPts[1]-2)*numPts[0]+numPts[0]]},
          counter++});

  for(uint i=1; i<numPts[2]-1; ++i)
  {
    for(uint j=0; j<numPts[1]-2; ++j)
    {
      id_T const base = blockOffset + i-1 + j*(numPts[2]-1);
      meshPtr->elementList.emplace_back(
        Quad{{&meshPtr->pointList[base],
              &meshPtr->pointList[base+1],
              &meshPtr->pointList[base+numPts[2]-1+1],
              &meshPtr->pointList[base+numPts[2]-1]},
              counter++});
    }
    id_T const base = blockOffset + i-1 +(numPts[1]-2)*(numPts[2]-1);
    meshPtr->elementList.emplace_back(
      Quad{{&meshPtr->pointList[base],
            &meshPtr->pointList[base+1],
            &meshPtr->pointList[numPts[0]*numPts[1]-1+numPts[0]*(i+1)],
            &meshPtr->pointList[numPts[0]*numPts[1]-1+numPts[0]*i]},
            counter++});
  }

  meshPtr->buildConnectivity();
  buildFacets(meshPtr);
  for(auto & f: meshPtr->facetList)
  {
    if(std::fabs(f.pointList[0]->coord[1]-origin[1]) < 1e-6*radius &&
       std::fabs(f.pointList[1]->coord[1]-origin[1]) < 1e-6*radius)
    {
      f.marker = side::TOP;
    }
    else if(std::fabs(f.pointList[0]->coord[0]-origin[0]) < 1e-6*radius &&
            std::fabs(f.pointList[1]->coord[0]-origin[0]) < 1e-6*radius)
    {
      f.marker = side::LEFT;
    }
    else if(std::fabs((f.pointList[0]->coord-origin).norm()-radius) < 1e-6*radius &&
            std::fabs((f.pointList[1]->coord-origin).norm()-radius) < 1e-6*radius)
    {
      f.marker = side::CIRCLE;
    }
  }
}

void buildMesh3D(std::shared_ptr<Mesh<Tetrahedron>> meshPtr,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 3> const numPts)
{
  Vec3 const h = {length(0) / (numPts[0]-1.), length(1) / (numPts[1]-1.), length(2) / (numPts[2]-1.)};
  meshPtr->pointList.reserve(numPts[0]*numPts[1]*numPts[2]);
  for(uint k=0; k<numPts[2]; ++k)
  {
    for(uint j=0; j<numPts[1]; ++j)
    {
      for(uint i=0; i<numPts[0]; ++i)
      {
        Vec3 p(i * h(0), j * h(1), k * h(2));
        auto const id = i + numPts[0]*j + numPts[0]*numPts[1]*k;
        meshPtr->pointList.emplace_back(origin + p, id);
      }
    }
  }

  for(uint si=0; si<numPts[0]; si++)
    for(uint sj=0; sj<numPts[1]; sj++)
    {
      auto const sb = si + numPts[0]*sj;
      meshPtr->pointList[sb].marker = side::BACK;
      auto const sf = si + numPts[0]*sj + numPts[1]*numPts[0]*(numPts[2]-1);
      meshPtr->pointList[sf].marker = side::FRONT;
    }
  for(uint si=0; si<numPts[0]; si++)
    for(uint sk=0; sk<numPts[2]; sk++)
    {
      auto const sb = si + numPts[1]*numPts[0]*sk;
      meshPtr->pointList[sb].marker = side::BOTTOM;
      auto const st = si + numPts[0]*(numPts[1]-1) + numPts[1]*numPts[0]*sk;
      meshPtr->pointList[st].marker = side::TOP;
    }
  for(uint sj=0; sj<numPts[1]; sj++)
    for(uint sk=0; sk<numPts[2]; sk++)
    {
      auto const sl = numPts[0]*sj + numPts[1]*numPts[0]*sk;
      meshPtr->pointList[sl].marker = side::LEFT;
      auto const sr = numPts[0]-1 + numPts[0]*sj + numPts[1]*numPts[0]*sk;
      meshPtr->pointList[sr].marker = side::RIGHT;
    }

  uint const numElems = 6*(numPts[0]-1)*(numPts[1]-1)*(numPts[2]-1);
  meshPtr->elementList.reserve(numElems);
  id_T counter = 0;
  auto const dx = 1;
  auto const dy = numPts[0];
  auto const dz = numPts[1]*numPts[0];
  for(uint k=0; k<numPts[2]-1; ++k)
    for(uint j=0; j<numPts[1]-1; ++j)
      for(uint i=0; i<numPts[0]-1; ++i)
      {
        id_T const base = i + j*numPts[0] + k*numPts[1]*numPts[0];
        meshPtr->elementList.emplace_back(
              Tetrahedron{{
                            &meshPtr->pointList[base], // 0
                            &meshPtr->pointList[base + dx], // 1
                            &meshPtr->pointList[base + dy], // 3
                            &meshPtr->pointList[base + dz]}, // 4
             counter++});
        meshPtr->elementList.emplace_back(
              Tetrahedron{{
                            &meshPtr->pointList[base + dx], // 1
                            &meshPtr->pointList[base + dx + dy], // 2
                            &meshPtr->pointList[base + dz], // 4
                            &meshPtr->pointList[base + dx + dz]}, // 5
             counter++});
        meshPtr->elementList.emplace_back(
              Tetrahedron{{
                            &meshPtr->pointList[base + dx], // 1
                            &meshPtr->pointList[base + dx + dy], // 2
                            &meshPtr->pointList[base + dy], // 3
                            &meshPtr->pointList[base + dz]}, // 4
             counter++});
        meshPtr->elementList.emplace_back(
              Tetrahedron{{
                            &meshPtr->pointList[base + dx + dy], // 2
                            &meshPtr->pointList[base + dx + dz], // 5
                            &meshPtr->pointList[base + dx + dy + dz], // 6
                            &meshPtr->pointList[base + dy + dz]}, // 7
             counter++});
        meshPtr->elementList.emplace_back(
              Tetrahedron{{
                            &meshPtr->pointList[base + dx + dy], // 2
                            &meshPtr->pointList[base + dz], // 4
                            &meshPtr->pointList[base + dx + dz], // 5
                            &meshPtr->pointList[base + dy + dz]}, // 7
             counter++});
        meshPtr->elementList.emplace_back(
              Tetrahedron{{
                            &meshPtr->pointList[base + dx + dy], // 2
                            &meshPtr->pointList[base + dy], // 3
                            &meshPtr->pointList[base + dz], // 4
                            &meshPtr->pointList[base + dy + dz]}, // 7
             counter++});
    }
  meshPtr->buildConnectivity();
  buildFacets(meshPtr);
  markFacetsCube(meshPtr, origin, length);
}
