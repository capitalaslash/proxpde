#include "mesh.hpp"

void buildMesh1D(Mesh<Line> & mesh,
                 Vec3 const& origin,
                 Vec3 const& length,
                 uint const numPts,
                 bool keepInternalFacets)
{
  assert(numPts > 1);
  assert(length[0] > 0.);
  Vec3 const h = length / (numPts-1);
  mesh.pointList.reserve(numPts);
  for(uint p=0; p<numPts; ++p)
  {
    mesh.pointList.emplace_back(origin + p * h, p);
  }

  uint const numElems = numPts-1;
  mesh.elementList.reserve(numElems);
  for(uint e=0; e<numElems; ++e)
  {
    Line line({&mesh.pointList[e], &mesh.pointList[e+1]}, e);
    mesh.elementList.push_back(std::move(line));
  }

  mesh.buildConnectivity();
  buildFacets(mesh, keepInternalFacets);
  mesh.facetList[0].marker = side::LEFT;
  mesh.facetList[1].marker = side::RIGHT;
}

template <typename Mesh>
void markFacetsCube(Mesh & mesh,
                    Vec3 const& origin,
                    Vec3 const& length)
{
  for(auto & f: mesh.facetList)
  {
    if (f.onBoundary())
    {
      auto const [min, max] = f.bbox();
      // BOTTOM -> ymin
      if(std::fabs(max[1]-origin[1]) < 1e-6*length[1])
      {
        f.marker = side::BOTTOM;
      }
      // RIGHT -> xmax
      else if(std::fabs(min[0]-origin[0]-length[0]) < 1e-6*length[0])
      {
        f.marker = side::RIGHT;
      }
      // TOP -> ymax
      else if(std::fabs(min[1]-origin[1]-length[1]) < 1e-6*length[1])
      {
        f.marker = side::TOP;
      }
      // LEFT -> xmin
      else if(std::fabs(max[0]-origin[0]) < 1e-6*length[0])
      {
        f.marker = side::LEFT;
      }
      // BACK -> zmin
      else if(std::fabs(max[2]-origin[2]) < 1e-6*length[2])
      {
        f.marker = side::BACK;
      }
      // FRONT -> zmax
      else if(std::fabs(min[2]-origin[2]-length[2]) < 1e-6*length[2])
      {
        f.marker = side::FRONT;
      }
    }
  }
}

void buildMesh2D(Mesh<Triangle> & mesh,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 2> const numPts,
                 bool keepInternalFacets)
{
  assert(numPts[0] > 1 && numPts[1] > 1);
  assert(length[0] > 0. && length[1] > 0.);
  Vec3 const h = {length(0) / (numPts[0]-1.), length(1) / (numPts[1]-1.), 1.};
  mesh.pointList.reserve(numPts[0]*numPts[1]);
  for(uint j=0; j<numPts[1]; ++j)
  {
    for(uint i=0; i<numPts[0]; ++i)
    {
      Vec3 p(i * h(0), j * h(1), 0.0);
      mesh.pointList.emplace_back(origin + p, i + numPts[0]*j);
    }
  }

  for(uint s=0; s<numPts[0]; s++)
  {
    mesh.pointList[s].marker = side::BOTTOM;
    mesh.pointList[(numPts[1]-1)*numPts[0]+s].marker = side::TOP;
  }
  for(uint s=0; s<numPts[1]; s++)
  {
    mesh.pointList[0+numPts[0]*s].marker = side::LEFT;
    mesh.pointList[numPts[0]-1+numPts[0]*s].marker = side::RIGHT;
  }

  uint const numElems = (numPts[0]-1)*(numPts[1]-1);
  mesh.elementList.reserve(numElems);
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
      mesh.elementList.emplace_back(
        Triangle{{&mesh.pointList[triplet_b[0]],
                  &mesh.pointList[triplet_b[1]],
                  &mesh.pointList[triplet_b[2]]},
                 counter++});
      mesh.elementList.emplace_back(
        Triangle{{&mesh.pointList[triplet_t[0]],
                  &mesh.pointList[triplet_t[1]],
                  &mesh.pointList[triplet_t[2]]},
                 counter++});
    }
  mesh.buildConnectivity();
  buildFacets(mesh, keepInternalFacets);
  markFacetsCube(mesh, origin, length);
}

void buildMesh2D(Mesh<Quad> & mesh,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 2> const numPts,
                 bool keepInternalFacets)
{
  assert(numPts[0] > 1 && numPts[1] > 1);
  assert(length[0] > 0. && length[1] > 0.);
  Vec3 const h = {length(0) / (numPts[0]-1.), length(1) / (numPts[1]-1.), 1.};
  mesh.pointList.reserve(numPts[0]*numPts[1]);
  for(uint j=0; j<numPts[1]; ++j)
  {
    for(uint i=0; i<numPts[0]; ++i)
    {
      Vec3 p(i * h(0), j * h(1), 0.0);
      mesh.pointList.emplace_back(origin + p, i + numPts[0]*j);
    }
  }

  for(uint s=0; s<numPts[0]; s++)
  {
    mesh.pointList[s].marker = side::BOTTOM;
    mesh.pointList[(numPts[1]-1)*numPts[0]+s].marker = side::TOP;
  }
  for(uint s=0; s<numPts[1]; s++)
  {
    mesh.pointList[0+numPts[0]*s].marker = side::LEFT;
    mesh.pointList[numPts[0]-1+numPts[0]*s].marker = side::RIGHT;
  }

  uint const numElems = (numPts[0]-1)*(numPts[1]-1);
  mesh.elementList.reserve(numElems);
  id_T counter = 0;
  for(uint j=0; j<numPts[1]-1; ++j)
    for(uint i=0; i<numPts[0]-1; ++i)
    {
      id_T const base = i + j*numPts[0];
      mesh.elementList.emplace_back(
        Quad{{&mesh.pointList[base],
              &mesh.pointList[base+1],
              &mesh.pointList[base+numPts[0]+1],
              &mesh.pointList[base+numPts[0]]},
             counter++});
    }
  mesh.buildConnectivity();
  buildFacets(mesh, keepInternalFacets);
  markFacetsCube(mesh, origin, length);
}

void buildCircleMesh(Mesh<Quad> & mesh,
                     Vec3 const& origin,
                     double const& radius,
                     array<uint, 3> const numPts)
{
  mesh.pointList.resize(
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
      mesh.pointList[id] = Point(origin + p, id);
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
      mesh.pointList[id] = Point(origin + p, id);
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
      mesh.pointList[id] = Point(origin + p, id);
    }
  }

  // for(uint s=0; s<numPts[0]; s++)
  // {
  //   mesh.pointList[s].marker = side::TOP;
  // }
  // for(uint s=0; s<numPts[1]; s++)
  // {
  //   mesh.pointList[0+numPts[0]*s].marker = side::LEFT;
  //   mesh.pointList[numPts[0]-1+numPts[0]*s].marker = side::RIGHT;
  // }

  uint const numElems =
    (numPts[0]-1)*(numPts[1]-1) +
    (numPts[0]-1)*(numPts[2]-1) +
    (numPts[1]-1)*(numPts[2]-1);

  mesh.elementList.reserve(numElems);
  id_T counter = 0;
  for(uint j=0; j<numPts[1]-1 + numPts[2]-1; ++j)
    for(uint i=0; i<numPts[0]-1; ++i)
    {
      id_T const base = i + j*numPts[0];
      mesh.elementList.emplace_back(
        Quad{{&mesh.pointList[base],
              &mesh.pointList[base+1],
              &mesh.pointList[base+numPts[0]+1],
              &mesh.pointList[base+numPts[0]]},
             counter++});
    }

  uint blockOffset = numPts[0]*numPts[1]+numPts[0]*(numPts[2]-1);
  for(uint j=0; j<numPts[1]-2; ++j)
  {
    id_T const base = numPts[0]-1 + j*numPts[0];
    mesh.elementList.emplace_back(
      Quad{{&mesh.pointList[base],
            &mesh.pointList[blockOffset+j*(numPts[2]-1)],
            &mesh.pointList[blockOffset+(j+1)*(numPts[2]-1)],
            &mesh.pointList[base+numPts[0]]},
            counter++});
  }
  mesh.elementList.emplace_back(
    Quad{{&mesh.pointList[numPts[0]-1 + (numPts[1]-2)*numPts[0]],
          &mesh.pointList[blockOffset+(numPts[1]-2)*(numPts[2]-1)],
          &mesh.pointList[numPts[0]-1 + (numPts[1]-2)*numPts[0]+numPts[0] + numPts[0]],
          &mesh.pointList[numPts[0]-1 + (numPts[1]-2)*numPts[0]+numPts[0]]},
          counter++});

  for(uint i=1; i<numPts[2]-1; ++i)
  {
    for(uint j=0; j<numPts[1]-2; ++j)
    {
      id_T const base = blockOffset + i-1 + j*(numPts[2]-1);
      mesh.elementList.emplace_back(
        Quad{{&mesh.pointList[base],
              &mesh.pointList[base+1],
              &mesh.pointList[base+numPts[2]-1+1],
              &mesh.pointList[base+numPts[2]-1]},
              counter++});
    }
    id_T const base = blockOffset + i-1 +(numPts[1]-2)*(numPts[2]-1);
    mesh.elementList.emplace_back(
      Quad{{&mesh.pointList[base],
            &mesh.pointList[base+1],
            &mesh.pointList[numPts[0]*numPts[1]-1+numPts[0]*(i+1)],
            &mesh.pointList[numPts[0]*numPts[1]-1+numPts[0]*i]},
            counter++});
  }

  mesh.buildConnectivity();
  buildFacets(mesh);
  for(auto & f: mesh.facetList)
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

void buildMesh3D(Mesh<Tetrahedron> & mesh,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 3> const numPts,
                 bool keepInternalFacets)
{
  assert(numPts[0] > 1 && numPts[1] > 1 && numPts[2] > 1);
  assert(length[0] > 0. && length[1] > 0. && length[2] > 0.);
  Vec3 const h = {length(0) / (numPts[0]-1.), length(1) / (numPts[1]-1.), length(2) / (numPts[2]-1.)};
  mesh.pointList.reserve(numPts[0]*numPts[1]*numPts[2]);
  for(uint k=0; k<numPts[2]; ++k)
  {
    for(uint j=0; j<numPts[1]; ++j)
    {
      for(uint i=0; i<numPts[0]; ++i)
      {
        Vec3 p(i * h(0), j * h(1), k * h(2));
        auto const id = i + numPts[0]*j + numPts[0]*numPts[1]*k;
        mesh.pointList.emplace_back(origin + p, id);
      }
    }
  }

  for(uint si=0; si<numPts[0]; si++)
    for(uint sj=0; sj<numPts[1]; sj++)
    {
      auto const sb = si + numPts[0]*sj;
      mesh.pointList[sb].marker = side::BACK;
      auto const sf = si + numPts[0]*sj + numPts[1]*numPts[0]*(numPts[2]-1);
      mesh.pointList[sf].marker = side::FRONT;
    }
  for(uint si=0; si<numPts[0]; si++)
    for(uint sk=0; sk<numPts[2]; sk++)
    {
      auto const sb = si + numPts[1]*numPts[0]*sk;
      mesh.pointList[sb].marker = side::BOTTOM;
      auto const st = si + numPts[0]*(numPts[1]-1) + numPts[1]*numPts[0]*sk;
      mesh.pointList[st].marker = side::TOP;
    }
  for(uint sj=0; sj<numPts[1]; sj++)
    for(uint sk=0; sk<numPts[2]; sk++)
    {
      auto const sl = numPts[0]*sj + numPts[1]*numPts[0]*sk;
      mesh.pointList[sl].marker = side::LEFT;
      auto const sr = numPts[0]-1 + numPts[0]*sj + numPts[1]*numPts[0]*sk;
      mesh.pointList[sr].marker = side::RIGHT;
    }

  uint const numElems = 5*(numPts[0]-1)*(numPts[1]-1)*(numPts[2]-1);
  mesh.elementList.reserve(numElems);
  id_T counter = 0;
  auto const dx = 1;
  auto const dy = numPts[0];
  auto const dz = numPts[1]*numPts[0];
  for(uint k=0; k<numPts[2]-1; ++k)
    for(uint j=0; j<numPts[1]-1; ++j)
      for(uint i=0; i<numPts[0]-1; ++i)
      {
        id_T const base = i + j*numPts[0] + k*numPts[1]*numPts[0];
        Eigen::Matrix<id_T, 5, 4> idTable;
        if ((i + j + k) % 2 )
        {
          idTable << base + dx, base + dy, base + dz, base + dx + dy + dz, // 1, 3, 4, 6
              base, base + dx, base + dy, base + dz, // 0, 1, 3, 4
              base + dx, base + dx + dy, base + dy, base + dx + dy + dz, // 1, 2, 3, 6
              base + dx, base + dz, base + dx + dz, base + dx + dy + dz, // 1, 4, 5, 6
              base + dy, base + dz, base + dx + dy + dz, base + dy + dz; // 3, 4, 6, 7
        }
        else
        {
          idTable << base, base + dx + dy, base + dx + dz, base + dy + dz, // 0, 2, 5, 7
              base, base + dx, base + dx + dy, base + dx + dz, // 0, 1, 2, 5
              base, base + dx + dy, base + dy, base + dy + dz, // 0, 2, 3, 7
              base, base + dz, base + dx + dz, base + dy + dz, // 0, 4, 5, 7
              base + dx + dy, base + dx + dz, base + dx + dy + dz, base + dy + dz; // 2, 5, 6, 7
        }
        for (uint r=0; r<5; ++r)
        {
          mesh.elementList.emplace_back(
                Tetrahedron{{
                              &mesh.pointList[idTable(r, 0)],
                              &mesh.pointList[idTable(r, 1)],
                              &mesh.pointList[idTable(r, 2)],
                              &mesh.pointList[idTable(r, 3)]},
                            counter++});
        }
      }
  mesh.buildConnectivity();
  buildFacets(mesh, keepInternalFacets);
  markFacetsCube(mesh, origin, length);
}

void buildMesh3D(Mesh<Hexahedron> & mesh,
                 Vec3 const& origin,
                 Vec3 const& length,
                 array<uint, 3> const numPts,
                 bool keepInternalFacets)
{
  assert(numPts[0] > 1 && numPts[1] > 1 && numPts[2] > 1);
  assert(length[0] > 0. && length[1] > 0. && length[2] > 0.);
  Vec3 const h = {length(0) / (numPts[0]-1.), length(1) / (numPts[1]-1.), length(2) / (numPts[2]-1.)};
  mesh.pointList.reserve(numPts[0]*numPts[1]*numPts[2]);
  for(uint k=0; k<numPts[2]; ++k)
  {
    for(uint j=0; j<numPts[1]; ++j)
    {
      for(uint i=0; i<numPts[0]; ++i)
      {
        Vec3 p(i * h(0), j * h(1), k * h(2));
        auto const id = i + numPts[0]*j + numPts[0]*numPts[1]*k;
        mesh.pointList.emplace_back(origin + p, id);
      }
    }
  }

  for(uint si=0; si<numPts[0]; si++)
    for(uint sj=0; sj<numPts[1]; sj++)
    {
      auto const sb = si + numPts[0]*sj;
      mesh.pointList[sb].marker = side::BACK;
      auto const sf = si + numPts[0]*sj + numPts[1]*numPts[0]*(numPts[2]-1);
      mesh.pointList[sf].marker = side::FRONT;
    }
  for(uint si=0; si<numPts[0]; si++)
    for(uint sk=0; sk<numPts[2]; sk++)
    {
      auto const sb = si + numPts[1]*numPts[0]*sk;
      mesh.pointList[sb].marker = side::BOTTOM;
      auto const st = si + numPts[0]*(numPts[1]-1) + numPts[1]*numPts[0]*sk;
      mesh.pointList[st].marker = side::TOP;
    }
  for(uint sj=0; sj<numPts[1]; sj++)
    for(uint sk=0; sk<numPts[2]; sk++)
    {
      auto const sl = numPts[0]*sj + numPts[1]*numPts[0]*sk;
      mesh.pointList[sl].marker = side::LEFT;
      auto const sr = numPts[0]-1 + numPts[0]*sj + numPts[1]*numPts[0]*sk;
      mesh.pointList[sr].marker = side::RIGHT;
    }

  uint const numElems = (numPts[0]-1)*(numPts[1]-1)*(numPts[2]-1);
  mesh.elementList.reserve(numElems);
  id_T counter = 0;
  auto const dx = 1;
  auto const dy = numPts[0];
  auto const dz = numPts[1]*numPts[0];
  for(uint k=0; k<numPts[2]-1; ++k)
    for(uint j=0; j<numPts[1]-1; ++j)
      for(uint i=0; i<numPts[0]-1; ++i)
      {
        id_T const base = i + j*numPts[0] + k*numPts[1]*numPts[0];
        mesh.elementList.emplace_back(
              Hexahedron{{
                           &mesh.pointList[base], // 0
                           &mesh.pointList[base + dx], // 1
                           &mesh.pointList[base + dx + dy], // 2
                           &mesh.pointList[base + dy], // 3
                           &mesh.pointList[base + dz], // 4
                           &mesh.pointList[base + dx + dz], // 5
                           &mesh.pointList[base + dx + dy + dz], // 6
                           &mesh.pointList[base + dy + dz]}, // 7
                         counter++});
      }
  mesh.buildConnectivity();
  buildFacets(mesh, keepInternalFacets);
  markFacetsCube(mesh, origin, length);
}

void refTriangleMesh(Mesh<Triangle> & mesh)
{
  mesh.pointList = {
    Point(Vec3(0., 0., 0.), 0),
    Point(Vec3(1., 0., 0.), 1),
    Point(Vec3(0., 1., 0.), 2),
  };
  mesh.elementList = {
    Triangle{{&mesh.pointList[0],
              &mesh.pointList[1],
              &mesh.pointList[2]},
              0},
  };
  mesh.buildConnectivity();
}

void hexagonMesh(Mesh<Triangle> & mesh)
{
  double const sr3o2 = 0.5*sqrt(3.);
  mesh.pointList = {
    Point(Vec3(0., 0., 0.), 0),
    Point(Vec3(    0.,  1.0, 0.), 1),
    Point(Vec3( sr3o2,  0.5, 0.), 2),
    Point(Vec3( sr3o2, -0.5, 0.), 3),
    Point(Vec3(    0., -1.0, 0.), 4),
    Point(Vec3(-sr3o2, -0.5, 0.), 5),
    Point(Vec3(-sr3o2,  0.5, 0.), 6),
  };
  mesh.elementList = {
    Triangle{{&mesh.pointList[0],
              &mesh.pointList[1],
              &mesh.pointList[2]},
              0},
    Triangle{{&mesh.pointList[0],
              &mesh.pointList[2],
              &mesh.pointList[3]},
              1},
    Triangle{{&mesh.pointList[0],
              &mesh.pointList[3],
              &mesh.pointList[4]},
              2},
    Triangle{{&mesh.pointList[0],
              &mesh.pointList[4],
              &mesh.pointList[5]},
              3},
    Triangle{{&mesh.pointList[0],
              &mesh.pointList[5],
              &mesh.pointList[6]},
              4},
    Triangle{{&mesh.pointList[0],
              &mesh.pointList[6],
              &mesh.pointList[1]},
              5},
  };
  mesh.buildConnectivity();
  buildFacets(mesh, true);
  for (auto & facet: mesh.facetList)
  {
    if (facet.onBoundary())
    {
      facet.marker = 1;
      for (auto point: facet.pointList)
      {
        point->marker = 1;
      }
    }
  }
  buildNormals(mesh);
  std::cout << mesh << std::endl;
}

void hexagonSquare(Mesh<Triangle> & mesh, bool keepInternalFacets)
{
  id_T counter = 0;
  double const h = 0.1;
  for (uint p=0; p<11; ++p) // id 1-10
  {
    mesh.pointList.emplace_back(Vec3(p*h, 0.0, 0.0), counter++);
  }
  double const oneOnSr3 = 1. / std::sqrt(3.);
  uint const rows = 8;
  for (uint r=0; r<rows; ++r)
  {
    for (uint p=0; p<5; ++p) // id 11-15
    {
      mesh.pointList.emplace_back(Vec3((1+2*p)*h, (1+2*r)*h*oneOnSr3, 0.0), counter++);
    }
    for (uint p=0; p<6; ++p) // id 16-21
    {
      mesh.pointList.emplace_back(Vec3(2*p*h, 2*(1+r)*h*oneOnSr3, 0.0), counter++);
    }
  }
  for (uint p=0; p<5; ++p) // id 22-26
  {
    mesh.pointList.emplace_back(Vec3((1+2*p)*h, (1+2*rows)*h*oneOnSr3, 0.0), counter++);
  }
  for (uint p=0; p<11; ++p) // id 27-37
  {
    mesh.pointList.emplace_back(Vec3(p*h, (2+2*rows)*h*oneOnSr3, 0.0), counter++);
  }

  counter = 0;
  for (uint e=0; e<5; ++e)
  {
    mesh.elementList.emplace_back(std::vector<Point*>{
      &mesh.pointList[0 + 2*e],
      &mesh.pointList[1 + 2*e],
      &mesh.pointList[11 + e]},
      counter++);
    mesh.elementList.emplace_back(std::vector<Point*>{
      &mesh.pointList[1 + 2*e],
      &mesh.pointList[2 + 2*e],
      &mesh.pointList[11 + e]},
      counter++);
  }
  for (uint e=0; e<5; ++e)
  {
    mesh.elementList.emplace_back(std::vector<Point*>{
                                    &mesh.pointList[0 + 2*e],
                                    &mesh.pointList[11 + e],
                                    &mesh.pointList[16 + e]},
                                  counter++);
    mesh.elementList.emplace_back(std::vector<Point*>{
                                    &mesh.pointList[2 + 2*e],
                                    &mesh.pointList[11 + e],
                                    &mesh.pointList[17 + e]},
                                  counter++);
  }
  for (uint r=0; r<rows-1; ++r)
  {
    for (uint e=0; e<5; ++e)
    {
      mesh.elementList.emplace_back(std::vector<Point*>{
                                      &mesh.pointList[16 + 11*r + e],
                                      &mesh.pointList[11 + 11*r + e],
                                      &mesh.pointList[22 + 11*r + e]},
                                    counter++);
      mesh.elementList.emplace_back(std::vector<Point*>{
                                      &mesh.pointList[17 + 11*r + e],
                                      &mesh.pointList[11 + 11*r + e],
                                      &mesh.pointList[22 + 11*r + e]},
                                    counter++);
    }
    for (uint e=0; e<5; ++e)
    {
      mesh.elementList.emplace_back(std::vector<Point*>{
                                      &mesh.pointList[16 + 11*r + e],
                                      &mesh.pointList[22 + 11*r + e],
                                      &mesh.pointList[27 + 11*r + e]},
                                    counter++);
      mesh.elementList.emplace_back(std::vector<Point*>{
                                      &mesh.pointList[17 + 11*r + e],
                                      &mesh.pointList[22 + 11*r + e],
                                      &mesh.pointList[28 + 11*r + e]},
                                    counter++);
    }
  }
  for (uint e=0; e<5; ++e)
  {
    mesh.elementList.emplace_back(std::vector<Point*>{
      &mesh.pointList[11 + 11*(rows-1) + e],
      &mesh.pointList[16 + 11*(rows-1) + e],
      &mesh.pointList[22 + 11*(rows-1) + e]},
      counter++);
    mesh.elementList.emplace_back(std::vector<Point*>{
      &mesh.pointList[11 + 11*(rows-1) + e],
      &mesh.pointList[17 + 11*(rows-1) + e],
      &mesh.pointList[22 + 11*(rows-1) + e]},
      counter++);
  }
  for (uint e=0; e<5; ++e)
  {
    mesh.elementList.emplace_back(std::vector<Point*>{
      &mesh.pointList[16 + 11*(rows-1) + e],
      &mesh.pointList[22 + 11*(rows-1) + e],
      &mesh.pointList[27 + 11*(rows-1) + 2*e]},
      counter++);
    mesh.elementList.emplace_back(std::vector<Point*>{
      &mesh.pointList[17 + 11*(rows-1) + e],
      &mesh.pointList[22 + 11*(rows-1) + e],
      &mesh.pointList[29 + 11*(rows-1) + 2*e]},
      counter++);
  }
  for (uint e=0; e<5; ++e)
  {
    mesh.elementList.emplace_back(std::vector<Point*>{
      &mesh.pointList[27 + 11*(rows-1) + 2*e],
      &mesh.pointList[28 + 11*(rows-1) + 2*e],
      &mesh.pointList[22 + 11*(rows-1) + e]},
      counter++);
    mesh.elementList.emplace_back(std::vector<Point*>{
      &mesh.pointList[28 + 11*(rows-1) + 2*e],
      &mesh.pointList[29 + 11*(rows-1) + 2*e],
      &mesh.pointList[22 + 11*(rows-1) + e]},
      counter++);
  }
  mesh.buildConnectivity();
  buildFacets(mesh, keepInternalFacets);
  markFacetsCube(mesh, {0., 0., 0.}, {1., (2+2*rows)*h*oneOnSr3, 0.});
  buildNormals(mesh);
}

void refTetrahedronMesh(Mesh<Tetrahedron> & mesh)
{
  mesh.pointList = {
    Point(Vec3(0., 0., 0.), 0),
    Point(Vec3(1., 0., 0.), 1),
    Point(Vec3(0., 1., 0.), 2),
    Point(Vec3(0., 0., 1.), 3),
  };
  mesh.elementList = {
    Tetrahedron {
      {&mesh.pointList[0],
       &mesh.pointList[1],
       &mesh.pointList[2],
       &mesh.pointList[3]},
      0},
  };
  mesh.buildConnectivity();
}
