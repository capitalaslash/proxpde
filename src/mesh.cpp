#include "mesh.hpp"

#include "geo.hpp"

namespace proxpde
{

// explicit instantiations
#ifdef PROXPDE_EXPLICIT_INSTANTIATION
template class Mesh<Line>;
template class Mesh<Triangle>;
template class Mesh<Quad>;
template class Mesh<Tetrahedron>;
template class Mesh<Hexahedron>;
#endif

void buildLine(
    Mesh<Line> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    uint const numElems,
    Bitmask<MeshFlags> flags)
{
  assert(numElems > 0);
  assert(length.norm() > 0.);
  Vec3 const h = length / numElems;
  mesh.pointList.reserve(numElems + 1);
  for (uint p = 0; p < numElems + 1; ++p)
  {
    mesh.pointList.emplace_back(origin + p * h, p);
  }

  mesh.elementList.reserve(numElems);
  for (uint e = 0; e < numElems; ++e)
  {
    Line line({&mesh.pointList[e], &mesh.pointList[e + 1]}, e);
    mesh.elementList.push_back(std::move(line));
  }

  mesh.buildConnectivity();
  buildFacets(mesh, flags);
  mesh.facetList[0].marker = side::LEFT;
  mesh.facetList[1].marker = side::RIGHT;
  mesh.physicalNames["left"] = side::LEFT;
  mesh.physicalNames["right"] = side::RIGHT;
}

void buildSquare(
    Mesh<Triangle> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 2> const numElems,
    Bitmask<MeshFlags> flags)
{
  assert(numElems[0] > 0 && numElems[1] > 0);
  // TODO: generalize for planes in 3D space
  assert(length[0] > 0. && length[1] > 0.);
  Vec3 const h = {length[0] / numElems[0], length[1] / numElems[1], 1.};
  std::array<uint, 2> const numPts = {numElems[0] + 1, numElems[1] + 1};
  mesh.pointList.reserve(numPts[0] * numPts[1]);
  for (uint j = 0; j < numPts[1]; ++j)
  {
    for (uint i = 0; i < numPts[0]; ++i)
    {
      Vec3 p(i * h(0), j * h(1), 0.0);
      mesh.pointList.emplace_back(origin + p, i + numPts[0] * j);
    }
  }

  for (uint s = 0; s < numPts[0]; s++)
  {
    mesh.pointList[s].marker = side::BOTTOM;
    mesh.pointList[numElems[1] * numPts[0] + s].marker = side::TOP;
  }
  for (uint s = 0; s < numPts[1]; s++)
  {
    mesh.pointList[0 + numPts[0] * s].marker = side::LEFT;
    mesh.pointList[numElems[0] + numPts[0] * s].marker = side::RIGHT;
  }

  mesh.elementList.reserve(numElems[0] * numElems[1]);
  id_T counter = 0;
  for (uint j = 0; j < numElems[1]; ++j)
  {
    bool flip = (j % 2 == 1);
    for (uint i = 0; i < numElems[0]; ++i)
    {
      id_T const base = i + j * numPts[0];
      std::array<std::array<id_T, 3>, 2> triplets;
      if (flip)
      {
        triplets[0][0] = base;
        triplets[0][1] = base + 1;
        triplets[0][2] = base + numPts[0];
        triplets[1][0] = base + 1;
        triplets[1][1] = base + numPts[0] + 1;
        triplets[1][2] = base + numPts[0];
        flip = false;
      }
      else
      {
        triplets[0][0] = base;
        triplets[0][1] = base + 1;
        triplets[0][2] = base + numPts[0] + 1;
        triplets[1][0] = base;
        triplets[1][1] = base + numPts[0] + 1;
        triplets[1][2] = base + numPts[0];
        flip = true;
      }
      for (uint e = 0; e < 2U; ++e)
        mesh.elementList.emplace_back(
            std::vector{
                {&mesh.pointList[triplets[e][0]],
                 &mesh.pointList[triplets[e][1]],
                 &mesh.pointList[triplets[e][2]]}},
            counter++);
    }
  }
  mesh.buildConnectivity();
  buildFacets(mesh, flags);
  markFacetsCube(mesh, origin, length);
}

void buildSquare(
    Mesh<Quad> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 2> const numElems,
    Bitmask<MeshFlags> flags)
{
  assert(numElems[0] > 0 && numElems[1] > 0);
  // TODO: generalize for planes in 3D space
  assert(length[0] > 0. && length[1] > 0.);
  Vec3 const h = {length[0] / numElems[0], length[1] / numElems[1], 1.};
  std::array<uint, 2> const numPts = {numElems[0] + 1, numElems[1] + 1};

  mesh.pointList.reserve(numPts[0] * numPts[1]);
  for (uint j = 0; j < numPts[1]; ++j)
  {
    for (uint i = 0; i < numPts[0]; ++i)
    {
      Vec3 p(i * h(0), j * h(1), 0.0);
      mesh.pointList.emplace_back(origin + p, i + numPts[0] * j);
    }
  }

  for (uint s = 0; s < numPts[0]; s++)
  {
    mesh.pointList[s].marker = side::BOTTOM;
    mesh.pointList[numElems[1] * numPts[0] + s].marker = side::TOP;
  }
  for (uint s = 0; s < numPts[1]; s++)
  {
    mesh.pointList[0 + numPts[0] * s].marker = side::LEFT;
    mesh.pointList[numElems[0] + numPts[0] * s].marker = side::RIGHT;
  }

  mesh.elementList.reserve(numElems[0] * numElems[1]);
  id_T counter = 0;
  for (uint j = 0; j < numElems[1]; ++j)
    for (uint i = 0; i < numElems[0]; ++i)
    {
      id_T const base = i + j * numPts[0];
      mesh.elementList.emplace_back(
          std::vector{
              {&mesh.pointList[base],
               &mesh.pointList[base + 1],
               &mesh.pointList[base + numPts[0] + 1],
               &mesh.pointList[base + numPts[0]]}},
          counter++);
    }
  mesh.buildConnectivity();
  buildFacets(mesh, flags);
  markFacetsCube(mesh, origin, length);
}

void buildCircleMesh(
    Mesh<Quad> & mesh,
    Vec3 const & origin,
    double const & radius,
    std::array<uint, 3> const numElems)
{
  std::array<uint, 3> const numPts = {
      numElems[0] + 1, numElems[1] + 1, numElems[2] + 1};
  mesh.pointList.resize(
      numPts[0] * numPts[1] + numPts[0] * numElems[2] + numElems[1] * numElems[2]);
  //    l0  r0
  //    + - + - +
  //    | A | C |
  // la + - +   |
  //    |   ra
  //      B   +
  //    + - -
  Vec3 r0(0.6 * radius, 0.0, 0.0);
  Vec3 ra(sqrt(2.) * 0.3 * radius, -sqrt(2.) * 0.3 * radius, 0.0);
  Vec3 l0(0.0, 0.0, 0.0);
  Vec3 la(0.0, -0.6 * radius, 0.0);
  // section A
  for (uint j = 0; j < numPts[1]; ++j)
  {
    double const sj = static_cast<double>(j) / numElems[1];
    // point between l0 and la
    Vec3 b = (1. - sj) * l0 + sj * la;
    // point between r0 and ra
    Vec3 e = (1. - sj) * r0 + sj * ra;
    for (uint i = 0; i < numPts[0]; ++i)
    {
      double si = static_cast<double>(i) / numElems[0];
      // point between b and e
      Vec3 p = (1. - si) * b + si * e;
      const id_T id = i + numPts[0] * j;
      mesh.pointList[id] = Point(origin + p, id);
    }
  }

  uint offset = numPts[0] * numPts[1];
  // section B
  for (uint i = 0; i < numPts[0]; ++i)
  {
    double const alpha = i * 0.25 * M_PI / numElems[0];
    double const si = static_cast<double>(i) / numElems[0];
    // point between la and ra
    Vec3 b = la * (1. - si) + ra * si;
    // point on circle
    Vec3 e(std::sin(alpha) * radius, -std::cos(alpha) * radius, 0.0);
    for (uint r = 1; r < numPts[2]; ++r)
    {
      double const sr = static_cast<double>(r) / numElems[2];
      // point between b and e
      Vec3 const p = b * (1 - sr) + e * sr;
      const id_T id = offset + i + numPts[0] * (r - 1);
      mesh.pointList[id] = Point(origin + p, id);
    }
  }

  offset += numPts[0] * (numPts[2] - 1);
  // section C
  for (uint j = 0; j < numPts[1] - 1; ++j)
  {
    double const alpha = j * 0.25 * M_PI / numElems[1];
    double const sj = static_cast<double>(j) / numElems[1];
    // point between r0 and ra
    Vec3 b = r0 * (1. - sj) + ra * sj;
    // point on circle
    Vec3 e(std::cos(alpha), -std::sin(alpha), 0.0);
    for (uint r = 1; r < numPts[2]; ++r)
    {
      double const sr = static_cast<double>(r) / numElems[2];
      // point between b  and e
      Vec3 p = b * (1. - sr) + e * sr;
      const id_T id = offset + r - 1 + (numPts[2] - 1) * j;
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

  uint const totalNumElems =
      numElems[0] * numElems[1] + numElems[0] * numElems[2] + numElems[1] * numElems[2];

  mesh.elementList.reserve(totalNumElems);
  id_T counter = 0;
  for (uint j = 0; j < numElems[1] + numElems[2]; ++j)
    for (uint i = 0; i < numElems[0]; ++i)
    {
      id_T const base = i + j * numPts[0];
      mesh.elementList.emplace_back(
          std::vector{
              {&mesh.pointList[base],
               &mesh.pointList[base + 1],
               &mesh.pointList[base + numPts[0] + 1],
               &mesh.pointList[base + numPts[0]]}},
          counter++);
    }

  uint blockOffset = numPts[0] * numPts[1] + numPts[0] * (numPts[2] - 1);
  for (uint j = 0; j < numElems[1] - 1; ++j)
  {
    id_T const base = numPts[0] - 1 + j * numPts[0];
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[base],
             &mesh.pointList[blockOffset + j * (numPts[2] - 1)],
             &mesh.pointList[blockOffset + (j + 1) * (numPts[2] - 1)],
             &mesh.pointList[base + numPts[0]]}},
        counter++);
  }
  mesh.elementList.emplace_back(
      std::vector{
          {&mesh.pointList[numPts[0] - 1 + (numPts[1] - 2) * numPts[0]],
           &mesh.pointList[blockOffset + (numPts[1] - 2) * (numPts[2] - 1)],
           &mesh.pointList
                [numPts[0] - 1 + (numPts[1] - 2) * numPts[0] + numPts[0] + numPts[0]],
           &mesh.pointList[numPts[0] - 1 + (numPts[1] - 2) * numPts[0] + numPts[0]]}},
      counter++);

  for (uint i = 1; i < numElems[2]; ++i)
  {
    for (uint j = 0; j < numElems[1] - 1; ++j)
    {
      id_T const base = blockOffset + i - 1 + j * (numPts[2] - 1);
      mesh.elementList.emplace_back(
          std::vector{
              {&mesh.pointList[base],
               &mesh.pointList[base + 1],
               &mesh.pointList[base + numPts[2] - 1 + 1],
               &mesh.pointList[base + numPts[2] - 1]}},
          counter++);
    }
    id_T const base = blockOffset + i - 1 + (numPts[1] - 2) * (numPts[2] - 1);
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[base],
             &mesh.pointList[base + 1],
             &mesh.pointList[numPts[0] * numPts[1] - 1 + numPts[0] * (i + 1)],
             &mesh.pointList[numPts[0] * numPts[1] - 1 + numPts[0] * i]}},
        counter++);
  }

  mesh.buildConnectivity();
  buildFacets(mesh);
  for (auto & f: mesh.facetList)
  {
    if (std::fabs(f.pts[0]->coord[1] - origin[1]) < 1e-6 * radius &&
        std::fabs(f.pts[1]->coord[1] - origin[1]) < 1e-6 * radius)
    {
      f.marker = side::TOP;
    }
    else if (
        std::fabs(f.pts[0]->coord[0] - origin[0]) < 1e-6 * radius &&
        std::fabs(f.pts[1]->coord[0] - origin[0]) < 1e-6 * radius)
    {
      f.marker = side::LEFT;
    }
    else if (
        std::fabs((f.pts[0]->coord - origin).norm() - radius) < 1e-6 * radius &&
        std::fabs((f.pts[1]->coord - origin).norm() - radius) < 1e-6 * radius)
    {
      f.marker = side::CIRCLE;
    }
  }
  mesh.physicalNames["left"] = side::LEFT;
  mesh.physicalNames["top"] = side::TOP;
  mesh.physicalNames["circle"] = side::CIRCLE;
}

void buildCube(
    Mesh<Tetrahedron> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 3> const numElems,
    Bitmask<MeshFlags> flags)
{
  assert(numElems[0] > 0 && numElems[1] > 0 && numElems[2] > 0);
  assert(length[0] > 0. && length[1] > 0. && length[2] > 0.);
  Vec3 const h = {
      length[0] / numElems[0], length[1] / numElems[1], length[2] / numElems[2]};
  std::array<uint, 3> const numPts = {
      numElems[0] + 1, numElems[1] + 1, numElems[2] + 1};
  mesh.pointList.reserve(numPts[0] * numPts[1] * numPts[2]);
  for (uint k = 0; k < numPts[2]; ++k)
  {
    for (uint j = 0; j < numPts[1]; ++j)
    {
      for (uint i = 0; i < numPts[0]; ++i)
      {
        Vec3 p(i * h(0), j * h(1), k * h(2));
        auto const id = i + numPts[0] * j + numPts[0] * numPts[1] * k;
        mesh.pointList.emplace_back(origin + p, id);
      }
    }
  }

  for (uint si = 0; si < numPts[0]; si++)
    for (uint sj = 0; sj < numPts[1]; sj++)
    {
      auto const sb = si + numPts[0] * sj;
      mesh.pointList[sb].marker = side::BACK;
      auto const sf = si + numPts[0] * sj + numPts[1] * numPts[0] * (numPts[2] - 1);
      mesh.pointList[sf].marker = side::FRONT;
    }
  for (uint si = 0; si < numPts[0]; si++)
    for (uint sk = 0; sk < numPts[2]; sk++)
    {
      auto const sb = si + numPts[1] * numPts[0] * sk;
      mesh.pointList[sb].marker = side::BOTTOM;
      auto const st = si + numPts[0] * (numPts[1] - 1) + numPts[1] * numPts[0] * sk;
      mesh.pointList[st].marker = side::TOP;
    }
  for (uint sj = 0; sj < numPts[1]; sj++)
    for (uint sk = 0; sk < numPts[2]; sk++)
    {
      auto const sl = numPts[0] * sj + numPts[1] * numPts[0] * sk;
      mesh.pointList[sl].marker = side::LEFT;
      auto const sr = numPts[0] - 1 + numPts[0] * sj + numPts[1] * numPts[0] * sk;
      mesh.pointList[sr].marker = side::RIGHT;
    }

  mesh.elementList.reserve(5U * numElems[0] * numElems[1] * numElems[2]);
  id_T counter = 0;
  auto const dx = 1;
  auto const dy = numPts[0];
  auto const dz = numPts[1] * numPts[0];
  for (uint k = 0; k < numElems[2]; ++k)
    for (uint j = 0; j < numElems[1]; ++j)
      for (uint i = 0; i < numElems[0]; ++i)
      {
        id_T const base = i + j * numPts[0] + k * numPts[1] * numPts[0];
        Eigen::Matrix<id_T, 5, 4> idTable;
        if ((i + j + k) % 2)
        {
          // clang-format off
          idTable <<
              base + dx, base + dy,      base + dz,           base + dx + dy + dz, // 1, 3, 4, 6
              base,      base + dx,      base + dy,           base + dz,           // 0, 1, 3, 4
              base + dx, base + dx + dy, base + dy,           base + dx + dy + dz, // 1, 2, 3, 6
              base + dx, base + dz,      base + dx + dz,      base + dx + dy + dz, // 1, 4, 5, 6
              base + dy, base + dz,      base + dx + dy + dz, base + dy + dz;      // 3, 4, 6, 7
          // clang-format on
        }
        else
        {
          // clang-format off
          idTable <<
              base,           base + dx + dy, base + dx + dz,      base + dy + dz, // 0, 2, 5, 7
              base,           base + dx,      base + dx + dy,      base + dx + dz, // 0, 1, 2, 5
              base,           base + dx + dy, base + dy,           base + dy + dz, // 0, 2, 3, 7
              base,           base + dz,      base + dx + dz,      base + dy + dz, // 0, 4, 5, 7
              base + dx + dy, base + dx + dz, base + dx + dy + dz, base + dy + dz; // 2, 5, 6, 7
          // clang-format on
        }
        for (uint r = 0; r < 5; ++r)
        {
          mesh.elementList.emplace_back(
              std::vector{
                  {&mesh.pointList[idTable(r, 0)],
                   &mesh.pointList[idTable(r, 1)],
                   &mesh.pointList[idTable(r, 2)],
                   &mesh.pointList[idTable(r, 3)]}},
              counter++);
        }
      }
  mesh.buildConnectivity();
  buildFacets(mesh, flags);
  markFacetsCube(mesh, origin, length);
}

void buildCube(
    Mesh<Hexahedron> & mesh,
    Vec3 const & origin,
    Vec3 const & length,
    std::array<uint, 3> const numElems,
    Bitmask<MeshFlags> flags)
{
  assert(numElems[0] > 0 && numElems[1] > 0 && numElems[2] > 0);
  assert(length[0] > 0. && length[1] > 0. && length[2] > 0.);
  Vec3 const h = {
      length[0] / numElems[0], length[1] / numElems[1], length[2] / numElems[2]};
  std::array<uint, 3> const numPts = {
      numElems[0] + 1, numElems[1] + 1, numElems[2] + 1};
  mesh.pointList.reserve(numPts[0] * numPts[1] * numPts[2]);
  for (uint k = 0; k < numPts[2]; ++k)
  {
    for (uint j = 0; j < numPts[1]; ++j)
    {
      for (uint i = 0; i < numPts[0]; ++i)
      {
        Vec3 p(i * h(0), j * h(1), k * h(2));
        auto const id = i + numPts[0] * j + numPts[0] * numPts[1] * k;
        mesh.pointList.emplace_back(origin + p, id);
      }
    }
  }

  for (uint si = 0; si < numPts[0]; si++)
    for (uint sj = 0; sj < numPts[1]; sj++)
    {
      auto const sb = si + numPts[0] * sj;
      mesh.pointList[sb].marker = side::BACK;
      auto const sf = si + numPts[0] * sj + numPts[1] * numPts[0] * (numPts[2] - 1);
      mesh.pointList[sf].marker = side::FRONT;
    }
  for (uint si = 0; si < numPts[0]; si++)
    for (uint sk = 0; sk < numPts[2]; sk++)
    {
      auto const sb = si + numPts[1] * numPts[0] * sk;
      mesh.pointList[sb].marker = side::BOTTOM;
      auto const st = si + numPts[0] * (numPts[1] - 1) + numPts[1] * numPts[0] * sk;
      mesh.pointList[st].marker = side::TOP;
    }
  for (uint sj = 0; sj < numPts[1]; sj++)
    for (uint sk = 0; sk < numPts[2]; sk++)
    {
      auto const sl = numPts[0] * sj + numPts[1] * numPts[0] * sk;
      mesh.pointList[sl].marker = side::LEFT;
      auto const sr = numPts[0] - 1 + numPts[0] * sj + numPts[1] * numPts[0] * sk;
      mesh.pointList[sr].marker = side::RIGHT;
    }

  mesh.elementList.reserve(numElems[0] * numElems[1] * numElems[2]);
  id_T counter = 0;
  auto const dx = 1;
  auto const dy = numPts[0];
  auto const dz = numPts[1] * numPts[0];
  for (uint k = 0; k < numElems[2]; ++k)
    for (uint j = 0; j < numElems[1]; ++j)
      for (uint i = 0; i < numElems[0]; ++i)
      {
        id_T const base = i + j * numPts[0] + k * numPts[1] * numPts[0];
        mesh.elementList.emplace_back(
            std::vector{
                {&mesh.pointList[base],                // 0
                 &mesh.pointList[base + dx],           // 1
                 &mesh.pointList[base + dx + dy],      // 2
                 &mesh.pointList[base + dy],           // 3
                 &mesh.pointList[base + dz],           // 4
                 &mesh.pointList[base + dx + dz],      // 5
                 &mesh.pointList[base + dx + dy + dz], // 6
                 &mesh.pointList[base + dy + dz]}},    // 7
            counter++);
      }

  mesh.buildConnectivity();
  buildFacets(mesh, flags);
  markFacetsCube(mesh, origin, length);
}

void buildWedge(
    Mesh<Triangle> & mesh,
    Vec3 const & origin,
    Vec3 const & radius,
    Vec3 const & normal,
    uint const numLayers,
    double const angle)
{
  mesh.pointList.reserve(numLayers * 3U + 1U);
  auto ptCounter = 0U;
  mesh.pointList.emplace_back(origin, ptCounter++);
  double const h = 1.0 / numLayers;
  RotationMatrix rotMinus{normal, -0.5 * angle};
  RotationMatrix rotPlus{normal, +0.5 * angle};
  for (uint k = 0U; k < numLayers; k++)
  {
    Vec3 midPoint = origin + (k + 1) * h * radius;
    mesh.pointList.emplace_back(rotMinus.apply(midPoint), ptCounter, side::RIGHT);
    mesh.pointList.emplace_back(midPoint, ptCounter + 1U, markerNotSet);
    mesh.pointList.emplace_back(rotPlus.apply(midPoint), ptCounter + 2U, side::LEFT);
    ptCounter += 3U;
  }

  mesh.elementList.reserve(4U * numLayers - 2U);
  auto elemCounter = 0U;
  mesh.elementList.emplace_back(
      Triangle::Pts_T{&mesh.pointList[0], &mesh.pointList[1], &mesh.pointList[2]},
      elemCounter++);
  mesh.elementList.emplace_back(
      Triangle::Pts_T{&mesh.pointList[0], &mesh.pointList[2], &mesh.pointList[3]},
      elemCounter++);
  for (uint k = 0U; k < numLayers - 1; k++)
  {
    if (k % 2)
    {
      mesh.elementList.emplace_back(
          Triangle::Pts_T{
              &mesh.pointList[3 * k + 1],
              &mesh.pointList[3 * k + 4],
              &mesh.pointList[3 * k + 2],
          },
          elemCounter);
      mesh.elementList.emplace_back(
          Triangle::Pts_T{
              &mesh.pointList[3 * k + 2],
              &mesh.pointList[3 * k + 4],
              &mesh.pointList[3 * k + 5],
          },
          elemCounter + 1U);
      mesh.elementList.emplace_back(
          Triangle::Pts_T{
              &mesh.pointList[3 * k + 2],
              &mesh.pointList[3 * k + 5],
              &mesh.pointList[3 * k + 6],
          },
          elemCounter + 2U);
      mesh.elementList.emplace_back(
          Triangle::Pts_T{
              &mesh.pointList[3 * k + 2],
              &mesh.pointList[3 * k + 6],
              &mesh.pointList[3 * k + 3],
          },
          elemCounter + 3U);
    }
    else
    {
      mesh.elementList.emplace_back(
          Triangle::Pts_T{
              &mesh.pointList[3 * k + 1],
              &mesh.pointList[3 * k + 4],
              &mesh.pointList[3 * k + 5],
          },
          elemCounter);
      mesh.elementList.emplace_back(
          Triangle::Pts_T{
              &mesh.pointList[3 * k + 1],
              &mesh.pointList[3 * k + 5],
              &mesh.pointList[3 * k + 2],
          },
          elemCounter + 1U);
      mesh.elementList.emplace_back(
          Triangle::Pts_T{
              &mesh.pointList[3 * k + 2],
              &mesh.pointList[3 * k + 5],
              &mesh.pointList[3 * k + 3],
          },
          elemCounter + 2U);
      mesh.elementList.emplace_back(
          Triangle::Pts_T{
              &mesh.pointList[3 * k + 3],
              &mesh.pointList[3 * k + 5],
              &mesh.pointList[3 * k + 6],
          },
          elemCounter + 3U);
    }
    elemCounter += 4U;
  }

  mesh.buildConnectivity();
  buildFacets(mesh, MeshFlags::NONE);
  for (auto & facet: mesh.facetList)
  {
    if (facet.pts[0]->marker == facet.pts[1]->marker)
      facet.marker = facet.pts[0]->marker;
    else if (facet.pts[0]->id == 0U)
      facet.marker = facet.pts[1]->marker;
    else if (facet.pts[1]->id == 0U)
      facet.marker = facet.pts[0U]->marker;
    else
      facet.marker = side::BACK;
  }
  mesh.physicalNames["back"] = side::BACK;
}

void extrude(
    Mesh<Triangle> const & mesh2d,
    Mesh<Tetrahedron> & mesh3d,
    uint const numLayers,
    Vec3 const & direction,
    double const distance)
{
  auto const pt2dSize = mesh2d.pointList.size();
  mesh3d.pointList.reserve(pt2dSize * (numLayers + 1));
  auto const h = distance / numLayers;
  auto ptCounter = 0U;
  for (auto k = 0U; k < numLayers + 1; k++)
  {
    for (auto const & p: mesh2d.pointList)
    {
      mesh3d.pointList.emplace_back(p.coord + k * h * direction, ptCounter);
      ptCounter++;
    }
  }

  mesh3d.elementList.reserve(mesh2d.elementList.size() * numLayers * 3);
  std::vector<std::array<id_T, 3U>> childIds(mesh2d.elementList.size());
  auto elemCounter = 0U;
  for (auto k = 0U; k < numLayers; k++)
  {
    for (auto const & elem2d: mesh2d.elementList)
    {
      mesh3d.elementList.emplace_back(
          GeoElem::Pts_T{
              &mesh3d.pointList[elem2d.pts[0]->id + k * pt2dSize],
              &mesh3d.pointList[elem2d.pts[1]->id + k * pt2dSize],
              &mesh3d.pointList[elem2d.pts[2]->id + k * pt2dSize],
              &mesh3d.pointList[elem2d.pts[1]->id + (k + 1) * pt2dSize]},
          elemCounter,
          elem2d.marker);

      mesh3d.elementList.emplace_back(
          GeoElem::Pts_T{
              &mesh3d.pointList[elem2d.pts[0]->id + k * pt2dSize],
              &mesh3d.pointList[elem2d.pts[0]->id + (k + 1) * pt2dSize],
              &mesh3d.pointList[elem2d.pts[1]->id + (k + 1) * pt2dSize],
              &mesh3d.pointList[elem2d.pts[2]->id + k * pt2dSize]},
          elemCounter + 1U,
          elem2d.marker);

      mesh3d.elementList.emplace_back(
          GeoElem::Pts_T{
              &mesh3d.pointList[elem2d.pts[2]->id + k * pt2dSize],
              &mesh3d.pointList[elem2d.pts[0]->id + (k + 1) * pt2dSize],
              &mesh3d.pointList[elem2d.pts[1]->id + (k + 1) * pt2dSize],
              &mesh3d.pointList[elem2d.pts[2]->id + (k + 1) * pt2dSize]},
          elemCounter + 2U,
          elem2d.marker);

      childIds[elem2d.id] = {elemCounter, elemCounter + 1, elemCounter + 2};
      elemCounter += 3U;
    }
  }

  mesh3d.facetList.reserve(mesh2d.facetList.size() * 2U);
  mesh3d.elemToFacet.resize(
      mesh3d.elementList.size(),
      std::array<id_T, 4U>{idNotSet, idNotSet, idNotSet, idNotSet});
  auto facetCounter = 0U;
  // bottom
  for (auto const & elem2d: mesh2d.elementList)
  {
    mesh3d.facetList.emplace_back(
        GeoElem::Pts_T{
            &mesh3d.pointList[elem2d.pts[0]->id],
            &mesh3d.pointList[elem2d.pts[1]->id],
            &mesh3d.pointList[elem2d.pts[2]->id]},
        facetCounter,
        side::BOTTOM /*elem2d.marker*/);
    mesh3d.facetList.back().facingElem[0] = {
        &mesh3d.elementList[childIds[elem2d.id][0]], 0U};
    mesh3d.elemToFacet[mesh3d.elementList[childIds[elem2d.id][0]].id][0] = facetCounter;
    facetCounter++;
  }
  // lateral
  for (auto k = 0U; k < numLayers; k++)
  {
    for (auto const & facet2d: mesh2d.facetList)
    {
      auto childs2d = childIds[facet2d.facingElem[0].ptr->id];

      if (facet2d.facingElem[0].side % 2)
      {
        mesh3d.facetList.emplace_back(
            GeoElem::Pts_T{
                &mesh3d.pointList[facet2d.pts[0]->id + k * pt2dSize],
                &mesh3d.pointList[facet2d.pts[1]->id + k * pt2dSize],
                &mesh3d.pointList[facet2d.pts[0]->id + (k + 1) * pt2dSize]},
            facetCounter,
            facet2d.marker);
        mesh3d.facetList.back().facingElem[0] = {&mesh3d.elementList[childs2d[0]], 1U};
        mesh3d.elemToFacet[mesh3d.elementList[childs2d[0]].id][1] = facetCounter;

        mesh3d.facetList.emplace_back(
            GeoElem::Pts_T{
                &mesh3d.pointList[facet2d.pts[1]->id + k * pt2dSize],
                &mesh3d.pointList[facet2d.pts[1]->id + (k + 1) * pt2dSize],
                &mesh3d.pointList[facet2d.pts[0]->id + (k + 1) * pt2dSize]},
            facetCounter + 1U,
            facet2d.marker);
        mesh3d.facetList.back().facingElem[0] = {&mesh3d.elementList[childs2d[1]], 0U};
        mesh3d.elemToFacet[mesh3d.elementList[childs2d[1]].id][0] = facetCounter + 1U;
      }
      else
      {
        mesh3d.facetList.emplace_back(
            GeoElem::Pts_T{
                &mesh3d.pointList[facet2d.pts[0]->id + k * pt2dSize],
                &mesh3d.pointList[facet2d.pts[1]->id + k * pt2dSize],
                &mesh3d.pointList[facet2d.pts[1]->id + (k + 1) * pt2dSize]},
            facetCounter,
            facet2d.marker);
        mesh3d.facetList.back().facingElem[0] = {&mesh3d.elementList[childs2d[0]], 1U};
        mesh3d.elemToFacet[mesh3d.elementList[childs2d[0]].id][1] = facetCounter;

        mesh3d.facetList.emplace_back(
            GeoElem::Pts_T{
                &mesh3d.pointList[facet2d.pts[0]->id + k * pt2dSize],
                &mesh3d.pointList[facet2d.pts[1]->id + (k + 1) * pt2dSize],
                &mesh3d.pointList[facet2d.pts[0]->id + (k + 1) * pt2dSize]},
            facetCounter + 1U,
            facet2d.marker);
        mesh3d.facetList.back().facingElem[0] = {&mesh3d.elementList[childs2d[1]], 0U};
        mesh3d.elemToFacet[mesh3d.elementList[childs2d[1]].id][0] = facetCounter + 1U;
      }

      facetCounter += 2U;
    }
  }
  // top
  for (auto const & elem2d: mesh2d.elementList)
  {
    mesh3d.facetList.emplace_back(
        GeoElem::Pts_T{
            &mesh3d.pointList[elem2d.pts[0]->id + numLayers * pt2dSize],
            &mesh3d.pointList[elem2d.pts[1]->id + numLayers * pt2dSize],
            &mesh3d.pointList[elem2d.pts[2]->id + numLayers * pt2dSize]},
        facetCounter,
        side::TOP /*elem2d.marker + 100*/);
    mesh3d.facetList.back().facingElem[0] = {
        &mesh3d.elementList[childIds[elem2d.id][2]], 3U};
    mesh3d.elemToFacet[mesh3d.elementList[childIds[elem2d.id][2]].id][3] = facetCounter;
    facetCounter++;
  }

  mesh3d.buildConnectivity();
  buildFacets(mesh3d, mesh2d.flags);
  if (mesh2d.flags & MeshFlags::NORMALS)
  {
    buildNormals(mesh3d, mesh2d.flags);
  }
  if (mesh2d.flags & MeshFlags::FACET_PTRS)
  {
    addElemFacetList(mesh3d);
  }

  mesh3d.flags = mesh2d.flags;
  // TODO: redefine physical names
  mesh3d.physicalNames = mesh2d.physicalNames;
  mesh3d.physicalNames["bottom"] = side::BOTTOM;
  mesh3d.physicalNames["top"] = side::TOP;
}

void refTriangleMesh(Mesh<Triangle> & mesh)
{
  mesh.pointList = {
      Point(Vec3(0., 0., 0.), 0),
      Point(Vec3(1., 0., 0.), 1),
      Point(Vec3(0., 1., 0.), 2),
  };
  mesh.elementList = {
      Triangle{{&mesh.pointList[0], &mesh.pointList[1], &mesh.pointList[2]}, 0},
  };
  mesh.buildConnectivity();
  buildFacets(mesh, MeshFlags::INTERNAL_FACETS);
  addElemFacetList(mesh);
}

void refQuadMesh(Mesh<Quad> & mesh)
{
  mesh.pointList = {
      Point(Vec3(-1., -1., 0.), 0),
      Point(Vec3(1., -1., 0.), 1),
      Point(Vec3(1., 1., 0.), 2),
      Point(Vec3(-1., 1., 0.), 3),
  };
  mesh.elementList = {
      Quad{
          {&mesh.pointList[0],
           &mesh.pointList[1],
           &mesh.pointList[2],
           &mesh.pointList[3]},
          0},
  };
  mesh.buildConnectivity();
  buildFacets(mesh, MeshFlags::INTERNAL_FACETS);
  addElemFacetList(mesh);
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
      Tetrahedron{
          {&mesh.pointList[0],
           &mesh.pointList[1],
           &mesh.pointList[2],
           &mesh.pointList[3]},
          0},
  };
  mesh.buildConnectivity();
  buildFacets(mesh, MeshFlags::INTERNAL_FACETS);
  addElemFacetList(mesh);
}

void refHexahedronMesh(Mesh<Hexahedron> & mesh)
{
  mesh.pointList = {
      Point(Vec3(-1., -1., -1.), 0),
      Point(Vec3(1., -1., -1.), 1),
      Point(Vec3(1., 1., -1.), 2),
      Point(Vec3(-1., 1., -1.), 3),
      Point(Vec3(-1., -1., 1.), 4),
      Point(Vec3(1., -1., 1.), 5),
      Point(Vec3(1., 1., 1.), 6),
      Point(Vec3(-1., 1., 1.), 7),
  };
  mesh.elementList = {
      Hexahedron{
          {&mesh.pointList[0],
           &mesh.pointList[1],
           &mesh.pointList[2],
           &mesh.pointList[3],
           &mesh.pointList[4],
           &mesh.pointList[5],
           &mesh.pointList[6],
           &mesh.pointList[7]},
          0},
  };
  mesh.buildConnectivity();
  buildFacets(mesh, MeshFlags::INTERNAL_FACETS);
  addElemFacetList(mesh);
}

void hexagonMesh(Mesh<Triangle> & mesh)
{
  double const sr3o2 = 0.5 * sqrt(3.);
  mesh.pointList = {
      Point(Vec3(0., 0., 0.), 0),
      Point(Vec3(0., 1.0, 0.), 1),
      Point(Vec3(sr3o2, 0.5, 0.), 2),
      Point(Vec3(sr3o2, -0.5, 0.), 3),
      Point(Vec3(0., -1.0, 0.), 4),
      Point(Vec3(-sr3o2, -0.5, 0.), 5),
      Point(Vec3(-sr3o2, 0.5, 0.), 6),
  };
  mesh.elementList = {
      Triangle{{&mesh.pointList[0], &mesh.pointList[1], &mesh.pointList[2]}, 0},
      Triangle{{&mesh.pointList[0], &mesh.pointList[2], &mesh.pointList[3]}, 1},
      Triangle{{&mesh.pointList[0], &mesh.pointList[3], &mesh.pointList[4]}, 2},
      Triangle{{&mesh.pointList[0], &mesh.pointList[4], &mesh.pointList[5]}, 3},
      Triangle{{&mesh.pointList[0], &mesh.pointList[5], &mesh.pointList[6]}, 4},
      Triangle{{&mesh.pointList[0], &mesh.pointList[6], &mesh.pointList[1]}, 5},
  };
  mesh.buildConnectivity();
  buildFacets(mesh, MeshFlags::INTERNAL_FACETS);
  for (auto & facet: mesh.facetList)
  {
    if (facet.onBoundary())
    {
      facet.marker = 1;
      for (auto point: facet.pts)
      {
        point->marker = 1;
      }
    }
  }
  buildNormals(mesh, MeshFlags::NONE);
  // std::cout << mesh << std::endl;
}

void hexagonSquare(Mesh<Triangle> & mesh, MeshFlags flags)
{
  id_T counter = 0;
  double const h = 0.1;
  for (uint p = 0; p < 11; ++p) // id 1-10
  {
    mesh.pointList.emplace_back(Vec3{p * h, 0.0, 0.0}, counter++);
  }
  double const oneOnSr3 = 1. / std::sqrt(3.);
  uint const rows = 8;
  for (uint r = 0; r < rows; ++r)
  {
    for (uint p = 0; p < 5; ++p) // id 11-15
    {
      mesh.pointList.emplace_back(
          Vec3{(1 + 2 * p) * h, (1 + 2 * r) * h * oneOnSr3, 0.0}, counter++);
    }
    for (uint p = 0; p < 6; ++p) // id 16-21
    {
      mesh.pointList.emplace_back(
          Vec3{2 * p * h, 2 * (1 + r) * h * oneOnSr3, 0.0}, counter++);
    }
  }
  for (uint p = 0; p < 5; ++p) // id 22-26
  {
    mesh.pointList.emplace_back(
        Vec3{(1 + 2 * p) * h, (1 + 2 * rows) * h * oneOnSr3, 0.0}, counter++);
  }
  for (uint p = 0; p < 11; ++p) // id 27-37
  {
    mesh.pointList.emplace_back(
        Vec3{p * h, (2 + 2 * rows) * h * oneOnSr3, 0.0}, counter++);
  }

  counter = 0;
  for (uint e = 0; e < 5; ++e)
  {
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[0 + 2 * e],
             &mesh.pointList[1 + 2 * e],
             &mesh.pointList[11 + e]}},
        counter++);
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[1 + 2 * e],
             &mesh.pointList[2 + 2 * e],
             &mesh.pointList[11 + e]}},
        counter++);
  }
  for (uint e = 0; e < 5; ++e)
  {
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[0 + 2 * e],
             &mesh.pointList[11 + e],
             &mesh.pointList[16 + e]}},
        counter++);
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[2 + 2 * e],
             &mesh.pointList[11 + e],
             &mesh.pointList[17 + e]}},
        counter++);
  }
  for (uint r = 0; r < rows - 1; ++r)
  {
    for (uint e = 0; e < 5; ++e)
    {
      mesh.elementList.emplace_back(
          std::vector{
              {&mesh.pointList[16 + 11 * r + e],
               &mesh.pointList[11 + 11 * r + e],
               &mesh.pointList[22 + 11 * r + e]}},
          counter++);
      mesh.elementList.emplace_back(
          std::vector{
              {&mesh.pointList[17 + 11 * r + e],
               &mesh.pointList[11 + 11 * r + e],
               &mesh.pointList[22 + 11 * r + e]}},
          counter++);
    }
    for (uint e = 0; e < 5; ++e)
    {
      mesh.elementList.emplace_back(
          std::vector{
              {&mesh.pointList[16 + 11 * r + e],
               &mesh.pointList[22 + 11 * r + e],
               &mesh.pointList[27 + 11 * r + e]}},
          counter++);
      mesh.elementList.emplace_back(
          std::vector{
              {&mesh.pointList[17 + 11 * r + e],
               &mesh.pointList[22 + 11 * r + e],
               &mesh.pointList[28 + 11 * r + e]}},
          counter++);
    }
  }
  for (uint e = 0; e < 5; ++e)
  {
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[11 + 11 * (rows - 1) + e],
             &mesh.pointList[16 + 11 * (rows - 1) + e],
             &mesh.pointList[22 + 11 * (rows - 1) + e]}},
        counter++);
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[11 + 11 * (rows - 1) + e],
             &mesh.pointList[17 + 11 * (rows - 1) + e],
             &mesh.pointList[22 + 11 * (rows - 1) + e]}},
        counter++);
  }
  for (uint e = 0; e < 5; ++e)
  {
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[16 + 11 * (rows - 1) + e],
             &mesh.pointList[22 + 11 * (rows - 1) + e],
             &mesh.pointList[27 + 11 * (rows - 1) + 2 * e]}},
        counter++);
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[17 + 11 * (rows - 1) + e],
             &mesh.pointList[22 + 11 * (rows - 1) + e],
             &mesh.pointList[29 + 11 * (rows - 1) + 2 * e]}},
        counter++);
  }
  for (uint e = 0; e < 5; ++e)
  {
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[27 + 11 * (rows - 1) + 2 * e],
             &mesh.pointList[28 + 11 * (rows - 1) + 2 * e],
             &mesh.pointList[22 + 11 * (rows - 1) + e]}},
        counter++);
    mesh.elementList.emplace_back(
        std::vector{
            {&mesh.pointList[28 + 11 * (rows - 1) + 2 * e],
             &mesh.pointList[29 + 11 * (rows - 1) + 2 * e],
             &mesh.pointList[22 + 11 * (rows - 1) + e]}},
        counter++);
  }
  mesh.buildConnectivity();
  buildFacets(mesh, flags);
  markFacetsCube(mesh, {0., 0., 0.}, {1., (2 + 2 * rows) * h * oneOnSr3, 0.});
  buildNormals(mesh, MeshFlags::NONE);
}

uint checkPlanarFacets(Mesh<Hexahedron> const & mesh)
{
  auto nonPlanarCount = 0U;
  std::ranges::for_each(
      mesh.facetList,
      [&nonPlanarCount](auto const & facet)
      {
        // check diagonal 0-2
        auto const v10 = facet.pts[1]->coord - facet.pts[0]->coord;
        auto const v20 = facet.pts[2]->coord - facet.pts[0]->coord;

        // plane equation: normal.dot(p) = d
        Vec3 const normal012 = v10.cross(v20);
        double const d = normal012.dot(facet.pts[0]->coord);

        assert(std::fabs(normal012.dot(facet.pts[1]->coord) - d) < 1e-12);
        assert(std::fabs(normal012.dot(facet.pts[2]->coord) - d) < 1e-12);

        // check if forth point is on the plane defined by the other three
        if (std::fabs(normal012.dot(facet.pts[3]->coord) - d) > 1e-12)
        {
          nonPlanarCount++;
          fmt::print(stderr, "facet not planar along diagonal 0-2!\n");
          return;
        }

        // check diagonal 1-3
        auto const v30 = facet.pts[3]->coord - facet.pts[0]->coord;

        // plane equation: normal.dot(p) = d
        Vec3 const normal013 = v10.cross(v30);
        double const d013 = normal013.dot(facet.pts[0]->coord);

        assert(std::fabs(normal013.dot(facet.pts[1]->coord) - d013) < 1e-12);
        assert(std::fabs(normal013.dot(facet.pts[2]->coord) - d013) < 1e-12);

        // check if forth point is on the plane defined by the other three
        if (std::fabs(normal013.dot(facet.pts[3]->coord) - d013) > 1e-12)
        {
          nonPlanarCount++;
          fmt::print(stderr, "facet not planar along diagonal 1-3!\n");
        }
      });
  return nonPlanarCount;
}

} // namespace proxpde
