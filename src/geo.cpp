#include "geo.hpp"

// -------------------------------------------------------------------------------------
GeoElem::~GeoElem() = default;

// -------------------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & out, Point const & p)
{
  out << "(" << p[0] << "," << p[1] << "," << p[2] << "), id: " << p.id
      << ", m: " << p.marker;
  return out;
}

// -------------------------------------------------------------------------------------
bool operator==(GeoElem const & e1, GeoElem const & e2)
{
  if (!geoEqual(e1, e2))
  {
    return false;
  }

  for (short_T f = 0; f < e1.facets.size(); ++f)
  {
    if (e1.facets[f] != e2.facets[f])
    {
      return false;
    }
  }

  if (e1.id != e2.id)
  {
    return false;
  }

  if (e1.marker != e2.marker)
  {
    return false;
  }

  if (e1.parent != e2.parent)
  {
    return false;
  }

  for (short_T c = 0; c < e1.children.size(); ++c)
  {
    if (e1.children[c] != e2.children[c])
    {
      return false;
    }
  }

  for (short_T f = 0; f < 2; ++f)
  {
    if (e1.facingElem[f] != e2.facingElem[f])
    {
      return false;
    }
  }

  if (e1._normal != e2._normal)
  {
    return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & out, GeoElem const & e)
{
  out << "pts: ";
  for (auto & p: e.pts)
  {
    out << p->id << " ";
  }
  out << "id: " << e.id << ", m: " << e.marker;
  auto const [inPtr, inSide] = e.facingElem[0];
  if (inPtr)
  {
    out << ", fe: (" << inPtr->id << ", " << static_cast<int>(inSide) << ")";
    auto const [outPtr, outSide] = e.facingElem[1];
    if (outPtr)
    {
      out << ", (" << outPtr->id << ", " << static_cast<int>(outSide) << ")";
    }
    else
    {
      out << ", on boundary";
    }
  }
  if (e.parent.ptr)
  {
    out << ", parent id: " << e.parent.ptr->id;
  }
  else
  {
    out << ", no parent";
  }
  out << ", children ids (" << e.children.size() << "): ";
  for (auto const ch: e.children)
  {
    out << ch.ptr->id << " ";
  }
  return out;
}

// -------------------------------------------------------------------------------------
// children dofs:
// coarse: 0 - 1
// fine: 0 - 2 - 1
// child: 0:  0 - 2  1: 2 - 1
// clang-format off
std::array<FMat<2, 2>, 2> const Line::embeddingMatrix =
    std::array<FMat<2, 2>, 2>{{
        (FMat<2, 2>{} << 1.0, 0.0,             // 0
                         0.5, 0.5).finished(), // 2
        (FMat<2, 2>{} << 0.5, 0.5,             // 2
                         0.0, 1.0).finished(), // 1
    }};
// clang-format on

// -------------------------------------------------------------------------------------
// children dofs:
// coarse: 2
//         | +
//         0 - 1
// fine: 2
//       | +
//       5 - 4
//       | + | +
//       0 - 3 - 1
// child: 0: 5      1: 4      2: 2      3: 5 - 4
//           | +       | +       | +         + |
//           0 - 3     3 - 1     5 - 4         3
// TODO; this MUST be consistent with refinement ordering!!!
// clang-format off
std::array<FMat<3, 3>, 4> const Triangle::embeddingMatrix =
    std::array<FMat<3, 3>, 4>{{
        (FMat<3, 3>{} << 1.0, 0.0, 0.0,             // 0
                         0.5, 0.5, 0.0,             // 3
                         0.5, 0.0, 0.5).finished(), // 5
        (FMat<3, 3>{} << 0.0, 1.0, 0.0,             // 1
                         0.0, 0.5, 0.5,             // 4
                         0.5, 0.5, 0.0).finished(), // 3
        (FMat<3, 3>{} << 0.0, 0.0, 1.0,             // 2
                         0.5, 0.0, 0.5,             // 5
                         0.0, 0.5, 0.5).finished(), // 4
        (FMat<3, 3>{} << 0.5, 0.5, 0.0,             // 3
                         0.0, 0.5, 0.5,             // 4
                         0.5, 0.0, 0.5).finished(), // 5
    }};
// clang-format on

// -------------------------------------------------------------------------------------
// children dofs:
// coarse: 3 - 2
//         |   |
//         0 - 1
// fine: 3 - 6 - 2
//       |   |   |
//       7 - 8 - 5
//       |   |   |
//       0 - 4 - 1
// child: 0: 7 - 8  1: 8 - 5  2: 6 - 2  3: 3 - 6
//           |   |     |   |     |   |     |   |
//           0 - 4     4 - 1     8 - 5     7 - 8
// clang-format off
std::array<FMat<4, 4>, 4> const Quad::embeddingMatrix = std::array<FMat<4, 4>, 4>{{
    (FMat<4, 4>{} << 1.0,  0.0,  0.0,  0.0,              // 0
                     0.5,  0.5,  0.0,  0.0,              // 4
                     0.25, 0.25, 0.25, 0.25,             // 8
                     0.5,  0.0,  0.0,  0.5).finished(),  // 7
    (FMat<4, 4>{} << 0.5,  0.5,  0.0,  0.0,              // 4
                     0.0,  1.0,  0.0,  0.0,              // 1
                     0.0,  0.5,  0.5,  0.0,              // 5
                     0.25, 0.25, 0.25, 0.25).finished(), // 8
    (FMat<4, 4>{} << 0.25, 0.25, 0.25, 0.25,             // 8
                     0.0,  0.5,  0.5,  0.0,              // 5
                     0.0,  0.0,  1.0,  0.0,              // 2
                     0.0,  0.0,  0.5,  0.5).finished(),  // 6
    (FMat<4, 4>{} << 0.5,  0.0,  0.0,  0.5,              // 7
                     0.25, 0.25, 0.25, 0.25,             // 8
                     0.0,  0.0,  0.5,  0.5,              // 6
                     0.0,  0.0,  0.0,  1.0).finished(),  // 3                     
}};
// clang-format on

// -------------------------------------------------------------------------------------
bool geoEqual(GeoElem const & e1, GeoElem const & e2)
{
  assert(e1.pts.size() == e2.pts.size());
  std::vector<id_T> ids1(e1.pts.size()), ids2(e2.pts.size());
  uint counter = 0;
  std::for_each(
      e1.pts.begin(),
      e1.pts.end(),
      [&ids1, &counter](Point const * p) { ids1[counter++] = p->id; });
  counter = 0;
  std::for_each(
      e2.pts.begin(),
      e2.pts.end(),
      [&ids2, &counter](Point const * p) { ids2[counter++] = p->id; });

  return std::is_permutation(ids1.begin(), ids1.end(), ids2.begin());
}
