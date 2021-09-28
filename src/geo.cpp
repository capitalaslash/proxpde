#include "geo.hpp"

GeoElem::~GeoElem() = default;

std::ostream & operator<<(std::ostream & out, Point const & p)
{
  out << "(" << p[0] << "," << p[1] << "," << p[2] << "), id: " << p.id
      << ", m: " << p.marker;
  return out;
}

std::ostream & operator<<(std::ostream & out, GeoElem const & e)
{
  out << "pts: ";
  for (auto & p: e.pointList)
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

// children indices:
// 0 - 2 - 1
//   0   1
// clang-format off
std::array<FMat<2, 2>, 2> const Line::embeddingMatrix =
    std::array<FMat<2, 2>, 2>{{
        (FMat<2, 2>{} << 1.0, 0.0,             // 0
                         0.5, 0.5).finished(), // 2
        (FMat<2, 2>{} << 0.5, 0.5,             // 2
                         0.0, 1.0).finished(), // 1
    }};
// clang-format on

// children indices:
// 2
// | 2 +
// 5   -   4
// | 0 + 3 | 1 +
// 0   -   3   -   1
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

// children indices:
// 3 - 6 - 2
// | 3 | 2 |
// 7 - 8 - 5
// | 0 | 1 |
// 0 - 4 - 1
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
