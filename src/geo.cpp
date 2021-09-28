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

