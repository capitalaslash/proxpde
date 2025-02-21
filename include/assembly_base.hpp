#pragma once

#include "def.hpp"

#include "geo.hpp"

namespace proxpde
{

struct ScalarCoef
{
  explicit ScalarCoef(double const c): coef{c} {}

  void reinit(GeoElem const & /*elem*/) {}

  Vec1 evaluate(short_T const /*q*/) const { return Vec1{coef}; }

  double coef;
};

struct AssemblyBase
{
  // TODO: convert to std::unordered_set
  using CompList = std::vector<short_T>;

  CompList const comp = {};

  bool hasComp(short_T c) const
  {
    return std::find(comp.begin(), comp.end(), c) != comp.end();
  }
};

} // namespace proxpde
