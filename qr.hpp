#pragma once

#include "def.hpp"
#include "geo.hpp"
#include "reffe.hpp"

double rectangleInt(Line const& e, id_T i, scalarFun_T rhs)
{
  // volume * shape fun midpoint * rhs fun midpoint
  return 2.*RefLineP1::phi[i](Vec3::Zero(3))*rhs(e.midpoint());
}

static Eigen::Vector2d const gauss2w = Eigen::Vector2d::Constant(1.L);
static double constexpr sqrt13rd = std::sqrt(1.L/3);
static std::array<Vec3,2> const gauss2n
{
  Vec3(-sqrt13rd, 0., 0.),
  Vec3( sqrt13rd, 0., 0.)
};

double gauss2Int(Line const& e, id_T i, scalarFun_T rhs)
{
  double sum = 0.;
  for(uint g=0; g<2; ++g)
  {
    Vec3 const node = e.midpoint() + 0.5*gauss2n[g]*e.volume();
    sum += gauss2w(g) * RefLineP1::phi[i](gauss2n[g]) * rhs(node);
  }
  return sum;
}

static Eigen::Vector3d const gauss3w = {5.L/9, 8.L/9, 5.L/9};
static double constexpr sqrt35th = std::sqrt(3.L/5);
static std::array<Vec3,3> const gauss3n
{
  Vec3(-sqrt35th, 0., 0.),
  Vec3(       0., 0., 0.),
  Vec3( sqrt35th, 0., 0.)
};

double gauss3Int(Line const& e, id_T i, scalarFun_T rhs)
{
  double sum = 0.;
  for(uint g=0; g<3; ++g)
  {
    Vec3 const node = e.midpoint() + 0.5*gauss3n[g]*e.volume();
    sum += gauss3w(g) * RefLineP1::phi[i](gauss3n[g]) * rhs(node);
  }
  return sum;
}
