#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename RefFE, typename QR>
struct CurFE
{
  typedef Eigen::Matrix2d LocalMat_T;
  typedef RefFE RefFE_T;
  typedef QR QR_T;

  CurFE()
  {
    for(uint i=0; i<RefFE::numPts; ++i)
    {
      for(uint q=0; q<QR::numPts; ++q)
      {
        phi(i, q) = RefFE::phiFun[i](QR::n[q]);
        dphi(i, q) = RefFE::dphiFun[i](QR::n[q]);
      }
    }
  }

  void reinit(GeoElem const & elem)
  {
    J = elem.volume() / RefFE::volume;
    Jm1 = 1. / J;
    map = [this, &elem] (Vec3 const & p)
    {
      return elem.midpoint() + J * p;
    };
    imap = [this, &elem] (Vec3 const & p)
    {
      return (p - elem.midpoint()) * Jm1;
    };
    for(uint q=0; q<QR::numPts; ++q)
    {
      JxW[q] = J * QR::w[q];
      qpoint[q] = map(QR::n[q]);
    }
    massMat = RefFE::massMat * Jm1;
    stiffMat = RefFE::gradMat * Jm1;
    for(uint i=0; i<RefFE::numPts; ++i)
    {
      phiFun[i]  = [this, i] (Vec3 const & p)
      {
        return RefFE::phiFun[i](imap(p));
      };
      dphiFun[i] = [this, i] (Vec3 const & p)
      {
        return RefFE::dphiFun[i](imap(p)) * Jm1;
      };
    }
  }

  vectorFun_T map;
  vectorFun_T imap;
  double J;
  double Jm1;
  std::array<double,QR::numPts> JxW;
  std::array<Vec3,QR::numPts> qpoint;
  Eigen::Array<double, RefFE::numPts, QR::numPts> phi;
  Eigen::Array<Vec3, RefFE::numPts, QR::numPts> dphi;
  std::array<scalarFun_T,RefFE::numPts> phiFun;
  std::array<vectorFun_T,RefFE::numPts> dphiFun;
  LocalMat_T massMat;
  LocalMat_T stiffMat;
};
