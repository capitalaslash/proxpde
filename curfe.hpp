#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename RefFE, typename QR>
struct CurFE
{
  typedef RefFE RefFE_T;
  typedef QR QR_T;
  typedef Eigen::Matrix<double,RefFE::numFuns,RefFE::numFuns> LocalMat_T;
  typedef Eigen::Matrix<double,RefFE::numFuns,1> LocalVec_T;

  CurFE()
  {
    for(uint i=0; i<RefFE::numFuns; ++i)
    {
      for(uint q=0; q<QR::numPts; ++q)
      {
        phiRef(i,q) = RefFE::phiFun[i](QR::n[q]);
        dphiRef(i,q) = RefFE::dphiFun[i](QR::n[q]);
      }
    }
  }

  CurFE(CurFE const &) = delete;
  CurFE operator=(CurFE const &) = delete;

  void reinit(GeoElem const & elem)
  {
    uint const dim = RefFE::dim;
    uint const codim = 3-dim;
    for(uint q=0; q<QR::numPts; ++q)
    {
      J[q] = Eigen::Matrix3d::Zero();
      for(uint n=0; n<RefFE::numFuns; ++n)
      {
        J[q] += dofPts[n] * dphiRef(n,q).transpose();
      }
      // fill extra dimension with identity
      J[q].block(dim,dim,codim,codim) = Eigen::Matrix<double,codim,codim>::Identity();

      detJ[q] = J[q].determinant();
      JmT[q] = J[q].inverse().transpose();

      JxW[q] = detJ[q] * QR::w[q];
      qpoint[q] = elem.origin() + J[q] * QR::n[q];

      for(uint i=0; i<RefFE::numFuns; ++i)
      {
        phi(i,q) = phiRef(i,q);
        dphi(i,q) = JmT[q] * dphiRef(i,q);
      }
    }
  }

  vectorFun_T map;
  // vectorFun_T imap;
  std::array<Vec3,RefFE::numFuns> dofPts;
  std::array<Eigen::Matrix3d,QR::numPts> J;
  std::array<double,QR::numPts> detJ;
  std::array<Eigen::Matrix3d,QR::numPts> JmT;
  std::array<double,QR::numPts> JxW;
  std::array<Vec3,QR::numPts> qpoint;
  Eigen::Matrix<double, RefFE::numFuns, QR::numPts> phiRef;
  Eigen::Matrix<Vec3, RefFE::numFuns, QR::numPts> dphiRef;
  Eigen::Matrix<double, RefFE::numFuns, QR::numPts> phi;
  Eigen::Matrix<Vec3, RefFE::numFuns, QR::numPts> dphi;
  // std::array<scalarFun_T,RefFE::numFuns> phiFun;
  // std::array<vectorFun_T,RefFE::numFuns> dphiFun;
  LocalMat_T massMat;
  LocalMat_T stiffMat;
};
