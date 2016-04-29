#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename RefFE, typename QR>
struct CurFE
{
  typedef RefFE RefFE_T;
  typedef QR QR_T;
  typedef FMat<RefFE::numFuns,RefFE::numFuns> LocalMat_T;
  typedef FVec<RefFE::numFuns> LocalVec_T;
  typedef FMat<3,RefFE::dim> JacMat_T;
  typedef FMat<RefFE::dim,3> JacTMat_T;
  typedef typename RefFE_T::Vec_T RefVec_T;

  CurFE()
  {
    for(uint i=0; i<RefFE::numFuns; ++i)
    {
      for(uint q=0; q<QR::numPts; ++q)
      {
        phiRef[q](i) = RefFE::phiFun[i](QR::node[q]);
        dphiRef[q].col(i) = RefFE::dphiFun[i](QR::node[q]);
      }
    }
  }

  CurFE(CurFE const &) = delete;
  CurFE & operator=(CurFE const &) = delete;

  constexpr static uint size() {return RefFE::numFuns;}

  void reinit(GeoElem const & elem)
  {
    // e = &elem;
    dofPts = RefFE::dofPts(elem);

    for(uint q=0; q<QR::numPts; ++q)
    {
      jac[q] = JacMat_T::Zero();
      for(uint n=0; n<RefFE::numFuns; ++n)
      {
        jac[q] += dofPts[n] * dphiRef[q].col(n).transpose();
      }

      // J^+ = (J^T J)^-1 J^T
      auto jTj = jac[q].transpose() * jac[q];
      auto jTjI = jTj.inverse();
      detJ[q] = std::sqrt(jTj.determinant());
      jacMT[q] = jac[q] * jTjI.transpose();

      JxW[q] = detJ[q] * QR::weight[q];
      qpoint[q] = elem.origin() + jac[q] * QR::node[q];

      for(uint i=0; i<RefFE::numFuns; ++i)
      {
        // phi values on qpoints are unaffected by the change of coords
        // this update can potentially be done in the constructor
        phi[q](i) = phiRef[q](i);
        dphi[q].col(i) = jacMT[q] * dphiRef[q].col(i);
      }
    }
  }

  // GeoElem const* e;
  // vectorFun_T map;
  // vectorFun_T imap;
  std::array<Vec3,RefFE::numFuns> dofPts;
  std::array<JacMat_T,QR::numPts> jac;
  std::array<JacMat_T,QR::numPts> jacMT;
  std::array<double,QR::numPts> detJ;
  std::array<double,QR::numPts> JxW;
  std::array<Vec3,QR::numPts> qpoint;
  std::array<FVec<RefFE::numFuns>,QR::numPts> phiRef;
  std::array<FMat<RefFE::dim, RefFE::numFuns>,QR::numPts> dphiRef;
  std::array<FVec<RefFE::numFuns>,QR::numPts> phi;
  std::array<FMat<3, RefFE::numFuns>,QR::numPts> dphi;
  // LocalMat_T massMat;
  // LocalMat_T stiffMat;
};
