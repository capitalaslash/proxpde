#pragma once

#include "def.hpp"
#include "reffe.hpp"
#include "qr.hpp"

template <typename RefFE, typename QR>
struct CurFE
{
  using RefFE_T = RefFE;
  using QR_T = QR;
  static uint constexpr numDOFs = RefFE::numFuns;
  using LocalMat_T = FMat<numDOFs,numDOFs>;
  using LocalVec_T = FVec<numDOFs>;
  using JacMat_T = FMat<3,RefFE::dim>;
  using JacTMat_T = FMat<RefFE::dim,3>;
  using RefVec_T = typename RefFE_T::Vec_T;

  CurFE()
  {
    for(uint q=0; q<QR::numPts; ++q)
    {
      for(uint i=0; i<RefFE::numFuns; ++i)
      {
        if constexpr (fedim_v<RefFE_T> == FEDimType::SCALAR)
        {
          phiRef[q](i) = RefFE::phiFun[i](QR::node[q]);
          dphiRef[q].row(i) = RefFE::dphiFun[i](QR::node[q]);
        }
        else if constexpr (fedim_v<RefFE_T> == FEDimType::VECTOR)
        {
          phiVectRef[q].row(i) = RefFE::phiVectFun[i](QR::node[q]);
          divphiRef[q](i) = RefFE::divphiFun[i](QR::node[q]);
        }
      }
      for(uint i=0; i<RefFE::numGeoFuns; ++i)
      {
        mapping[q].row(i) = RefFE::mapping[i](QR::node[q]);
      }
    }
  }

  // CurFE(CurFE const &) = delete;
  // CurFE & operator=(CurFE const &) = delete;

  void reinit(GeoElem const & e)
  {
    // no need to recompute everything if the element is the same
    // TODO: watch out for GeoElems from different classes
    if (!elem || elem->id != e.id)
    {
      elem = &e;
      dofPts = RefFE::dofPts(*elem);

      auto const mappingPts = RefFE::mappingPts(*elem);

      for(uint q=0; q<QR::numPts; ++q)
      {
        jac[q] = JacMat_T::Zero();
        for(uint n=0; n<RefFE::numGeoFuns; ++n)
        {
          jac[q] += mappingPts[n] * mapping[q].row(n);
        }

        // J^+ = (J^T J)^-1 J^T
        auto const jTj = jac[q].transpose() * jac[q];
        auto const jTjI = jTj.inverse();
        detJ[q] = std::sqrt(jTj.determinant());
        jacPlus[q] = jTjI * jac[q].transpose();

        JxW[q] = detJ[q] * QR::weight[q];
        // forward mapping on qpoint
        qpoint[q] = elem->origin() + jac[q] * QR::node[q];

        if constexpr (fedim_v<RefFE_T> == FEDimType::SCALAR)
        {
          // TODO: phi values on qpoints are unaffected by the change of coords
          // this update can potentially be done in the constructor
          phi[q] = phiRef[q];
          dphi[q] = dphiRef[q] * jacPlus[q];
        }
        else if constexpr (fedim_v<RefFE_T> == FEDimType::VECTOR)
        {
          // Piola transformation
          double detJInv = 1. / detJ[q];
          phiVect[q] = detJInv * phiVectRef[q] * jac[q].transpose();
          divphi[q] = detJInv * divphiRef[q];
          // adjust signs based on normal going from lower id to greater id
          for (uint f=0; f<RefFE_T::GeoElem_T::numFacets; ++f)
          {
            if (elem->facetList[f]->facingElem[0].ptr->id != elem->id)
            {
              phiVect[q].row(f) *= -1.;
              divphi[q](f) *= -1.;
            }
          }
        }
      }
    }
    // else
    // {
    //   std::cout << "warning: no reinit, element coincide" << std::endl;
    // }
  }

  Vec3 approxMap(FVec<RefFE::dim> const & pt)
  {
    // jac is constant on the element only for linear mappings (not for bi-linear)
    return elem->origin() + jac[QR_T::bestPt] * pt;
  }

  Vec3 map(FVec<RefFE::dim> const & pt)
  {
    JacMat_T jac = JacMat_T::Zero();
    auto const mappingPts = RefFE::mappingPts(*elem);
    for(uint n=0; n<RefFE::numGeoFuns; ++n)
    {
      jac += mappingPts[n] * (RefFE::mapping[n](pt)).transpose();
    }
    return elem->origin() + jac * pt;
  }

  FVec<RefFE::dim> approxInverseMap(Vec3 const & pt)
  {
    // TODO: jacPlus should be computed on the point
    // here we assume that jacPlus does not change on the element
    // (this is true only for linear mappings)
    FVec<RefFE::dim> const ptHat = jacPlus[QR_T::bestPt] * (pt - elem->origin());
    // check that the approximate inverse mapping is close enough to the
    // requested point
    assert((pt - approxMap(ptHat)).norm() < 1.e-14);
    return ptHat;
  }

  FVec<RefFE::dim> inverseMap(Vec3 const & pt)
  {
    FVec<RefFE::dim> approxPt = approxInverseMap(pt);
    auto const mappingPts = RefFE::mappingPts(*elem);
    auto const origin = elem->origin();
    int iter = 0;
    while (iter < 200)
    {
      JacMat_T jac = JacMat_T::Zero();
      for(uint n=0; n<RefFE::numGeoFuns; ++n)
      {
        jac += mappingPts[n] * (RefFE::mapping[n](approxPt)).transpose();
      }
      auto const predictedPt = origin + jac * approxPt;
      // filelog << iter << " " << std::setprecision(16) << approxPt.transpose() << " " << (pt - predictedPt).norm() << std::endl;
      // TODO: this error is computed on the real elem, maybe move it to reference elem using detJac?
      if (!RefFE_T::inside(approxPt) || (pt - predictedPt).norm() < 3.e-15)
        break;
      auto const jacPlus = (jac.transpose() * jac).inverse() * jac.transpose();
      // fixed-point iteration
      approxPt = jacPlus * (pt - origin);
      iter++;
    }
    // if we reach the maximum number of iterations for a point inside the search failed
    assert(iter != 200 || !RefFE_T::inside(approxPt));
    return approxPt;
  }

  bool approxInside(Vec3 const & pt)
  {
    return RefFE_T::inside(approxInverseMap(pt));
  }

  bool inside(Vec3 const & pt)
  {
    return RefFE_T::inside(inverseMap(pt));
  }

  GeoElem const * elem = nullptr;
  array<Vec3,RefFE::numFuns> dofPts;
  array<JacMat_T,QR::numPts> jac;
  array<JacTMat_T,QR::numPts> jacPlus;
  array<double,QR::numPts> detJ;
  array<double,QR::numPts> JxW;
  array<Vec3,QR::numPts> qpoint;
  array<FVec<RefFE::numFuns>,QR::numPts> phiRef;
  array<FMat<RefFE::numFuns,RefFE::dim>,QR::numPts> dphiRef;
  array<FMat<RefFE::numFuns,RefFE::dim>,QR::numPts> phiVectRef;
  array<FVec<RefFE::numFuns>,QR::numPts> divphiRef;
  array<FMat<RefFE::numGeoFuns,RefFE::dim>,QR::numPts> mapping;
  // vectorFun_T inverseMapping;
  array<FVec<RefFE::numFuns>,QR::numPts> phi;
  array<FMat<RefFE::numFuns,3>,QR::numPts> dphi;
  // TODO: box phiVect in external data struct that is specialized on the FEDimType
  array<FMat<RefFE::numFuns,3>,QR::numPts> phiVect;
  array<FVec<RefFE::numFuns>,QR::numPts> divphi;
};

template <typename QR>
struct CurFE<RefPointP1,QR>
{
  using RefFE_T = RefPointP1;

  void reinit(GeoElem const & e)
  {
    elem = &e;
    dofPts = RefFE_T::dofPts(*elem);

    qpoint = {{elem->origin()}};
  }

  GeoElem const * elem;
  array<Vec3,RefFE_T::numFuns> dofPts;
  array<double,QR::numPts> JxW = {{1.L}};
  array<Vec3,QR::numPts> qpoint;
  array<FVec<1>,QR::numPts> phi = {{FVec<1>(1.L)}};
};
