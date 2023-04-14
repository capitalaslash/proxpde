#pragma once

#include "def.hpp"

#include "qr.hpp"
#include "reffe.hpp"

template <typename RefFE, typename QR>
struct CurFE
{
  using RefFE_T = RefFE;
  using QR_T = QR;
  static uint constexpr numDOFs = RefFE_T::numDOFs;
  static uint constexpr dim = RefFE_T::dim;
  using JacMat_T = FMat<3, dim>;
  using RefVec_T = typename RefFE_T::Vec_T;

  CurFE()
  {
    for (uint q = 0; q < QR::numPts; ++q)
    {
      for (uint i = 0; i < numDOFs; ++i)
      {
        phiRef[q](i) = RefFE::phiFun[i](QR::node[q]);
        dphiRef[q].row(i) = RefFE::dphiFun[i](QR::node[q]);

#ifdef MINIFEM_ENABLE_SECONDDERIV
        if constexpr (enableSecondDeriv_v<RefFE_T>)
        {
          for (uint d = 0; d < dim; ++d)
          {
            d2phiRef[q].row(i + d * numDOFs) =
                RefFE::d2phiFun[i + d * numDOFs](QR::node[q]);
          }
        }
#endif
      }
      for (uint i = 0; i < RefFE::numGeoDOFs; ++i)
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

      for (uint q = 0; q < QR::numPts; ++q)
      {
        jac[q] = JacMat_T::Zero();
        for (uint n = 0; n < RefFE::numGeoDOFs; ++n)
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

        // TODO: phi values on qpoints are unaffected by the change of coords
        // this update can potentially be done in the constructor
        phi[q] = phiRef[q];
        dphi[q] = dphiRef[q] * jacPlus[q];

#ifdef MINIFEM_ENABLE_SECONDDERIV
        if constexpr (enableSecondDeriv_v<RefFE_T>)
        {
          d2phi[q] = FMat<numDOFs * 3, 3>::Zero();
          for (uint d1 = 0; d1 < 3; ++d1)
          {
            for (uint d2 = 0; d2 < dim; ++d2)
            {
              d2phi[q].template block<numDOFs, 3>(d1 * numDOFs, 0) +=
                  jacPlus[q](d2, d1) *
                  d2phiRef[q].template block<numDOFs, dim>(d2 * numDOFs, 0) *
                  jacPlus[q];
            }
          }
        }
#endif
      }
    }
    // else
    // {
    //   std::cout << "warning: no reinit, element coincide" << std::endl;
    // }
  }

  Vec3 approxMap(FVec<dim> const & pt)
  {
    // jac is constant on the element only for linear mappings (not for bi-linear)
    return elem->origin() + jac[QR_T::bestPt] * pt;
  }

  Vec3 map(FVec<dim> const & pt)
  {
    JacMat_T jac = JacMat_T::Zero();
    auto const mappingPts = RefFE::mappingPts(*elem);
    for (uint n = 0; n < RefFE::numGeoDOFs; ++n)
    {
      jac += mappingPts[n] * (RefFE::mapping[n](pt)).transpose();
    }
    return elem->origin() + jac * pt;
  }

  FVec<dim> approxInverseMap(Vec3 const & pt)
  {
    // TODO: jacPlus should be computed on the point
    // here we assume that jacPlus does not change on the element
    // (this is true only for linear mappings)
    FVec<dim> const ptHat = jacPlus[QR_T::bestPt] * (pt - elem->origin());
    // check that the approximate inverse mapping is close enough to the
    // requested point
    assert((pt - approxMap(ptHat)).norm() < 1.e-14);
    return ptHat;
  }

  std::tuple<FVec<dim>, int> inverseMap(Vec3 const & pt, double const toll = 1.e-4)
  {
    FVec<dim> approxPt = approxInverseMap(pt);
    auto const mappingPts = RefFE::mappingPts(*elem);
    auto const origin = elem->origin();
    int iter = 0;
    static constexpr int maxIter = 30;
    while (iter < maxIter)
    {
      JacMat_T jac = JacMat_T::Zero();
      for (uint n = 0; n < RefFE::numGeoDOFs; ++n)
      {
        jac += mappingPts[n] * (RefFE::mapping[n](approxPt)).transpose();
      }
      Vec3 const predictedPt = origin + jac * approxPt;
      // filelog << iter << " " << std::setprecision(16) << approxPt.transpose() << " "
      // << (pt - predictedPt).norm() << std::endl;
      Vec3 const deltaReal = pt - predictedPt;
      auto const jacPlus = (jac.transpose() * jac).inverse() * jac.transpose();
      // Newton iteration
      FVec<dim> deltaRef = jacPlus * deltaReal;
      approxPt += deltaRef;
      iter++;
      if (!RefFE_T::inside(approxPt) || deltaRef.norm() < toll)
        break;
    }
    // if we reach the maximum number of iterations for a point inside the search failed
    assert(iter < maxIter || !RefFE_T::inside(approxPt));
    return std::tie(approxPt, iter);
  }

  bool approxInside(Vec3 const & pt) { return RefFE_T::inside(approxInverseMap(pt)); }

  bool inside(Vec3 const & pt) { return RefFE_T::inside(inverseMap(pt)); }

  GeoElem const * elem = nullptr;
  std::array<Vec3, numDOFs> dofPts;
  std::array<JacMat_T, QR::numPts> jac;
  std::array<Transpose_T<JacMat_T>, QR::numPts> jacPlus;
  std::array<double, QR::numPts> detJ;
  std::array<double, QR::numPts> JxW;
  std::array<Vec3, QR::numPts> qpoint;
  std::array<FVec<numDOFs>, QR::numPts> phiRef;
  std::array<FVec<numDOFs>, QR::numPts> phi;
  std::array<FMat<numDOFs, dim>, QR::numPts> dphiRef;
  std::array<FMat<numDOFs, 3>, QR::numPts> dphi;
#ifdef MINIFEM_ENABLE_SECONDDERIV
  std::array<FMat<numDOFs * dim, dim>, QR::numPts> d2phiRef;
  std::array<FMat<numDOFs * 3, 3>, QR::numPts> d2phi;
#endif
  std::array<FMat<RefFE::numGeoDOFs, dim>, QR::numPts> mapping;
};

template <typename RefFE, typename QR>
struct CurFE0D
{
  using RefFE_T = RefFE;

  void reinit(GeoElem const & e)
  {
    elem = &e;
    dofPts = RefFE_T::dofPts(*elem);

    qpoint = {{elem->origin()}};
  }

  GeoElem const * elem;
  std::array<Vec3, RefFE_T::numDOFs> dofPts;
  std::array<double, QR::numPts> JxW = {{1.L}};
  std::array<Vec3, QR::numPts> qpoint;
  std::array<FVec<1>, QR::numPts> phi = {{FVec<1>(1.L)}};
};

template <typename RefFE, typename QR>
struct VectorCurFE
{
  using RefFE_T = RefFE;
  using QR_T = QR;
  static uint constexpr numDOFs = RefFE_T::numDOFs;
  static uint constexpr dim = RefFE_T::dim;
  using JacMat_T = FMat<3, dim>;
  using RefVec_T = typename RefFE_T::Vec_T;

  VectorCurFE()
  {
    for (uint q = 0; q < QR::numPts; ++q)
    {
      for (uint i = 0; i < numDOFs; ++i)
      {
        phiVectRef[q].row(i) = RefFE::phiVectFun[i](QR::node[q]);
        divphiRef[q](i) = RefFE::divphiFun[i](QR::node[q]);
      }
      for (uint i = 0; i < RefFE::numGeoDOFs; ++i)
      {
        mapping[q].row(i) = RefFE::mapping[i](QR::node[q]);
      }
    }
    for (uint i = 0; i < RefFE::numGeoDOFs; ++i)
    {
      mapping0.row(i) = RefFE::mapping[i](FVec<dim>::Zero());
    }
  }

  // VectorCurFE(VectorCurFE const &) = delete;
  // VectorCurFE & operator=(VectorCurFE const &) = delete;

  void reinit(GeoElem const & e)
  {
    // no need to recompute everything if the element is the same
    // TODO: watch out for GeoElems from different classes
    if (!elem || elem->id != e.id)
    {
      elem = &e;

      for (uint f = 0; f < RefFE_T::GeoElem_T::numFacets; ++f)
      {
        // check if the internal element is the master one for the facet
        if (elem->facets[f]->facingElem[0].ptr->id == elem->id)
        {
          facetSign[f] = 1;
        }
        else
        {
          facetSign[f] = -1;
        }
      }

      dofPts = RefFE::dofPts(*elem);

      auto const mappingPts = RefFE::mappingPts(*elem);

      jac0 = JacMat_T::Zero();
      for (uint n = 0; n < RefFE::numGeoDOFs; ++n)
      {
        jac0 += mappingPts[n] * mapping0.row(n);
      }
      // J^+ = (J^T J)^-1 J^T
      auto const jTj0 = jac0.transpose() * jac0;
      detJ0 = std::sqrt(jTj0.determinant());
      double const detJ0Inv = 1. / detJ0;

      for (uint q = 0; q < QR::numPts; ++q)
      {
        jac[q] = JacMat_T::Zero();
        for (uint n = 0; n < RefFE::numGeoDOFs; ++n)
        {
          jac[q] += mappingPts[n] * mapping[q].row(n);
        }

        // J^+ = (J^T J)^-1 J^T
        auto const jTj = jac[q].transpose() * jac[q];
        detJ[q] = std::sqrt(jTj.determinant());
        auto const jTjI = jTj.inverse();
        jacPlus[q] = jTjI * jac[q].transpose();

        JxW[q] = detJ[q] * QR::weight[q];
        // forward mapping on qpoint
        qpoint[q] = elem->origin() + jac[q] * QR::node[q];

        // Piola transformation
        double detJInv = 1. / detJ[q];
        phiVect[q] = detJInv * phiVectRef[q] * jac[q].transpose();
        phiVect0[q] = detJ0Inv * phiVectRef[q] * jac0.transpose();
        divphi[q] = detJInv * divphiRef[q];
        divphi0[q] = detJ0Inv * divphiRef[q];
        // adjust signs based on normal going from lower to higher id
        for (uint f = 0; f < RefFE_T::GeoElem_T::numFacets; ++f)
        {
          phiVect[q].row(f) *= facetSign[f];
          phiVect0[q].row(f) *= facetSign[f];
          divphi[q](f) *= facetSign[f];
          divphi0[q](f) *= facetSign[f];
        }
      }
    }
    // else
    // {
    //   std::cout << "warning: no reinit, element coincide" << std::endl;
    // }
  }

  GeoElem const * elem = nullptr;
  std::array<int, RefFE_T::GeoElem_T::numFacets> facetSign;
  std::array<Vec3, numDOFs> dofPts;
  std::array<JacMat_T, QR::numPts> jac;
  // TODO: implement the method for general quadrangles with a dedicated MODIFIED_RT
  // CurFE
  JacMat_T jac0;
  std::array<Transpose_T<JacMat_T>, QR::numPts> jacPlus;
  std::array<double, QR::numPts> detJ;
  double detJ0;
  std::array<double, QR::numPts> JxW;
  std::array<Vec3, QR::numPts> qpoint;
  std::array<FMat<numDOFs, dim>, QR::numPts> phiVectRef;
  std::array<FMat<numDOFs, 3>, QR::numPts> phiVect;
  std::array<FMat<numDOFs, 3>, QR::numPts> phiVect0;
  std::array<FVec<numDOFs>, QR::numPts> divphiRef;
  std::array<FVec<numDOFs>, QR::numPts> divphi;
  std::array<FVec<numDOFs>, QR::numPts> divphi0;
  std::array<FMat<RefFE::numGeoDOFs, dim>, QR::numPts> mapping;
  FMat<RefFE::numGeoDOFs, dim> mapping0;
};

template <
    typename RefFE,
    typename QR,
    int dim = RefFE::dim,
    FamilyType F = family_v<RefFE>>
struct CurFETraits
{};

template <typename RefFE, typename QR>
struct CurFETraits<RefFE, QR, 0, FamilyType::NONE>
{
  using type = CurFE0D<RefFE, QR>;
};

template <typename RefFE, typename QR, int dim>
struct CurFETraits<RefFE, QR, dim, FamilyType::LAGRANGE>
{
  using type = CurFE<RefFE, QR>;
};

template <typename RefFE, typename QR, int dim>
struct CurFETraits<RefFE, QR, dim, FamilyType::CROUZEIX_RAVIART>
{
  using type = CurFE<RefFE, QR>;
};

template <typename RefFE, typename QR, int dim>
struct CurFETraits<RefFE, QR, dim, FamilyType::RAVIART_THOMAS>
{
  using type = VectorCurFE<RefFE, QR>;
};
