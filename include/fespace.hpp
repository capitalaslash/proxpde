#pragma once

#include "def.hpp"

#include "curfe.hpp"
#include "dof.hpp"
#include "qr.hpp"

namespace proxpde
{

template <typename FESpace>
struct FEVar;

template <
    typename Mesh,
    typename RefFE,
    typename QR,
    uint Dimension = 1,
    DofType t = DofType::CONTINUOUS,
#ifdef PROXPDE_DOF_INTERLEAVED
    DofOrdering o = DofOrdering::INTERLEAVED>
#else
    DofOrdering o = DofOrdering::BLOCK>
#endif
struct FESpace
{
  using Mesh_T = Mesh;
  using RefFE_T = RefFE;
  using QR_T = QR;
  static constexpr DofType type = t;
  static constexpr DofOrdering ordering = o;
  using DOF_T = DOF<Mesh, RefFE, Dimension, type, ordering>;
  using CurFE_T = typename CurFETraits<RefFE, QR>::type;
  static short_T constexpr dim = Dimension;

  FESpace(Mesh const & m, uint const os): mesh{&m}, dof{m}, offset{os}
  {
    _compatibilityCheck();
  }

  explicit FESpace(Mesh const & m): FESpace{m, 0} {}

  FESpace() { _compatibilityCheck(); }

  void init(Mesh const & m, uint const os)
  {
    mesh = &m;
    dof.init(m);
    offset = os;
  }

  void init(Mesh const & m) { init(m, 0); }

  void _compatibilityCheck()
  {
    static_assert(
        std::is_same_v<typename Mesh_T::Elem_T, typename RefFE_T::GeoElem_T>,
        "mesh element and reference element are not compatible.");
    static_assert(
        std::is_same_v<typename RefFE_T::GeoElem_T, typename QR_T::GeoElem_T>,
        "reference element and quad rule are not compatible.");
    if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
    {
      static_assert(dim == 1, "Vector FESpaces can only be scalar on the dof.");
      // RT0 requires internal facets and facet ptrs
      assert(
          (mesh->flags & (MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS))
              .count() == 2);
    }
  }

  static uint constexpr physicalDim()
  {
    if constexpr (fedim_v<RefFE_T> == FEDimType::SCALAR)
    {
      return Dimension;
    }
    else if constexpr (fedim_v<RefFE_T> == FEDimType::VECTOR)
    {
      // return Mesh_T::Elem_T::dim;
      // when using RT elements, all physical structures are 3d
      return 3;
    }
  }

  static uint constexpr refDim()
  {
    if constexpr (fedim_v<RefFE_T> == FEDimType::SCALAR)
    {
      return Dimension;
    }
    else if constexpr (fedim_v<RefFE_T> == FEDimType::VECTOR)
    {
      return Mesh_T::Elem_T::dim;
    }
  }

  FVec<physicalDim()> evaluate(GeoElem const & elem, Vec const & data, Vec3 pt)
  {
    this->curFE.reinit(elem);

    // TODO: this needs inverse mapping
    // FMat<RefFE_T::numDOFs, RefFE_T::numDOFs> jac = FMat<RefFE_T::numDOFs,
    // RefFE_T::numDOFs>::Zero();
    // for(uint n=0; n<RefFE_T::numGeoDOFs; ++n)
    // {
    //   jac += curFE.dofPts[n] * RefFE::mapping[n](pt);
    // }

    FMat<RefFE_T::numDOFs, dim> localValue;
    FMat<RefFE_T::numDOFs, feDimValue<RefFE_T>()> phi;
    for (uint n = 0; n < CurFE_T::RefFE_T::numDOFs; ++n)
    {
      for (uint d = 0; d < dim; ++d)
      {
        id_T const dofId = this->dof.getId(elem.id, n, d);
        localValue(n, d) = data[dofId];
      }
      // TODO: jacPlus should be computed on the point
      // here we assume that jacPlus does not change on the element
      // (this is true only for affine mappings)
      auto const ptRef = this->curFE.jacPlus[QR_T::bestPt] * (pt - elem.origin());
      // check that the approximated inverse mapping is close enough
      assert(
          (pt - (elem.origin() + this->curFE.jac[QR_T::bestPt] * ptRef)).norm() <
          1.e-15);
      if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
      {
        // scalar functions do not change from reffe to curfe
        phi[n] = RefFE_T::phiFun[n](ptRef);
      }
      else
      {
        abort();
      }
      // instead of moving the point back to the ref element, we could
      // move the shape functions to the real element
    }

    FVec<physicalDim()> value;
    if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
    {
      value = localValue.transpose() * phi;
    }
    else
    {
      abort();
    }
    return value;
  }

  Vec3 findCoords(DOFid_T const id)
  {
    for (auto const & elem: mesh->elementList)
    {
      for (uint d = 0; d < DOF_T::clms; ++d)
      {
        if (dof.getId(elem.id, d) == id)
        {
          curFE.reinit(elem);
          return curFE.dofPts[d];
        }
      }
    }
    // we should never reach this point
    abort();
    return Vec3(0., 0., 0.);
  }

  Mesh const * mesh;
  CurFE_T mutable curFE;
  DOF_T dof;
  uint offset = 0U;
};

template <typename FESpace>
void interpolateAnalyticFunction(
    Fun<FESpace::physicalDim(), 3> const & fun, FESpace & feSpace, Vec & v)
{
  // set the vector data to the appropriate dimension if it comes with length 0
  if (v.size() == 0)
  {
    v = Vec::Zero(feSpace.dof.size * FESpace::dim);
  }

  if constexpr (family_v<typename FESpace::RefFE_T> == FamilyType::RAVIART_THOMAS)
  {
    assert((feSpace.mesh->flags & MeshFlags::NORMALS).any());
  }

  for (auto const & elem: feSpace.mesh->elementList)
  {
    feSpace.curFE.reinit(elem);
    for (uint i = 0; i < FESpace::RefFE_T::numDOFs; ++i)
    {
      auto const value = fun(feSpace.curFE.dofPts[i]);
      if constexpr (family_v<typename FESpace::RefFE_T> == FamilyType::LAGRANGE)
      {
        // the value of the dof is the value of the function
        // u_k = u phi_k
        // for vector fespaces this is a vector
        for (uint d = 0; d < FESpace::dim; ++d)
        {
          auto const dofId = feSpace.dof.getId(elem.id, i, d) + feSpace.offset;
          v[dofId] = value[d];
        }
      }
      else if constexpr (
          family_v<typename FESpace::RefFE_T> == FamilyType::RAVIART_THOMAS)
      {
        // the value of the dof is the flux through the facet
        // u_k = \int_{f_k} u.dot(n_k)
        // TODO: this should be an integral, we are taking the mean value at the
        // center of facet (ok for linear and piecewise constant functions only)
        id_T const facetId = feSpace.mesh->elemToFacet[elem.id][i];
        auto const & facet = feSpace.mesh->facetList[facetId];
        auto const dofId = feSpace.dof.getId(elem.id, i) + feSpace.offset;
        v[dofId] = value.dot(facet.normal()) * facet.volume();
      }
      else if constexpr (
          family_v<typename FESpace::RefFE_T> == FamilyType::CROUZEIX_RAVIART)
      {
        // the value of the dof is the integral on the opposite facet
        // u = \sum_i u_i \phi^CR_i
        // \frac{1}{|f_k|} \int_{f_k} u ds =
        // \frac{1}{|f_k|} \int_{f_k} sum_i u_i \phi^CR_i ds =
        // \sum_i \frac{u_i}{|f_k|} \int_{f_k} \phi^CR_i ds =
        // \frac{u_k}{|f_k|} \int_{f_k} \phi^CR_k
        // => u_k = \frac{1}{|f_k|} \int_{f_k} u ds
        // TODO: do a real integral.
        // approximating it with the mean value changes the element to a a Lagrange type
        // with different interpolating space (see Guermond TAMU 2015 ch. 3)
        for (uint d = 0; d < FESpace::dim; ++d)
        {
          auto const dofId = feSpace.dof.getId(elem.id, i, d) + feSpace.offset;
          v[dofId] = value[d];
        }
      }
      else
      {
        std::abort();
      }

      // auto const value = evaluate(f, feSpace, i);
      // for (uint d=0; d<FESpace::dim; ++d)
      // {
      //   auto const dofId = feSpace.dof.getId(e.id, i, d);
      //   v[offset + dofId] = value[d];
      // }
    }
  }
}

template <typename FESpace>
void interpolateAnalyticFunction(scalarFun_T const & f, FESpace & feSpace, Vec & v)
{
  interpolateAnalyticFunction([f](Vec3 const & p) { return Vec1(f(p)); }, feSpace, v);
}

template <typename FESpace>
void integrateAnalyticFunction(
    Fun<FESpace::dim, 3> const & f, FESpace & feSpace, Vec & v)
{
  auto & curFE = feSpace.curFE;

  // set the vector data to the appropriate dimension if it comes with length 0
  if (v.size() == 0)
  {
    v = Vec::Zero(feSpace.dof.size * FESpace::dim);
  }

  for (auto const & e: feSpace.mesh->elementList)
  {
    curFE.reinit(e);

    // avoid costly integration when 0 or 1
    FVec<FESpace::dim> sum = FVec<FESpace::dim>::Zero();
    for (auto const & pt: e.pts)
    {
      sum += f(pt->coord);
    }
    sum /= FESpace::Mesh_T::Elem_T::numPts;

    for (uint i = 0; i < FESpace::RefFE_T::numDOFs; ++i)
    {
      FVec<FESpace::dim> value = FVec<FESpace::dim>::Zero();
      double const sum_norm = sum.norm();
      if (sum_norm < 1.e-15 || sum_norm > 1.0 - 1.e-15)
      {
        value = sum;
      }
      else
      {
        for (uint q = 0; q < FESpace::QR_T::numPts; ++q)
        {
          value += curFE.JxW[q] * f(curFE.qpoint[q]);
        }
        // TODO: this is required only for color functions, should be optional
        value /= e.volume();
      }

      // set the value in the output vector
      for (uint d = 0; d < FESpace::dim; ++d)
      {
        auto const dofId = feSpace.dof.getId(e.id, i, d);
        if constexpr (family_v<typename FESpace::RefFE_T> == FamilyType::LAGRANGE)
        {
          // the value of the dof is the value of the function
          // u_k = u phi_k
          v[feSpace.offset + dofId] = value[d];
        }
        else if constexpr (
            family_v<typename FESpace::RefFE_T> == FamilyType::RAVIART_THOMAS)
        {
          std::abort();
          // the value of the dof is the flux through the face
          // u_k = u.dot(n_k)
          v[feSpace.offset + dofId] = value.dot(FESpace::RefFE_T::normal(e)[i]);
        }
        else
        {
          std::abort();
        }
      }
    }
  }
}

template <typename FESpace>
void integrateAnalyticFunction(scalarFun_T const & f, FESpace & feSpace, Vec & v)
{
  integrateAnalyticFunction([f](Vec3 const & p) { return Vec1(f(p)); }, feSpace, v);
}

template <typename FESpaceData, typename FESpaceGrad>
void reconstructGradient(
    Vec & grad,
    FESpaceGrad & feSpaceGrad,
    Vec const & data,
    FESpaceData & feSpaceData,
    std::vector<short_T> const & comp = allComp<FESpaceGrad>(),
    uint const offset = 0)
{
  auto const size = feSpaceData.dof.size;
  assert(data.size() == size);
  assert(FESpaceGrad::dim == comp.size());
  grad = Vec::Zero(feSpaceGrad.dof.size * FESpaceGrad::dim);
  Eigen::VectorXi numberOfPasses = Eigen::VectorXi::Zero(size);
  for (auto const & elem: feSpaceData.mesh->elementList)
  {
    feSpaceData.curFE.reinit(elem);
    feSpaceGrad.curFE.reinit(elem);

    for (uint k = 0; k < FESpaceData::RefFE_T::numDOFs; ++k)
    {
      auto const baseDof = feSpaceData.dof.getId(elem.id, k);
      auto const value = data[baseDof];
      numberOfPasses[baseDof] += 1;
      for (uint i = 0; i < FESpaceGrad::RefFE_T::numDOFs; ++i)
      {
        for (auto const d: comp)
        {
          if constexpr (family_v<typename FESpaceData::RefFE_T> == FamilyType::LAGRANGE)
          {
            // grad u (x_i) = sum_k u_k dphi_k (x_i)
            grad[offset + feSpaceGrad.dof.getId(elem.id, i, d)] +=
                value * feSpaceGrad.curFE.dphi[i](k, d);
          }
          else if constexpr (
              family_v<typename FESpaceData::RefFE_T> == FamilyType::RAVIART_THOMAS)
          {
            abort();
            // the value of the dof is the flux through the face
            // u_k = u.dot(n_k)
            // v[offset + baseDof + d*size] =
            // value[d].dot(FESpace::RefFE_T::normal(e)[i]);
          }
        }
      }
    }
  }
  // divide by the number of elements which provided a gradient to get the mean value of
  // the gradient
  for (uint i = 0; i < size; ++i)
  {
    for (auto const d: comp)
    {
      grad[i + d * size] /= numberOfPasses[i];
    }
  }
}

template <typename FESpaceDest, typename FESpaceOrig>
void getComponent(
    Vec & dest,
    FESpaceDest const & feSpaceDest,
    Vec const & orig,
    [[maybe_unused]] FESpaceOrig const & feSpaceOrig,
    uint component)
{
  static_assert(FESpaceOrig::dim > 1, "should not be used on scalar fespaces");
  auto const size = feSpaceDest.dof.size;
  auto constexpr dim = FESpaceOrig::dim;
  assert(size == feSpaceOrig.dof.size);
  if (FESpaceOrig::DOF_T::ordering == DofOrdering::BLOCK)
  {
    dest = orig.segment(component * size, size);
  }
  else // FESpaceVel_T::DOF_T::ordering == DofOrdering::INTERLEAVED
  {
    dest = orig(Eigen::seqN(component, size, dim));
  }
}

template <typename FESpaceOrig>
void getComponents(
    std::array<Vec &, FESpaceOrig::dim> dest,
    Vec & orig,
    FESpaceOrig const & feSpaceOrig)
{
  using FESpaceDest = FESpace<
      typename FESpaceOrig::Mesh_T,
      typename FESpaceOrig::RefFE_T,
      typename FESpaceOrig::QR_T,
      1>;
  FESpaceDest feSpaceDest{*feSpaceOrig.mesh};
  for (uint d = 0; d < feSpaceOrig; ++d)
  {
    getComponent(dest[d], feSpaceDest, orig, feSpaceOrig, d);
  }
}

template <typename FESpaceDest, typename FESpaceOrig>
void setComponent(
    Vec & dest,
    [[maybe_unused]] FESpaceDest const & feSpaceDest,
    Vec const & orig,
    FESpaceOrig const & feSpaceOrig,
    uint component)
{
  static_assert(FESpaceDest::dim > 1, "should not be used on scalar fespaces");
  auto const size = feSpaceOrig.dof.size;
  auto constexpr dim = FESpaceDest::dim;
  assert(size == feSpaceDest.dof.size);
  assert(component < dim);
  if (FESpaceDest::DOF_T::ordering == DofOrdering::BLOCK)
  {
    dest.segment(component * size, size) = orig;
  }
  else // FESpaceVel_T::DOF_T::ordering == DofOrdering::INTERLEAVED
  {
    for (uint k = 0; k < size; ++k)
    {
      dest[k * dim + component] = orig[k];
    }
  }
}

template <typename FESpaceT>
double integrateOnBoundary(Vec const & u, FESpaceT const & feSpace, marker_T const m)
{
  using FESpace_T = FESpaceT;
  using Mesh_T = typename FESpace_T::Mesh_T;
  using Elem_T = typename Mesh_T::Elem_T;
  using RefFE_T = typename FESpace_T::RefFE_T;
  using FacetFE_T = typename RefFE_T::FacetFE_T;
  using FacetQR_T = SideQR_T<typename FESpace_T::QR_T>;
  using FacetCurFE_T = typename CurFETraits<FacetFE_T, FacetQR_T>::type;
  using FacetFESpace_T =
      FESpace<Mesh_T, RefFE_T, SideGaussQR<Elem_T, FacetQR_T::numPts>>;

  FacetCurFE_T facetCurFE;
  FacetFESpace_T facetFESpace{*feSpace.mesh};
  FEVar uFacet{facetFESpace};
  uFacet.data = u;

  double integral = 0.;
  for (auto & facet: feSpace.mesh->facetList)
  {
    if (facet.marker == m)
    {
      // std::cout << "facet " << facet.id << std::endl;
      facetCurFE.reinit(facet);
      auto elem = facet.facingElem[0].ptr;
      auto const side = facet.facingElem[0].side;
      uFacet.reinit(*elem);
      for (uint q = 0; q < FacetQR_T::numPts; ++q)
      {
        double value;
        if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
        {
          // TODO: this assumes that the variable is scalar
          value = uFacet.evaluate(side * FacetQR_T::numPts + q)[0];
        }
        else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
        {
          value = uFacet.evaluate(side * FacetQR_T::numPts + q).dot(facet._normal);
        }
        integral += facetCurFE.JxW[q] * value;
      }
    }
  }
  return integral;
}

template <typename FESpaceVec>
using Scalar_T = FESpace<
    typename FESpaceVec::Mesh_T,
    typename FESpaceVec::RefFE_T,
    typename FESpaceVec::QR_T,
    1,
    FESpaceVec::type,
    FESpaceVec::ordering>;

} // namespace proxpde
