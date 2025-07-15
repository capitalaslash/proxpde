#pragma once

#include "def.hpp"

// stl
#include <filesystem>

// local
#include "fespace.hpp"
#include "qr.hpp"

namespace proxpde
{

template <typename FESpace>
struct IOManager;

// ---------------------------------------------------------------------
struct Var
{
  explicit Var(std::string_view const n = "none", unsigned long const size = 0):
      name{n},
      data(size)
  {}

  Var(std::string_view const n, Vec const vec): name{n}, data(std::move(vec)) {}

  Var(std::string_view const n, Vec const & vec, uint offset, uint size):
      name{n},
      data(vec.segment(offset, size))
  {}

  std::string name;
  Vec data;
};

std::vector<uint> offsetInit(std::vector<uint> blocks);

// ---------------------------------------------------------------------
struct BlockVar: public Var
{
  BlockVar(std::string_view const n, std::vector<uint> const & bs):
      Var{n, std::accumulate(bs.begin(), bs.end(), 0U)},
      offsets{offsetInit(bs)},
      blocks{bs}
  {}

  auto block(uint i) { return this->data.segment(offsets[i], blocks[i]); }

  std::vector<uint> const offsets;
  std::vector<uint> const blocks;
};

// ---------------------------------------------------------------------
template <typename FESpace>
struct FEVar
{
  using FESpace_T = FESpace;
  using RefFE_T = typename FESpace_T::RefFE_T;
  using Vec_T = FVec<FESpace_T::physicalDim()>;

  FEVar(std::string_view n, FESpace_T const & fe):
      name{n},
      feSpace{&fe},
      data{fe.dof.size * FESpace_T::dim}
  {}

  explicit FEVar(FESpace_T const & fe): FEVar("none", fe) {}

  FEVar() = default;

  void init(std::string_view n, FESpace_T const & fe)
  {
    name = n;
    feSpace = &fe;
    data.resize(fe.dof.size * FESpace_T::dim);
  }

  FEVar<FESpace_T> & operator<<(Fun<FESpace_T::physicalDim(), 3> const & f)
  {
    interpolateAnalyticFunction(f, *this->feSpace, this->data);
    return *this;
  }

  FEVar<FESpace_T> & operator<<(scalarFun_T const & f)
  {
    static_assert(FESpace_T::dim == 1, "this is available only on scalar fe spaces");
    return operator<<([f](Vec3 const & p) { return Vec1::Constant(f(p)); });
  }

  FEVar<FESpace_T> & operator<<(Vec_T const & v)
  {
    return operator<<([&v](Vec3 const &) { return v; });
  }

  FEVar<FESpace_T> & operator<<(double const v)
  {
    static_assert(FESpace_T::dim == 1, "this is available only on scalar fe spaces");
    return operator<<([&v](Vec3 const &) { return Vec1::Constant(v); });
  }

  FEVar<FESpace_T> & operator>>(std::filesystem::path const & path)
  {
    IOManager io{*feSpace, path};
    io.print({*this});
    return *this;
  }

  double operator[](id_T const id) const { return data[id]; }

  void reinit(GeoElem const & elem)
  {
    feSpace->curFE.reinit(elem);
    setLocalData(elem.id);
  }

  void setLocalData(id_T const elemId)
  {
    for (uint d = 0; d < FESpace_T::dim; ++d)
    {
      for (uint n = 0; n < FESpace_T::RefFE_T::numDOFs; ++n)
      {
        auto const id = feSpace->dof.getId(elemId, n, d);
        dataLocal(n, d) = data[id];
      }
    }
  }

  auto getFacetMeanValue(uint const side) const
  {
    assert(side < FESpace_T::Mesh_T::Elem_T::numFacets);
    FVec<FESpace_T::dim> sum = FVec<FESpace_T::dim>::Zero();
    for (auto const dofFacet: RefFE_T::dofOnFacet[side])
    {
      sum += dataLocal.row(dofFacet).transpose();
    }
    sum /= RefFE_T::dofPerFacet;
    return sum;
  }

  auto integrate()
  {
    FVec<FESpace_T::dim> integral = FVec<FESpace_T::dim>::Zero();
    for (auto const & elem: feSpace->mesh->elementList)
    {
      this->reinit(elem);
      for (uint q = 0; q < FESpace_T::QR_T::numPts; ++q)
      {
        integral += feSpace->curFE.JxW[q] * this->evaluate(q);
      }
    }
    return integral;
  }

  auto l2NormSquared()
  {
    double integral = 0.0;
    for (auto const & elem: feSpace->mesh->elementList)
    {
      this->reinit(elem);
      for (uint q = 0; q < FESpace_T::QR_T::numPts; ++q)
      {
        auto const localValue = this->evaluate(q);
        integral += feSpace->curFE.JxW[q] * localValue.dot(localValue);
      }
    }
    return integral;
  }

  auto l2ErrorSquared(Fun<FESpace_T::physicalDim(), 3> const & exact)
  {
    double error = 0.0;
    for (auto const & elem: feSpace->mesh->elementList)
    {
      this->reinit(elem);
      for (uint q = 0; q < FESpace_T::QR_T::numPts; ++q)
      {
        auto const localValue = this->evaluate(q);
        auto const exactQ = exact(feSpace->curFE.qpoint[q]);
        error +=
            feSpace->curFE.JxW[q] * ((localValue - exactQ).dot(localValue - exactQ));
      }
    }
    return error;
  }

  auto l2ErrorSquared(scalarFun_T const & exact)
  {
    return l2ErrorSquared([exact](Vec3 const & p) { return Vec1{exact(p)}; });
  }

  auto l2ErrorSquared0(Fun<FESpace_T::physicalDim(), 3> const & exact)
  {
    double error = 0.0;
    for (auto const & elem: feSpace->mesh->elementList)
    {
      this->reinit(elem);
      for (uint q = 0; q < FESpace_T::QR_T::numPts; ++q)
      {
        auto const localValue = this->evaluate0(q);
        auto const exactQ = exact(feSpace->curFE.qpoint[q]);
        error +=
            feSpace->curFE.JxW[q] * ((localValue - exactQ).dot(localValue - exactQ));
      }
    }
    return error;
  }

  auto h1SemiNormSquared()
  {
    double integral = 0.0;
    for (auto const & elem: feSpace->mesh->elementList)
    {
      this->reinit(elem);
      for (uint q = 0; q < FESpace_T::QR_T::numPts; ++q)
      {
        auto const localGrad = this->evaluateGrad(q);
        integral += feSpace->curFE.JxW[q] * (localGrad.transpose() * localGrad).trace();
      }
    }
    return integral;
  }

  auto divL2ErrorSquared(scalarFun_T const & exact)
  {
    double error = 0.0;
    for (auto const & elem: feSpace->mesh->elementList)
    {
      this->reinit(elem);
      for (uint q = 0; q < FESpace_T::QR_T::numPts; ++q)
      {
        auto const localValue = this->evaluateDiv(q);
        auto const exactQ = exact(feSpace->curFE.qpoint[q]);
        error +=
            feSpace->curFE.JxW[q] * ((localValue - exactQ) * (localValue - exactQ));
      }
    }
    return error;
  }

  auto DIVL2ErrorSquared(scalarFun_T const & exact)
  {
    double error = 0.0;
    for (auto const & elem: feSpace->mesh->elementList)
    {
      this->reinit(elem);
      for (uint q = 0; q < FESpace_T::QR_T::numPts; ++q)
      {
        auto const localValue = this->evaluateDIV(q);
        auto const exactQ = exact(feSpace->curFE.qpoint[q]);
        error +=
            feSpace->curFE.JxW[q] * ((localValue - exactQ) * (localValue - exactQ));
      }
    }
    return error;
  }

  auto evaluate(uint const q) const
  {
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
    {
      FVec<FESpace_T::dim> dataQ = feSpace->curFE.phi[q].transpose() * dataLocal;
      // if constexpr (FESpace_T::dim == 1U)
      //   return dataQ[0]; // return a double directly
      return dataQ;
    }
    else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
    {
      Vec3 const dataQ = feSpace->curFE.phiVect[q].transpose() * dataLocal;
      return dataQ;
    }
    else
    {
      std::cerr << "only Lagrange and Raviart-Thomas implemented." << std::endl;
      std::exit(PROXPDE_NOT_IMPLEMENTED);
    }
  }

  auto evaluate0(uint const q) const
  {
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
    {
      Vec3 const dataQ = feSpace->curFE.phiVect0[q].transpose() * dataLocal;
      return dataQ;
    }
    else
    {
      std::cerr << "only Raviart-Thomas should require evaluate0()." << std::endl;
      std::exit(PROXPDE_NOT_IMPLEMENTED);
    }
  }

  auto evaluateGrad(uint const q) const
  {
    static_assert(
        family_v<RefFE_T> == FamilyType::LAGRANGE,
        "gradient is available only for Lagrange elements.");
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    return feSpace->curFE.dphi[q].transpose() * dataLocal;
  }

  auto evaluateDiv(uint const q) const
  {
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
    {
      auto const gradQ = evaluateGrad(q);
      return gradQ[0] + gradQ[4] + gradQ[8];
    }
    else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
    {
      double const dataQ = feSpace->curFE.divphi[q].transpose() * dataLocal;
      return dataQ;
    }
  }

  auto evaluateDIV(uint const q) const
  {
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
    {
      double const dataQ = feSpace->curFE.divphi0[q].transpose() * dataLocal;
      return dataQ;
    }
    else
    {
      std::cerr << "only Raviart-Thomas should require evaluateDIV()." << std::endl;
      std::exit(PROXPDE_NOT_IMPLEMENTED);
    }
  }

  // expensive version for points that are not the ones defined by the fixed qr rule
  auto evaluateOnRef(FVec<RefFE_T::dim> const & p) const
  {
    if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
    {
      Vec_T value = Vec_T::Zero();
      for (uint k = 0; k < FESpace_T::RefFE_T::numDOFs; ++k)
      {
        value += FESpace_T::RefFE_T::phiFun[k](p) * dataLocal.row(k);
      }
      return value;
    }
    else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
    {
      using QRTmp_T = DynamicQR<typename RefFE_T::GeoElem_T, 1U>;
      QRTmp_T::node = {p};
      VectorCurFE<RefFE_T, QRTmp_T> curFETmp;
      curFETmp.reinit(*(feSpace->curFE.elem));
      Vec3 const value = curFETmp.phiVect[0].transpose() * dataLocal;
      return value;
    }
    else
    {
      abort();
    }
  }

  // expensive version for points that are not the ones defined by the fixed qr rule
  auto evaluateGradOnRef(FVec<RefFE_T::dim> const & p) const
  {
    // available only for Lagrange fe types currently
    static_assert(family_v<RefFE_T> == FamilyType::LAGRANGE);
    using QRTmp_T = DynamicQR<typename RefFE_T::GeoElem_T, 1U>;
    QRTmp_T::node = {p};
    CurFE<RefFE_T, QRTmp_T> curFETmp;
    curFETmp.reinit(*(feSpace->curFE.elem));
    FMat<3, FESpace_T::dim> gradLocal = curFETmp.dphi[0].transpose() * dataLocal;
    return gradLocal;
  }

  // expensive version for points that are not the ones defined by the fixed qr rule
  auto evaluateDivOnRef(FVec<RefFE_T::dim> const & p) const
  {
    // available only for raviart-Thomas fe types currently
    static_assert(family_v<RefFE_T> == FamilyType::RAVIART_THOMAS);
    using QRTmp_T = DynamicQR<typename RefFE_T::GeoElem_T, 1U>;
    QRTmp_T::node = {p};
    VectorCurFE<RefFE_T, QRTmp_T> curFETmp;
    curFETmp.reinit(*(feSpace->curFE.elem));
    double const divOnPt = curFETmp.divphi[0].transpose() * dataLocal;
    return divOnPt;
  }

  // super-expensive version for points that are not the ones defined by the qr rule
  // and that we need to trace back to the ref element
  auto evaluateOnReal(Vec3 const & p) const
  {
    auto const pRef = feSpace->curFE.approxInverseMap(p);
    return evaluateOnRef(pRef);
  }

  // super-expensive version for points that are not the ones defined by the qr rule
  // and that we need to trace back to the ref element
  auto evaluateGradOnReal(Vec3 const & p) const
  {
    auto const pRef = feSpace->curFE.approxInverseMap(p);
    return evaluateGradOnRef(pRef);
  }

  // super-expensive version for points that are not the ones defined by the qr rule
  // and that we need to trace back to the ref element
  auto evaluateDivOnReal(Vec3 const & p) const
  {
    auto const pRef = feSpace->curFE.approxInverseMap(p);
    return evaluateDivOnRef(pRef);
  }

  template <typename FESpaceVec>
  void setData(Vec const & v, FESpaceVec const & feSpaceVec, uint component)
  {
    setComponent(data, feSpace, v, feSpaceVec, component);
  }

  template <typename FESpaceVec>
  void setFromComponent(Vec & v, FESpaceVec const & feSpaceVec, uint component)
  {
    getComponent(data, feSpace, v, feSpaceVec, component);
  }

  std::string name = "none";
  FESpace_T const * feSpace;
  Vec data;
  FMat<FESpace_T::RefFE_T::numDOFs, FESpace_T::dim> dataLocal;
};

// template <typename FEList>
// struct BlockFEVar
// {
//   BlockFEVar(std::string_view n, FEList & fe):
//     name(n),
//     feList(fe)
//   {
//     uint sum = 0;
//     static_for(feList, [&sum] (const auto /*i*/, const auto & feSpace){
//       sum += feSpace->dof.size;
//     });
//     data.resize(sum);
//   }
//
//   explicit BlockFEVar(FEList & fe):
//     feList{fe},
//   {
//     uint sum = 0;
//     static_for(feList, [&sum] (const auto /*i*/, const auto & feSpace){
//       sum += feSpace->dof.size;
//     });
//     data.resize(sum);
//   }
//
//   std::string const name = "";
//   FEList & feList;
//   Vec data;
// };

} // namespace proxpde
