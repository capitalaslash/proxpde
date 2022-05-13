#pragma once

#include "def.hpp"

#include "fespace.hpp"

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

template <typename FESpace>
struct FEVar
{
  using FESpace_T = FESpace;
  using RefFE_T = typename FESpace_T::RefFE_T;
  using Vec_T = FVec<FESpace_T::physicalDim()>;

  FEVar(std::string_view n, FESpace_T const & fe):
      name(n),
      feSpace(fe),
      data(fe.dof.size * FESpace_T::dim)
  {}

  explicit FEVar(FESpace_T const & fe): feSpace(fe), data(fe.dof.size * FESpace_T::dim)
  {}

  FEVar<FESpace_T> & operator<<(Fun<FESpace_T::physicalDim(), 3> const & f)
  {
    interpolateAnalyticFunction(f, this->feSpace, this->data);
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

  FEVar<FESpace_T> & operator<<(double const & v)
  {
    static_assert(FESpace_T::dim == 1, "this is available only on scalar fe spaces");
    return operator<<([&v](Vec3 const &) { return Vec1::Constant(v); });
  }

  double operator[](id_T const id) const { return data[id]; }

  void reinit(GeoElem const & elem)
  {
    feSpace.curFE.reinit(elem);
    setLocalData(elem.id);
  }

  void setLocalData(id_T const elemId)
  {
    for (uint d = 0; d < FESpace_T::dim; ++d)
    {
      for (uint n = 0; n < FESpace_T::RefFE_T::numDOFs; ++n)
      {
        auto const id = feSpace.dof.getId(elemId, n, d);
        _localData(n, d) = data[id];
      }
    }
  }

  auto getFacetMeanValue(uint const side)
  {
    assert(side < FESpace_T::Mesh_T::Elem_T::numFacets);
    FVec<FESpace_T::dim> sum = FVec<FESpace_T::dim>::Zero();
    for (auto const dofFacet: RefFE_T::dofOnFacet[side])
    {
      sum += _localData.row(dofFacet).transpose();
    }
    sum /= RefFE_T::dofPerFacet;
    return sum;
  }

  auto integrate()
  {
    FVec<FESpace_T::dim> integral = FVec<FESpace_T::dim>::Zero();
    for (auto const & elem: feSpace.mesh->elementList)
    {
      this->reinit(elem);
      for (uint q = 0; q < FESpace_T::QR_T::numPts; ++q)
      {
        integral += feSpace.curFE.JxW[q] * this->evaluate(q);
      }
    }
    return integral;
  }

  auto evaluate(uint const q) const
  {
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    if constexpr (family_v<RefFE_T> == FamilyType::LAGRANGE)
    {
      FVec<FESpace_T::dim> dataQ = feSpace.curFE.phi[q].transpose() * _localData;
      return dataQ;
    }
    else if constexpr (family_v<RefFE_T> == FamilyType::RAVIART_THOMAS)
    {
      Vec3 const dataQ = feSpace.curFE.phiVect[q].transpose() * _localData;
      return dataQ;
    }
  }

  auto evaluateGrad(uint const q) const
  {
    static_assert(
        family_v<RefFE_T> == FamilyType::LAGRANGE,
        "gradient is available only for Lagrange elements.");
    // check that the qr is compatible
    assert(q < FESpace_T::QR_T::numPts);
    return feSpace.curFE.dphi[q].transpose() * _localData;
  }

  // expensive version for points that are not the ones defined by the qr rule
  auto evaluateOnRef(FVec<FESpace_T::RefFE_T::dim> const & p) const
  {
    FVec<FESpace_T::dim> value = FVec<FESpace_T::dim>::Zero();
    for (uint k = 0; k < FESpace_T::RefFE_T::numDOFs; ++k)
    {
      value += FESpace_T::RefFE_T::phiFun[k](p) * _localData.row(k);
    }
    return value;
  }

  // super-expensive version for points that are not the ones defined by the qr rule
  // and that we need to trace back to the ref element
  auto evaluateOnReal(Vec3 const & p) const
  {
    auto const pRef = feSpace.curFE.approxInverseMap(p);
    return evaluateOnRef(pRef);
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

  std::string const name = "";
  FESpace_T const & feSpace;
  Vec data;

private:
  FMat<FESpace_T::RefFE_T::numDOFs, FESpace_T::dim> _localData;
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
//       sum += feSpace.dof.size;
//     });
//     data.resize(sum);
//   }
//
//   explicit BlockFEVar(FEList & fe):
//     feList{fe},
//   {
//     uint sum = 0;
//     static_for(feList, [&sum] (const auto /*i*/, const auto & feSpace){
//       sum += feSpace.dof.size;
//     });
//     data.resize(sum);
//   }
//
//   std::string const name = "";
//   FEList & feList;
//   Vec data;
// };
