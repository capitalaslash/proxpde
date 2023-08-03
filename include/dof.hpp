#pragma once

#include "def.hpp"

#include "mesh.hpp"
#include "reffe.hpp"

namespace proxpde
{

enum class DofOrdering : uint8_t
{
  BLOCK,
  INTERLEAVED,
};

enum class DofType : uint8_t
{
  CONTINUOUS,
  DISCONTINUOUS,
};

template <typename Mesh, typename RefFE, uint dimension, DofType type, DofOrdering ord>
struct DOF;

template <typename Mesh, typename RefFE, uint dimension, DofType type, DofOrdering ord>
std::ostream &
operator<<(std::ostream &, DOF<Mesh, RefFE, dimension, type, ord> const &);

template <typename Mesh, typename RefFE, uint dimension, DofType t, DofOrdering ord>
struct DOF
{
  static DofType constexpr type = t;
  static uint constexpr dim = dimension;
  static DofOrdering constexpr ordering = ord;
  static uint constexpr clms = numDOFs<RefFE>();
  using RefFE_T = RefFE;
  // using ElemMap_T = std::vector<std::array<DOFid_T,clms*dim>>;
  using ElemMap_T = Table<DOFid_T, clms * dim>;
  // using GeoMap_T = std::vector<std::array<DOFid_T,RefFE::numGeoDOFs>>;
  using GeoMap_T = Table<DOFid_T, RefFE::numGeoDOFs>;
  using PtMap_T = std::vector<DOFid_T>;
  // we need ordered sets in order to compare them avoiding permutation issues
  using edgeIdList_T = std::set<id_T>;
  using faceIdList_T = std::set<id_T>;

  DOF() = default;

  explicit DOF(Mesh const & mesh):
      rows{mesh.elementList.size()},
      elemMap(rows, clms * dim),
      ptMap(mesh.pointList.size(), dofIdNotSet),
      geoMap(rows, RefFE::numGeoDOFs)
  {
    setupElemMap(mesh);
    setupGeoMap(mesh);
  }

  void init(Mesh const & mesh)
  {
    rows = mesh.elementList.size();
    elemMap.resize(rows, clms * dim);
    ptMap.resize(mesh.pointList.size(), dofIdNotSet);
    geoMap.resize(rows, RefFE::numGeoDOFs);

    setupElemMap(mesh);
    setupGeoMap(mesh);
  }

  void setupElemMap(Mesh const & mesh);

  void setupGeoMap(Mesh const & mesh);

  DOFid_T getId(id_T elemId, id_T pos = 0, uint d = 0) const
  {
    return elemMap(elemId, pos + d * clms);
  }

  friend std::ostream & operator<< <>(std::ostream & out, DOF const & dof);

  std::tuple<id_T, id_T> findPos(DOFid_T const id) const
  {
    typename ElemMap_T::Index minRow, minCol;
    auto const min = (elemMap - ElemMap_T::Constant(rows, clms * dim, id))
                         .abs()
                         .minCoeff(&minRow, &minCol);
    assert(min == 0);
    return std::tie(minRow, minCol);
  }

  std::size_t rows;
  ElemMap_T elemMap;
  PtMap_T ptMap;
  GeoMap_T geoMap;
  // total number of DOFs for single component
  uint size;
  uint mapSize;
};

template <typename Mesh, typename RefFE, uint dim, DofType type, DofOrdering ordering>
inline std::ostream &
operator<<(std::ostream & out, DOF<Mesh, RefFE, dim, type, ordering> const & dof)
{
  out << "DOF map\n";
  out << "elemMap: " << dof.rows << "x" << dof.clms << "\n";
  out << dof.elemMap << "\n";
  out << "geoMap: " << dof.rows << "x" << RefFE::numGeoDOFs << "\n";
  out << dof.geoMap;
  return out;
}

template <typename Mesh, typename RefFE, uint dimension, DofType t, DofOrdering ord>
void DOF<Mesh, RefFE, dimension, t, ord>::setupElemMap(Mesh const & mesh)
{
  std::vector<DOFid_T> elemDOFs;
  if constexpr (RefFE::dofPlace[0])
  {
    elemDOFs.resize(mesh.elementList.size(), dofIdNotSet);
  }
  // mappings from a list of ids that identifies uniquely a face/edge and the dof
  // associated to it
  std::map<edgeIdList_T, DOFid_T> edgeDOFs;
  std::map<faceIdList_T, DOFid_T> faceDOFs;

  size = 0;
  for (auto const & e: mesh.elementList)
  {
    uint localDofCount = 0;

    // dofs on points
    if constexpr (RefFE::dofPlace[3])
    {
      for (auto & p: e.pts)
      {
        // check if dofs have already been assigned to this point
        // TODO: the DISCONTINUOUS check can be static
        if (ptMap[p->id] == dofIdNotSet || type == DofType::DISCONTINUOUS)
        {
          if constexpr (ordering == DofOrdering::BLOCK)
          {
            elemMap(e.id, localDofCount) = size;
          }
          else // ordering == DofOrdering::INTERLEAVED
          {
            for (uint d = 0; d < dim; ++d)
            {
              elemMap(e.id, localDofCount + d * clms) = size * dim + d;
            }
          }
          ptMap[p->id] = size;
          size++;
        }
        else
        {
          if constexpr (ordering == DofOrdering::BLOCK)
          {
            elemMap(e.id, localDofCount) = ptMap[p->id];
          }
          else // ordering == DofOrdering::INTERLEAVED
          {
            for (uint d = 0; d < dim; ++d)
            {
              elemMap(e.id, localDofCount + d * clms) = ptMap[p->id] * dim + d;
            }
          }
        }
        localDofCount++;
      }
    }

    // dofs on edges
    if constexpr (RefFE::dofPlace[2])
    {
      for (uint i = 0; i < Mesh::Elem_T::numEdges; i++)
      {
        edgeIdList_T edgeIDs;
        // edges can always be identified by 2 points
        for (uint j = 0; j < 2; j++)
        {
          edgeIDs.insert(e.pts[Mesh::Elem_T::elemToEdge[i][j]]->id);
        }

        // check if dofs have already been assigned to this edge
        if (!edgeDOFs.contains(edgeIDs) || type == DofType::DISCONTINUOUS)
        {
          if constexpr (ordering == DofOrdering::BLOCK)
          {
            elemMap(e.id, localDofCount) = size;
          }
          else // ordering == DofOrdering::INTERLEAVED
          {
            for (uint d = 0; d < dim; ++d)
            {
              elemMap(e.id, localDofCount + d * clms) = size * dim + d;
            }
          }
          edgeDOFs[edgeIDs] = size;
          size++;
        }
        else
        {
          if constexpr (ordering == DofOrdering::BLOCK)
          {
            elemMap(e.id, localDofCount) = edgeDOFs[edgeIDs];
          }
          else // ordering == DofOrdering::INTERLEAVED
          {
            for (uint d = 0; d < dim; ++d)
            {
              elemMap(e.id, localDofCount + d * clms) = edgeDOFs[edgeIDs] * dim + d;
            }
          }
        }
        localDofCount++;
      }
    }

    // dofs on faces
    if constexpr (RefFE::dofPlace[1])
    {
      for (uint i = 0; i < Mesh::Elem_T::numFaces; i++)
      {
        faceIdList_T faceIDs;
        for (uint j = 0; j < Mesh::Elem_T::Face_T::numPts; j++)
        {
          faceIDs.insert(e.pts[Mesh::Elem_T::elemToFace[i][j]]->id);
        }

        // check if dofs have already been assigned to this face
        if (!faceDOFs.contains(faceIDs) || type == DofType::DISCONTINUOUS)
        {
          if constexpr (ordering == DofOrdering::BLOCK)
          {
            elemMap(e.id, localDofCount) = size;
          }
          else // ordering == DofOrdering::INTERLEAVED
          {
            for (uint d = 0; d < dim; ++d)
            {
              elemMap(e.id, localDofCount + d * clms) = size * dim + d;
            }
          }
          faceDOFs[faceIDs] = size;
          size++;
        }
        else
        {
          if constexpr (ordering == DofOrdering::BLOCK)
          {
            elemMap(e.id, localDofCount) = faceDOFs[faceIDs];
          }
          else // ordering == DofOrdering::INTERLEAVED
          {
            for (uint d = 0; d < dim; ++d)
            {
              elemMap(e.id, localDofCount + d * clms) = faceDOFs[faceIDs] * dim + d;
            }
          }
        }
        localDofCount++;
      }
    }

    // dofs on element
    if constexpr (RefFE::dofPlace[0])
    {
      // check if dofs have already been assigned to this element
      if (elemDOFs[e.id] == idNotSet || type == DofType::DISCONTINUOUS)
      {
        if constexpr (ordering == DofOrdering::BLOCK)
        {
          elemMap(e.id, localDofCount) = size;
        }
        else // ordering == DofOrdering::INTERLEAVED
        {
          for (uint d = 0; d < dim; ++d)
          {
            elemMap(e.id, localDofCount + d * clms) = size * dim + d;
          }
        }
        elemDOFs[e.id] = size;
        size++;
      }
      else
      {
        if constexpr (ordering == DofOrdering::BLOCK)
        {
          elemMap(e.id, localDofCount) = elemDOFs[e.id];
        }
        else // ordering == DofOrdering::INTERLEAVED
        {
          for (uint d = 0; d < dim; ++d)
          {
            elemMap(e.id, localDofCount + d * clms) = elemDOFs[e.id] * dim + d;
          }
        }
      }
      localDofCount++;
    }
    assert(localDofCount == clms);
  }

  if constexpr (ordering == DofOrdering::BLOCK)
  {
    for (uint d = 0; d < dim - 1; d++)
    {
      elemMap.block(0, clms * (d + 1), rows, clms).setConstant(size * (d + 1));
      elemMap.block(0, clms * (d + 1), rows, clms) += elemMap.block(0, 0, rows, clms);
      // TODO: use eigen blocks
      // for(uint e=0; e<rows; ++e)
      // {
      //   for(uint i=0; i<clms; i++)
      //   {
      //      elemMap(e, clms*(d+1)+i) = size * (d+1) + elemMap(e, i);
      //   }
      // }
    }
  }
}

template <typename Mesh, typename RefFE, uint dimension, DofType t, DofOrdering ord>
void DOF<Mesh, RefFE, dimension, t, ord>::setupGeoMap(Mesh const & mesh)
{
  PtMap_T geoPtMap;
  if constexpr (mappingIsSeparate_v<RefFE>)
  {
    geoPtMap.resize(mesh.pointList.size(), dofIdNotSet);
  }

  mapSize = 0;
  if constexpr (mappingIsSeparate_v<RefFE>)
  {
    for (auto const & e: mesh.elementList)
    {
      // geometric mapping should always have points only
      static_assert(RefFE::geoPlace[3] == 1);
      {
        uint localMapCount = 0;
        for (auto & p: e.pts)
        {
          // check if dofs have already been assigned to this point
          if (geoPtMap[p->id] == dofIdNotSet)
          {
            geoMap(e.id, localMapCount) = mapSize;
            geoPtMap[p->id] = mapSize;
            mapSize++;
          }
          else
          {
            geoMap(e.id, localMapCount) = geoPtMap[p->id];
          }
          localMapCount++;
        }
      }
    }
  }

  if constexpr (!mappingIsSeparate_v<RefFE>)
  {
    if constexpr (ordering == DofOrdering::BLOCK)
    {
      geoMap = elemMap.block(0, 0, rows, RefFE::numGeoDOFs);
    }
    else
    {
      for (uint r = 0; r < geoMap.rows(); ++r)
      {
        for (uint c = 0; c < geoMap.cols(); ++c)
        {
          // TODO: build it in a more eccifient way
          geoMap(r, c) = this->getId(r, c) / dim;
        }
      }
    }
    mapSize = size;
  }
}

} // namespace proxpde
