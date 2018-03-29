#pragma once

#include "def.hpp"
#include "mesh.hpp"
#include "reffe.hpp"

#include <set>

template <typename Mesh, typename RefFE, uint dim>
struct DOF
{
  static uint constexpr clms = numDOFs<RefFE>();
  using RefFE_T = RefFE;
  using ElemMap_T = std::vector<array<DOFid_T,clms*dim>>;
  using GeoMap_T = std::vector<array<DOFid_T,RefFE::numGeoFuns>>;
  using PtMap_T = std::vector<DOFid_T>;
  // we need ordered sets in order to compare them avoiding permutation issues
  using edgeIdList_T = std::set<id_T>;
  using faceIdList_T = std::set<id_T>;

  explicit DOF(Mesh & mesh):
    rows(mesh.elementList.size())
  {
    elemMap.resize(rows);
    geoMap.resize(rows);
    PtMap_T ptMap(mesh.pointList.size(), DOFidNotSet);
    PtMap_T geoPtMap(mesh.pointList.size(), DOFidNotSet);
    std::vector<DOFid_T> elemDOFs(mesh.elementList.size(), DOFidNotSet);
    // mappings from a list of ids that identifies uniquely a face/edge and the dof associated to it
    std::map<edgeIdList_T,DOFid_T> edgeDOFs;
    std::map<faceIdList_T,DOFid_T> faceDOFs;

    mapSize = 0;
    size = 0;
    for(auto const & e: mesh.elementList)
    {
      uint localDofCount = 0;

      // dofs on points
      if constexpr (RefFE::dofPlace[3])
      {
        for(auto & p: e.pointList)
        {
          // check if dofs have already been assigned to this point
          if(ptMap[p->id] == DOFidNotSet)
          {
            elemMap[e.id][localDofCount] = size;
            ptMap[p->id] = size;
            size++;
            localDofCount++;
          }
          else
          {
            elemMap[e.id][localDofCount] = ptMap[p->id];
            localDofCount++;
          }
        }
      }

      // dofs on edges
      if constexpr (RefFE::dofPlace[2])
      {
        for(uint i=0; i<Mesh::Elem_T::numEdges; i++)
        {
          edgeIdList_T edgeIDs;
          // edges can always be identified by 2 points
          for(uint j=0; j<2; j++)
          {
            edgeIDs.insert(e.pointList[Mesh::Elem_T::elemToEdge[i][j]]->id);
          }

          // check if dofs have already been assigned to this edge
          if (edgeDOFs.find(edgeIDs) == edgeDOFs.end())
          {
            edgeDOFs[edgeIDs] = size;
            elemMap[e.id][localDofCount] = size;
            size++;
            localDofCount++;
          }
          else
          {
            elemMap[e.id][localDofCount] = edgeDOFs[edgeIDs];
            localDofCount++;
          }
        }
      }

      // dofs on faces
      if constexpr (RefFE::dofPlace[1])
      {
        for(uint i=0; i<Mesh::Elem_T::numFaces; i++)
        {
          faceIdList_T faceIDs;
          for(uint j=0; j<Mesh::Elem_T::Face_T::numPts; j++)
          {
            faceIDs.insert(e.pointList[Mesh::Elem_T::elemToFace[i][j]]->id);
          }

          // check if dofs have already been assigned to this face
          if(faceDOFs.find(faceIDs) == faceDOFs.end())
          {
            faceDOFs[faceIDs] = size;
            elemMap[e.id][localDofCount] = size;
            size++;
            localDofCount++;
          }
          else
          {
            elemMap[e.id][localDofCount] = faceDOFs[faceIDs];
            localDofCount++;
          }
        }
      }

      // dofs on element
      if constexpr (RefFE::dofPlace[0])
      {
        // check if dofs have already been assigned to this element
        if(elemDOFs[e.id] == DOFidNotSet)
        {
          elemDOFs[e.id] = size;
          elemMap[e.id][localDofCount] = size;
          size++;
          localDofCount++;
        }
        else
        {
          elemMap[e.id][localDofCount] = elemDOFs[e.id];
          localDofCount++;
        }
      }
      assert(localDofCount == clms);

      if constexpr (MappingIsSeparate<RefFE>::value)
      {
        // geometric mapping has points
        if constexpr (RefFE::geoPlace[3])
        {
          uint localMapCount = 0;
          for(auto & p: e.pointList)
          {
            // check if dofs have already been assigned to this point
            if(geoPtMap[p->id] == DOFidNotSet)
            {
              geoMap[e.id][localMapCount] = mapSize;
              geoPtMap[p->id] = mapSize;
              mapSize++;
            }
            else
            {
              geoMap[e.id][localMapCount] = geoPtMap[p->id];
            }
            localMapCount++;
          }
        }
      }
    }

    for(uint d=0; d<dim-1; d++)
    {
      for(uint e=0; e<rows; ++e)
      {
        for(uint i=0; i<clms; i++)
        {
          elemMap[e][clms*(d+1)+i] = size * (d+1) + elemMap[e][i];
        }
      }
    }

    if constexpr (!MappingIsSeparate<RefFE>::value)
    {
      // check sizes
      assert(geoMap.size() == elemMap.size() && geoMap[0].size() == RefFE::numGeoFuns);
      // no point in building the same matrix twice
      // geoMap = elemMap;
      for (uint r=0; r<geoMap.size(); ++r)
      {
        for (uint c=0; c<RefFE::numGeoFuns; ++c)
        {
          geoMap[r][c] = elemMap[r][c];
        }
      }
      mapSize = size;
    }
  }

  uint rows;
  ElemMap_T elemMap;
  GeoMap_T geoMap;
  // total number of DOFs for single component
  uint size;
  uint mapSize;
};

template <typename Mesh, typename RefFE, uint dim>
inline std::ostream & operator<<(std::ostream & out, DOF<Mesh,RefFE,dim> const & dof)
{
  out << "DOF map\n";
  out << "elemMap: " << dof.rows << "x" << dof.clms << "\n";
  for(auto & row: dof.elemMap)
  {
    for(auto & id: row)
      out << id << " ";
    out << "\n";
  }
  out << "geoMap: " << dof.rows << "x" << RefFE::numGeoFuns << "\n";
  for(auto & row: dof.geoMap)
  {
    for(auto & id: row)
      out << id << " ";
    out << "\n";
  }
  return out;
}
