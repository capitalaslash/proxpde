#pragma once

#include "def.hpp"
#include "reffe.hpp"

#include <set>

template <typename Mesh, typename RefFE, uint dim>
struct DOF
{
  static uint constexpr _clms = numDOFs<RefFE>();
  using ElemMap_T = std::vector<array<DOFid_T,_clms*dim>>;
  using PtMap_T = std::vector<DOFid_T>;

  explicit DOF(Mesh & mesh):
    _rows(mesh.elementList.size())
  {
    elemMap.resize(_rows);
    ptMap.resize(mesh.pointList.size(), DOFidNotSet);
    std::vector<DOFid_T> elemDOFs(mesh.elementList.size(), DOFidNotSet);
    std::map<std::set<id_T>,DOFid_T> edgeDOFs;
    std::map<std::set<id_T>,DOFid_T> faceDOFs;

    uint dof_count = 0;
    for(auto const & e: mesh.elementList)
    {
      uint local_dof_count = 0;

      // check if there are dofs on the points
      if(RefFE::dof_place[3])
      {
        for(auto & p: e.pointList)
        {
          // check if dofs have already been assigned to this point
          if(ptMap[p->id] == DOFidNotSet)
          {
            elemMap[e.id][local_dof_count] = dof_count;
            ptMap[p->id] = dof_count;
            dof_count++;
            local_dof_count++;
          }
          else
          {
            elemMap[e.id][local_dof_count] = ptMap[p->id];
            local_dof_count++;
          }
        }
      }

      // check if there dofs on the edges
      if(RefFE::dof_place[2])
      {
        for(uint i=0; i<Mesh::Elem_T::numEdges; i++)
        {
          std::set<id_T> edgeIDs;
          for(uint j=0; j<Mesh::Elem_T::Edge_T::numPts; j++)
          {
            edgeIDs.insert(e.pointList[Mesh::Elem_T::elemToEdge[i][j]]->id);
          }

          // check if dofs have already been assigned to this edge
          auto it = edgeDOFs.find(edgeIDs);
          if(it == edgeDOFs.end())
          {
            edgeDOFs[edgeIDs] = dof_count;
            elemMap[e.id][local_dof_count] = dof_count;
            dof_count++;
            local_dof_count++;
          }
          else
          {
            elemMap[e.id][local_dof_count] = edgeDOFs[edgeIDs];
            local_dof_count++;
          }
        }
      }

      // check if there dofs on the faces
      if(RefFE::dof_place[1])
      {
        for(uint i=0; i<Mesh::Elem_T::numFaces; i++)
        {
          std::set<id_T> faceIDs;
          for(uint j=0; j<Mesh::Elem_T::Face_T::numPts; j++)
          {
            faceIDs.insert(e.pointList[Mesh::Elem_T::elemToFace[i][j]]->id);
          }

          // check if dofs have already been assigned to this face
          auto it = faceDOFs.find(faceIDs);
          if(it == faceDOFs.end())
          {
            faceDOFs[faceIDs] = dof_count;
            elemMap[e.id][local_dof_count] = dof_count;
            dof_count++;
            local_dof_count++;
          }
          else
          {
            elemMap[e.id][local_dof_count] = faceDOFs[faceIDs];
            local_dof_count++;
          }
        }
      }

      // check if there dofs on the element
      if(RefFE::dof_place[0])
      {
        // check if dofs have already been assigned to this element
        if(elemDOFs[e.id] == DOFidNotSet)
        {
          elemDOFs[e.id] = dof_count;
          elemMap[e.id][local_dof_count] = dof_count;
          dof_count++;
          local_dof_count++;
        }
        else
        {
          elemMap[e.id][local_dof_count] = elemDOFs[e.id];
          local_dof_count++;
        }
      }
      assert(local_dof_count == _clms);
    }
    totalNum = dof_count;

    for(uint d=0; d<dim-1; d++)
    {
      for(uint e=0; e<_rows; ++e)
      {
        for(uint i=0; i<_clms; i++)
        {
          elemMap[e][_clms*(d+1)+i] = totalNum * (d+1) + elemMap[e][i];
        }
      }
    }
  }

  uint _rows;
  ElemMap_T elemMap;
  PtMap_T ptMap;
  // total number of DOFs for single component
  uint totalNum;
};

template <typename Mesh, typename RefFE, uint dim>
inline std::ostream & operator<<(std::ostream & out, DOF<Mesh,RefFE,dim> const & dof)
{
  out << "DOF map\n";
  out << "elemMap: " << dof._rows << "x" << dof._clms << "\n";
  for(auto & row: dof.elemMap)
  {
    for(auto & id: row)
      out << id << " ";
    out << "\n";
  }
  out << "ptMap: " << dof.ptMap.size() << "\n";
  for(auto & i: dof.ptMap)
  {
    out << i << " ";
  }
  out << "\n";
  return out;
}
