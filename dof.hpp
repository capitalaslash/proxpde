#pragma once

#include "def.hpp"

template <typename Mesh, typename RefFE>
struct DOF
{
  explicit DOF(Mesh & mesh):
    _rows(mesh.elementList.size())
  {
    std::cout << "_clms = " << _clms << std::endl;
    _map.resize(_rows);
    uint dof_count = 0;
    for(auto & e: mesh.elementList)
    {
      uint local_dof_count = 0;

      // check if there are dofs on the points
      if(RefFE::dof_place[3])
      {
        for(auto & p: e.pointList)
        {
          // check if dofs have already been assigned to this point
          if(p->dof_ids.empty())
          {
            p->dof_ids.resize(RefFE::dof_place[3]);
            for(uint i=0; i<RefFE::dof_place[3]; i++)
            {
              p->dof_ids[i] = dof_count;
              _map[e.id][local_dof_count] = dof_count;
              dof_count++;
              local_dof_count++;
            }
          }
          else
          {
            for(uint i=0; i<RefFE::dof_place[3]; i++)
            {
              _map[e.id][local_dof_count] = p->dof_ids[i];
              local_dof_count++;
            }
          }
        }
      }

      // check if there dofs on the element
      if(RefFE::dof_place[0])
      {
        // check if dofs have already been assigned to this element
        if(e.dof_ids.empty())
        {
          e.dof_ids.resize(RefFE::dof_place[0]);
          for(uint i=0; i<RefFE::dof_place[0]; i++)
          {
            e.dof_ids[i] = dof_count;
            _map[e.id][local_dof_count] = dof_count;
            dof_count++;
            local_dof_count++;
          }
        }
        else
        {
          for(uint i=0; i<RefFE::dof_place[0]; i++)
          {
            _map[e.id][local_dof_count] = e.dof_ids[i];
            local_dof_count++;
          }
        }
      }
    }
    num = dof_count;
  }

  static uint constexpr _clms = numDOFs<RefFE>();
  uint _rows;
  std::vector<std::array<id_T,_clms>> _map;
  uint num;
};

template <typename Mesh, typename RefFE>
inline std::ostream & operator<<(std::ostream & out, DOF<Mesh, RefFE> const & dof)
{
  out << "DOF map with " << dof._rows << " and " << dof._clms << "\n";
  for(auto & row: dof._map)
  {
    for(auto & id: row)
      out << id << " ";
    out << "\n";
  }
  return out;
}
