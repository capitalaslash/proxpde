#pragma once

#include "def.hpp"
#include "xdmf_traits.hpp"

#include <fstream>
#include <tinyxml2.h>

template <typename FESpace>
struct IOManager
{
  typedef typename FESpace::Mesh_T::Elem_T Elem_T;
  typedef std::shared_ptr<typename FESpace::Mesh_T> MeshPtr_T;
  typedef XDMFTraits<Elem_T> Traits_T;

  explicit IOManager(FESpace const & space):
    feSpace(space)
  {}

  void print(Vec const& sol);

  FESpace const & feSpace;
};

// implementation --------------------------------------------------------------

template <typename FESpace>
void IOManager<FESpace>::print(Vec const& sol)
{
  typename FESpace::Mesh_T const & mesh = *(feSpace.meshPtr);
  uint const numPts = mesh.pointList.size();
  uint const numElems = mesh.elementList.size();

  //  system("mkdir -p output");
  const std::string fileName = "sol.xmf";

  using namespace tinyxml2;

  XMLDocument doc;
  doc.InsertEndChild(doc.NewDeclaration());
  doc.InsertEndChild(doc.NewUnknown("DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []"));

  auto xdmf_el = doc.NewElement("Xdmf");
  xdmf_el->SetAttribute("xmlns:xi", "http://www.w3.org/2003/XInclude");
  xdmf_el->SetAttribute("Version", "2.2");
  auto xdmf = doc.InsertEndChild(xdmf_el);

  auto domain = xdmf->InsertEndChild(doc.NewElement("Domain"));

  auto grid_el = doc.NewElement("Grid");
  grid_el->SetAttribute("GridType", "Uniform");
  auto grid = domain->InsertEndChild(grid_el);

  auto time_el = doc.NewElement("Time");
  time_el->SetAttribute("Value", 0.0);
  grid->InsertEndChild(time_el);

  auto topo_el = doc.NewElement("Topology");
  topo_el->SetAttribute("TopologyType", Traits_T::shape_name);
  topo_el->SetAttribute("Dimensions", numElems);
  auto topo = grid->InsertEndChild(topo_el);

  auto topodata_el = doc.NewElement("DataItem");
  topodata_el->SetAttribute("Dimensions",
    (std::to_string(numElems) + " " + std::to_string(Elem_T::numPts)).c_str());
  topodata_el->SetAttribute("NumberType", "Int");
  topodata_el->SetAttribute("Precision", 8);
  topodata_el->SetAttribute("Format", "XML");
  std::stringstream buf;
  buf << std::endl;
  for(auto& e: mesh.elementList)
  {
    for(auto n: e.pointList)
      buf << n->id << " ";
    buf << std::endl;
  }
  // buf << _data.name << ".mesh.h5:connectivity " << std::endl;
  topodata_el->SetText(buf.str().c_str());
  topo->InsertEndChild(topodata_el);

  auto geometry_el = doc.NewElement("Geometry");
  geometry_el->SetAttribute("GeometryType", "XYZ");
  auto geometry = grid->InsertEndChild(geometry_el);

  auto geodata_el = doc.NewElement("DataItem");
  geodata_el->SetAttribute("Dimensions",
    (std::to_string(numPts) + " 3").c_str());
  geodata_el->SetAttribute("NumberType", "Float");
  geodata_el->SetAttribute("Precision", 8);
  geodata_el->SetAttribute("Format", "XML");
  buf.str("");
  buf << std::endl;
  for(auto &p: mesh.pointList)
  {
    buf << p.coord[0] << " " << p.coord[1] << " " << p.coord[2] << std::endl;
  }
  // buf << _data.name << ".mesh.h5:coords " << std::endl;
  geodata_el->SetText(buf.str().c_str());
  geometry->InsertEndChild(geodata_el);

  for(auto& data: {sol})
  {
    Vec print_data(numPts);
    for(auto const & e: mesh.elementList)
    {
      uint count = 0;
      for(auto const & p: e.pointList)
      {
        print_data(p->id) = data(feSpace.dof.elemMap[e.id][count]);
        count++;
      }
    }
    std::string const varName = "sol";
    auto var_el = doc.NewElement("Attribute");
    var_el->SetAttribute("Name", varName.c_str());
    var_el->SetAttribute("Active", 1);
    var_el->SetAttribute("AttributeType", "Scalar");
    var_el->SetAttribute("Center", "Node");
    auto var = grid->InsertEndChild(var_el);

    auto vardata_el = doc.NewElement("DataItem");
    vardata_el->SetAttribute("Dimensions",
      (std::to_string(numPts) + " 1").c_str());
    vardata_el->SetAttribute("NumberType", "Float");
    vardata_el->SetAttribute("Precision", 8);
    vardata_el->SetAttribute("Format", "XML");
    buf.str("");
    buf << std::endl;
    for(uint p=0; p<numPts; ++p)
    {
      buf << print_data(p) << "\n";
    }
    // buf << _data.name << "." << step << ".h5:" << varName << std::endl;
    vardata_el->SetText(buf.str().c_str());
    var->InsertEndChild(vardata_el);
  }

  doc.SaveFile(fileName.c_str());
}
