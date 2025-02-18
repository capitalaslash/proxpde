#include "xdmf_doc.hpp"

namespace proxpde
{

XDMFDoc::XDMFDoc(std::filesystem::path const fp, XDMFGridType gridType): filePath{fp}
{
  assert(fp.extension() == ".xmf");

  auto decl = doc.prepend_child(pugi::node_declaration);
  decl.append_attribute("version") = "1.0";
  decl.append_attribute("encoding") = "UTF-8";

  doc.append_child(pugi::node_doctype).set_value("Xdmf SYSTEM \"Xdmf.dtd\" []");

  auto xdmfNode = doc.append_child("Xdmf");
  xdmfNode.append_attribute("xmlns:xi") = "http://www.w3.org/2003/XInclude";
  xdmfNode.append_attribute("Version") = "2.2";

  auto domainNode = xdmfNode.append_child("Domain");

  gridNode = domainNode.append_child("Grid");
  gridNode.append_attribute("GridType") = to_string(gridType).data();
  if (gridType == XDMFGridType::COLLECTION)
  {
    gridNode.append_attribute("CollectionType") = "Temporal";
  }
}

XDMFDoc::~XDMFDoc() { doc.save_file(filePath.c_str()); }

void XDMFDoc::setTime(double const time)
{
  auto timeNode = gridNode.append_child("Time");
  timeNode.append_attribute("Value") = time;
}

void XDMFDoc::setGeometry(
    std::filesystem::path meshPath,
    uint const mapSize,
    std::string_view const coordName)
{
  assert(meshPath.extension() == ".h5");
  auto geoNode = gridNode.append_child("Geometry");
  geoNode.append_attribute("GeometryType") = "XYZ";

  createDataItem(
      geoNode,
      {mapSize, 3},
      XDMFNumberType::FLOAT,
      8,
      XDMFFormat::HDF,
      meshPath.string() + ":/" + coordName.data());
}

void XDMFDoc::setVar(std::filesystem::path dataPath, XDMFVar const & var)
{
  assert(dataPath.extension() == ".h5");
  auto varNode = gridNode.append_child("Attribute");
  varNode.append_attribute("Name") = var.name.data();
  varNode.append_attribute("Active") = 1;
  varNode.append_attribute("AttributeType") = to_string(var.type).data();
  varNode.append_attribute("Center") = to_string(var.center).data();

  // auto const buf = filepath.filename().string() + ".time.h5:/" + name + "." + iter;
  createDataItem(
      varNode,
      {var.size, var.dim},
      var.numberType,
      8,
      XDMFFormat::HDF,
      dataPath.string() + ":/" + var.name);
}

void XDMFDoc::setTimeSeries(
    std::filesystem::path basePath,
    std::vector<std::pair<uint, double>> const & timeSeries)
{
  // TODO: validate that the timeSeries contains unique time values
  std::string timesStr = "";
  for (auto const & step: timeSeries)
  {
    timesStr += std::to_string(step.second) + " ";
    auto stepNode = gridNode.append_child("xi:include");
    std::filesystem::path p = basePath;
    p += fmt::format(".{}.xmf", step.first);
    stepNode.append_attribute("href") = p.c_str();
    stepNode.append_attribute("xpointer") = "element(/1/1/1)";
  }
  auto timeNode = gridNode.append_child("Time");
  timeNode.append_attribute("TimeType") = "List";
  createDataItem(
      timeNode,
      {1, timeSeries.size()},
      XDMFNumberType::FLOAT,
      8,
      XDMFFormat::INLINE,
      timesStr);
}

pugi::xml_node XDMFDoc::createDataItem(
    pugi::xml_node & parent,
    std::vector<ulong> const & dims,
    XDMFNumberType const type,
    uint const precision,
    XDMFFormat const format,
    std::string const & content)
{
  auto node = parent.append_child("DataItem");
  node.append_attribute("Dimensions") = join(dims).c_str();
  node.append_attribute("NumberType") = to_string(type).data();
  node.append_attribute("Precision") = precision;
  node.append_attribute("Format") = to_string(format).data();
  node.text() = ("\n" + content + "\n").data();
  return node;
}

} // namespace proxpde
