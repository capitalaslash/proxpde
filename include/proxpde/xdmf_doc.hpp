#pragma once

#include <cassert>     // for assert
#include <cstdlib>     // for size_t, uint
#include <filesystem>  // for path, operator==
#include <string>      // for allocator, basic_string, char_traits
#include <string_view> // for string_view
#include <utility>     // for pair
#include <vector>      // for vector

#include <pugixml.hpp> // for xml_node, xml_attribute, xml_docu...
#include <sys/types.h> // for ulong

#include <proxpde/xdmf_traits.hpp> // for XDMFGridType, XDMFFormat, XDMFNum...

namespace proxpde
{

template <typename T>
std::string join(std::vector<T> const & v, std::string_view const divider = " ")
{
  std::string buf;
  for (auto const & i: v)
  {
    buf += std::to_string(i);
    buf += divider;
  }
  return buf.substr(0, buf.size() - 1);
}

class XDMFDoc
{
public:
  explicit XDMFDoc(
      std::filesystem::path const fp, XDMFGridType gridType = XDMFGridType::SINGLE);

  ~XDMFDoc();

  void setTime(double const time);

  template <typename RefFE>
  void setTopology(
      std::filesystem::path meshPath,
      size_t const numElems,
      std::string_view const connName = "connectivity");

  void setGeometry(
      std::filesystem::path meshPath,
      size_t const mapSize,
      std::string_view const coordName = "coords");

  void setVar(std::filesystem::path dataPath, XDMFVar const & var);

  void setTimeSeries(
      std::filesystem::path basePath,
      std::vector<std::pair<uint, double>> const & timeSeries);

private:
  pugi::xml_node createDataItem(
      pugi::xml_node & parent,
      std::vector<ulong> const & dims,
      XDMFNumberType const type,
      uint const precision,
      XDMFFormat const format,
      std::string const & content);

  std::filesystem::path filePath;
  pugi::xml_document doc;
  pugi::xml_node gridNode;
};

template <typename RefFE>
void XDMFDoc::setTopology(
    std::filesystem::path meshPath,
    size_t const numElems,
    std::string_view const connName)
{
  assert(meshPath.extension() == ".h5");
  auto topoNode = gridNode.append_child("Topology");
  topoNode.append_attribute("TopologyType") = XDMFTraits<RefFE>::shapeName;
  topoNode.append_attribute("Dimensions") = numElems;

  createDataItem(
      topoNode,
      {numElems, RefFE::numGeoDOFs},
      XDMFNumberType::INT,
      8,
      XDMFFormat::HDF,
      meshPath.string() + ":/" + connName.data());
}

} // namespace proxpde
