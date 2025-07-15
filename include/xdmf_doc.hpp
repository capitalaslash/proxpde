#pragma once

// defines
#include "def.hpp"

// stl
#include <filesystem>

// external libs
#include <pugixml.hpp>

// local
#include "xdmf_traits.hpp"

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
