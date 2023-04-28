#include "def.hpp"

#include <yaml-cpp/yaml.h>

int main(int /*argc*/, char * /*argv*/[])
{
  using namespace proxpde;

  // default config
  ParameterDict config;
  config["a"] = 2;
  config["b"] = 4.0;
  config["nested"]["c"] = 6;
  config["nested"]["d"] = 8;

  std::cout << "before override:\n" << config << std::endl;

  config.override("config.yaml");

  assert(config["a"].as<int>() == 3);
  assert(config["nested"]["c"].as<int>() == 5);
  assert(config["nested"]["d"].as<int>() == 8);

  std::cout << "after override:\n" << config << std::endl;

  return 0;
}
