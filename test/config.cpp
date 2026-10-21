#include <cassert> // for assert

#include <fmt/base.h>      // for println
#include <yaml-cpp/yaml.h> // for Node, operator<<

#include <proxpde/def.hpp> // for ParameterDict

int main(int /*argc*/, char * /*argv*/[])
{
  using namespace proxpde;

  // default config
  ParameterDict config;
  config["a"] = 2;
  config["b"] = 4.0;
  config["nested"]["c"] = 6;
  config["nested"]["d"] = 8;

  fmt::println("before override:\n{}", config);

  config.override("config.yaml");

  assert(config["a"].as<int>() == 3);
  assert(config["nested"]["c"].as<int>() == 5);
  assert(config["nested"]["d"].as<int>() == 8);

  fmt::println("after override:\n{}", config);

  return 0;
}
