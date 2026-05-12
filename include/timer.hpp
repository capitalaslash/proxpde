#pragma once

#include "def.hpp"

#include <chrono>

namespace proxpde
{

template <typename TimeUnit>
struct TimerTraits
{};

template <>
struct TimerTraits<std::chrono::milliseconds>
{
  static constexpr char const * uom = "ms";
};

template <typename TimeUnit>
class Timer
{
public:
  using Clock_T = std::chrono::high_resolution_clock;

  Timer() = default;

  void start(std::string_view const name = "no name")
  {
    _start = std::chrono::time_point_cast<TimeUnit>(Clock_T::now());
    _tmpName = name;
    if (!_data.contains(_tmpName))
    {
      _data[_tmpName].id = _currentId;
      _currentId++;
    }
  }

  void stop() { this->elapsed(); }

  TimeUnit elapsed()
  {
    auto const elapsed = std::chrono::duration_cast<TimeUnit>(Clock_T::now() - _start);
    _data[_tmpName].time += elapsed.count();
    return elapsed;
  }

  void print(FILE * out = stdout) const
  {
    uint totalTime = 0u;
    uint nameLength = 9u;
    // sort by id
    std::vector<std::pair<std::string, typename TimeUnit::rep>> sortedData(
        _data.size());
    for (auto const & [name, data]: _data)
    {
      nameLength = std::max(nameLength, static_cast<uint>(name.length()));
      totalTime += static_cast<uint>(data.time);
      sortedData[data.id] = {name, data.time};
    }

    fmt::print(out, "{}", Utils::separator);
    fmt::println(
        "| section{} | time ({}) |      % |",
        std::string(nameLength - 7u, ' '),
        TimerTraits<TimeUnit>::uom);
    fmt::print(out, "{}", Utils::separator);
    for (auto const & [name, time]: sortedData)
    {
      fmt::println(
          out,
          "| {:{}s} | {:9d} | {:6.2f} |",
          name,
          nameLength,
          time,
          100. * static_cast<double>(time) / totalTime);
    }
    fmt::print(out, "{}", Utils::separator);
    fmt::println(
        out,
        "| total{} | {:9d} |        |",
        std::string(nameLength - 5u, ' '),
        totalTime);
    fmt::print(out, "{}", Utils::separator);
  }

private:
  struct Data
  {
    int id;
    // int level;
    typename TimeUnit::rep time;
  };

  Clock_T::time_point _start;
  std::string _tmpName;
  std::unordered_map<std::string, Data> _data;
  int _currentId = 0;
  // int _currentLevel = 0;
};

template <typename TimeUnit>
std::ostream & operator<<(std::ostream & out, Timer<TimeUnit> & timer)
{
  out << timer.elapsed().count();
  return out;
}

using MilliTimer = Timer<std::chrono::milliseconds>;

} // namespace proxpde

template <>
struct fmt::formatter<proxpde::MilliTimer>: fmt::formatter<uint>
{
  auto format(proxpde::MilliTimer t, format_context & ctx) const
  {
    return formatter<uint>::format(static_cast<uint>(t.elapsed().count()), ctx);
  }
};
