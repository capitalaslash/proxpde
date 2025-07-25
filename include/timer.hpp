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

  void print(std::ostream & out = std::cout) const
  {
    int totalTime = 0;
    int strLength = 9;
    // sort by id
    std::vector<std::pair<std::string, typename TimeUnit::rep>> sortedData(
        _data.size());
    for (auto const & [name, data]: _data)
    {
      strLength = std::max(strLength, static_cast<int>(name.length()));
      totalTime += static_cast<int>(data.time);
      sortedData[data.id] = {name, data.time};
    }

    auto const curSettings = out.flags();
    auto const curPrecision = out.precision();
    out << Utils::separator << "| " << std::setw(strLength) << "section"
        << " | time (" << TimerTraits<TimeUnit>::uom << ") |     % |\n"
        << Utils::separator;
    for (auto const & [name, time]: sortedData)
    {
      out << "| " << std::setw(strLength) << name << " | " << std::setw(9) << time
          << " | " << std::fixed << std::setprecision(2) << std::setw(5)
          << 100. * static_cast<double>(time) / totalTime << " |\n";
    }
    out << Utils::separator;
    out << "| " << std::setw(static_cast<int>(strLength)) << "total"
        << " | " << std::setw(9) << totalTime << " |       |\n";
    out << Utils::separator;
    out.flags(curSettings);
    out.precision(curPrecision);
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
