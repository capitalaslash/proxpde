#pragma once

#include "def.hpp"

#include <chrono>

template <typename TimeUnit>
struct TimerTraits {};

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
    _tmpName = name;
    _start = std::chrono::time_point_cast<TimeUnit>(Clock_T::now());
  }

  void stop()
  {
    this->elapsed();
  }

  TimeUnit elapsed()
  {
    auto const elapsed = std::chrono::duration_cast<TimeUnit>(Clock_T::now() - _start);
    _times.push_back(std::pair(_tmpName, elapsed.count()));
    return elapsed;
  }

  void print(std::ostream & out = std::cout) const
  {
    double totalTime = 0;
    size_t strLength = 9;
    for (auto const [name, time]: _times)
    {
      strLength = std::max(strLength, name.length());
      totalTime += time;
    }

    auto const curSettings = out.flags();
    auto const curPrecision = out.precision();
    out << separator << "| "  << std::string(strLength-8, ' ') << " section | time (" << TimerTraits<TimeUnit>::uom << ") |     % |\n" << separator;
    for (auto const [name, time]: _times)
    {
      out << "| " << std::setw(strLength) << name
          << " | " << std::setw(9) << time
          << " | " << std::fixed << std::setprecision(2) << std::setw(5) << 100. * time / totalTime << " |" << std::endl;
    }
    out << separator;
    out.flags(curSettings);
    out.precision(curPrecision);
  }

private:
  Clock_T::time_point _start;
  std::string _tmpName;
  std::list<std::pair<std::string, typename TimeUnit::rep>> _times;
};

template <typename TimeUnit>
std::ostream & operator<<(std::ostream & out, Timer<TimeUnit> & timer)
{
  out << timer.elapsed().count();
  return out;
}

using MilliTimer = Timer<std::chrono::milliseconds>;
