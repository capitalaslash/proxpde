#pragma once

#include "def.hpp"

#include <chrono>

template <typename TimeUnit>
class Timer
{
public:
  typedef std::chrono::high_resolution_clock Clock_T;

  explicit Timer() {}

  void start()
  {
    _start = std::chrono::time_point_cast<TimeUnit>(Clock_T::now());
  }

  TimeUnit elapsed() const
  {
    return std::chrono::duration_cast<TimeUnit>(Clock_T::now() - _start);
  }

private:
  Clock_T::time_point _start;
};

template <typename TimeUnit>
std::ostream & operator<<(std::ostream & out, Timer<TimeUnit> const & timer)
{
  out << timer.elapsed().count();
  return out;
}

typedef Timer<std::chrono::milliseconds> MilliTimer;
