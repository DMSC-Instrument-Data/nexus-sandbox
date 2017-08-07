#ifndef TIMER_H
#define TIMER_H

class Timer {
public:
  /// Create timer and store initial time stamp.
  Timer() { checkpoint(); }

  /// Store timestap and return elapsed time since initial time stamp.
  double checkpoint() {
    m_times.push_back(std::chrono::system_clock::now());
    return seconds(m_times.back() - m_times.front());
  }

  /// Return vector of deltas between all stored time stamps.
  std::vector<double> deltas() const {
    std::vector<double> result;
    for (size_t i = 1; i < m_times.size(); ++i)
      result.push_back(seconds(m_times[i] - m_times[i - 1]));
    return result;
  }

private:
  double seconds(std::chrono::duration<double> elapsed) const {
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
               .count() /
           (1000000000);
  }

  std::vector<std::chrono::system_clock::time_point> m_times;
};

#endif // TIMER_H
