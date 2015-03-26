#ifndef _SIMPLETIMER_H_
#define _SIMPLETIMER_H_

#include <ctime>

class SimpleTimer {
 public:
  SimpleTimer() {
    _elapsed = 0.;
    start();
  };
  void start() { _start = time(0); }
  double stop() {
    _stop = time(0);
    _elapsed += (int)(_stop - _start);
    return _elapsed;
  }
  void restart() { _start = time(0); }
  double getSeconds() { return _elapsed; }

 private:
  time_t _start;
  time_t _stop;
  int _elapsed;
};

#ifdef USE_ACCURATE_TIMER

#include <ctime>
#include <ratio>
#include <chrono>

class AccurateTimer {
 public:
  AccurateTimer() { start(); };
  void start() { _start = std::chrono::high_resolution_clock::now(); }
  double stop() {
    _stop = std::chrono::high_resolution_clock::now();
    _elapsed +=
        0.001 *
        (std::chrono::duration_cast<std::chrono::milliseconds>(_stop - _start))
            .count();
    return _elapsed;
  }
  void restart() { _start = std::chrono::high_resolution_clock::now(); }
  double getSeconds() { return _elapsed; }

 private:
  std::chrono::high_resolution_clock::time_point _start;
  std::chrono::high_resolution_clock::time_point _stop;

  double _elapsed;
};
#endif

#endif /* _SIMPLETIMER_H_ */
