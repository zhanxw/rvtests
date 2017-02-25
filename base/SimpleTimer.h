#ifndef _SIMPLETIMER_H_
#define _SIMPLETIMER_H_

#include <ctime>

/**
 * A simple timer
 *
 * As time() returns the seconds since the Epoch,
 * the precision here is second
 *
 */
class SimpleTimer {
 public:
  SimpleTimer() { start(); };
  void start() {
    _start = time(0);
    _stop = _start;
    _elapsed = 0;
  }
  double stop() {
    _stop = time(0);
    _elapsed = (int)(_stop - _start);
    return _elapsed;
  }
  double getSeconds() const { return _elapsed; }

 private:
  time_t _start;
  time_t _stop;
  int _elapsed;
};

//////////////////////////////////////////////////
// Accurate timer is platform dependent and also language (C++11) dependent
// We use a series of preprocessors below.
// make sure at most one of the following preprocessor is set
// _USE_CXX11, _WIN32, __linux, __APPLE__

#if !defined(_UES_CXX11) && !defined(_WIN32) && !defined(__linux) && \
    !defined __APPLE__
#pragma message "Use crude timer"
class AccurateTimer : public SimpleTimer {};
#endif

#ifdef _USE_CXX11
//////////////////////////////////////////////////
// C++11 standard version of accurate timer
#include <chrono>
#include <ctime>
#include <ratio>

class AccurateTimer {
 public:
  AccurateTimer() { start(); }
  void start() {
    _start = std::chrono::high_resolution_clock::now();
    _elapsed = 0.;
  }
  double stop() {
    _stop = std::chrono::high_resolution_clock::now();
    _elapsed =
        0.001 *
        (std::chrono::duration_cast<std::chrono::milliseconds>(_stop - _start))
            .count();
    return _elapsed;
  }
  double getSeconds() const { return _elapsed; }

 private:
  std::chrono::high_resolution_clock::time_point _start;
  std::chrono::high_resolution_clock::time_point _stop;
  double _elapsed;
};
#endif

#if !defined(_USE_CXX11) && defined(__APPLE__)
#pragma message "Use native Apple Timer"
#include <mach/clock.h>
#include <mach/mach.h>
#include <sys/time.h>
#include <time.h>

class AccurateTimer {
 public:
  AccurateTimer() { start(); }
  void start() {
    clock_gettime(&_start);
    _elapsed = 0.;
  }
  double stop() {
    clock_gettime(&_end);
    _elapsed = (_end.tv_sec - _start.tv_sec) +
               1.0e-9 * (_end.tv_nsec - _start.tv_nsec);
    return _elapsed;
  }
  double getSeconds() const { return _elapsed; }

 private:
  void clock_gettime(struct timespec* ts) {
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts->tv_sec = mts.tv_sec;
    ts->tv_nsec = mts.tv_nsec;
  }

 private:
  struct timespec _start;
  struct timespec _end;
  double _elapsed;
};

#endif

#if !defined(_USE_CXX11) && defined(__linux)
//////////////////////////////////////////////////
// Linux/BSD high accuracy timer
#include <time.h>  // clock_gettime(), need to link with -lrt

class AccurateTimer {
 public:
  AccurateTimer() { start(); }
  void start() {
    clock_gettime(CLOCK_MONOTONIC, &_start);
    _end = _start;
    _elapsed = 0.;
  }
  double stop() {
    clock_gettime(CLOCK_MONOTONIC, &_end);
    _elapsed = (_end.tv_sec - _start.tv_sec) +
               1.0e-9 * (_end.tv_nsec - _start.tv_nsec);
    return _elapsed;
  }
  double getSeconds() const { return _elapsed; }

 private:
  struct timespec _start;
  struct timespec _end;

  double _elapsed;
};
#endif

#if !defined(_USE_CXX11) && defined(_WIN32)
//////////////////////////////////////////////////
// Windows version of accurate timer

// from:
// https://msdn.microsoft.com/en-us/library/windows/desktop/dn553408(v=vs.85).aspx

#include <windows.h>

class AccurateTimer {
 public:
  AccurateTimer() { start(); }
  void start() {
    QueryPerformanceFrequency(&_freq);
    QueryPerformanceCounter(&_start);
    _elapsed = 0.;
  }
  double stop() {
    QueryPerformanceCounter(&_end);
    _elapsed += (_end.QuadPart - _start.QuadPart);
    _elapsed /= _freq.QuadPart;
    return _elapsed;
  }
  double getSeconds() const { return _elapsed; }

 private:
  LARGE_INTEGER _start;
  LARGE_INTEGER _end;
  LARGE_INTEGER _freq;

  double _elapsed;
};

#endif  // end #ifdef _WIN32

#endif /* _SIMPLETIMER_H_ */
