#ifndef _PROFILER_H_
#define _PROFILER_H_

#include <unordered_map>
#include <string>

#include "SimpleTimer.h"

class AccurateTimer;

class Profiler {
 public:
  Profiler();
  virtual ~Profiler();
  static void addTimer(const char* func);
  static void deleteTimer(const char* func);
  static void dump();

  struct Timer {
    Timer(const char* func) {
      this->func = func;
      Profiler::addTimer(func);
    }
    ~Timer() {
      Profiler::deleteTimer(func);
    }
    const char* func;
  };

  struct Metric {
    Metric() : nHits(0) {}
    int nHits;
    AccurateTimer timer;
  };

 private:
  static std::unordered_map<std::string, Metric> data;
};

#define PROFILE_DUMP() Profiler::dump();
#define CURRENT_FUNCTION() __PRETTY_FUNCTION__
#define PROFILE_FUNCTION() \
  Profiler::Timer tempProfilerTimer(CURRENT_FUNCTION());
#define PROFILE_SCOPE(text) \
  Profiler::Timer tempProfilerTimer##__LINE__(text);
#define PROFILE_NAME_START(text) \
  Profiler::addTimer(text)
#define PROFILE_NAME_STOP(text) \
  Profiler::deleteTimer(text)  
#endif /* _PROFILER_H_ */
