#include "Profiler.h"

#include <vector>
#include <algorithm>

#include "SimpleTimer.h"

std::unordered_map<std::string, Profiler::Metric> Profiler::data;

void Profiler::addTimer(const char* func) {
  Metric& m = data[func];
  m.nHits++;
  m.timer.start();
}

void Profiler::deleteTimer(const char* func) {
  Metric& m = data[func];
  m.timer.stop();
}

struct FlatMetric{
  FlatMetric(const std::string& func, int nHits, double elapsed) {
    this->func = func;
    this->nHits = nHits;
    this->avgElapsed = elapsed / nHits;
    this->totalElapsed = elapsed;
  }
  std::string func;
  int nHits;
  double avgElapsed;
  double totalElapsed;
};

void Profiler::dump() {
  std::vector< FlatMetric > v;
  for (auto& x : Profiler::data) {
    v.push_back(FlatMetric(x.first, x.second.nHits, x.second.timer.getSeconds()));
  }
  std::sort(v.begin(), v.end(),
            [](const FlatMetric& a, const FlatMetric& b) -> bool {
              return a.avgElapsed > b.avgElapsed;
            });
  for (auto& x: v) {
    fprintf(
        stderr,
        "Function [ %s ] hit [ %d ] times, total elapsed time [ %g ] seconds, avg elapsed time [ %g ] seconds\n",
        x.func.c_str(), x.nHits, x.totalElapsed, x.avgElapsed);
  }
}
