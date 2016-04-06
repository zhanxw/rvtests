#include "Indexer.h"

Indexer::Indexer(const std::vector<std::string>& a) {
  for (size_t i = 0; i < a.size(); ++i) {
    if (m.count(a[i]) == 0) {
      m[a[i]] = i;
    } else {
      duplication.insert(a[i]);
    }
  }
}

int Indexer::translate(const std::vector<std::string>& input,
                       std::vector<int>* out) {
  if (!out) {
    return -1;
  }

  out->clear();
  for (size_t i = 0; i != input.size(); ++i) {
    if (m.count(input[i])) {
      out->push_back(m[input[i]]);
    } else {
      out->push_back(-1);
    }
  }
  return 0;
}
