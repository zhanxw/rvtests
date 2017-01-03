#ifndef _INDEXER_H_
#define _INDEXER_H_

#include <map>
#include <set>
#include <string>
#include <vector>

/**
 * Indexer class is convenient tool to index a vector class
 *
 * v = ["a", "b", "c"]
 * indexer will store {"a": 0,  "b": 1, "c": 2}
 *
 * NOTE: duplicated element is ignored by default
 * in another word, only the first unique elements are used.
 */
class Indexer {
 public:
  Indexer(const std::vector<std::string>& a);
  bool hasDuplication() { return this->duplication.size() > 0; }
  int operator[](const std::string& s) const {
    if (m.count(s) == 0) return -1;
    return m.find(s)->second;
  }
  int translate(const std::vector<std::string>& input, std::vector<int>* out);

 private:
  std::map<std::string, int> m;
  std::set<std::string> duplication;
};

#endif /* _INDEXER_H_ */
