#include <vector>
#include <cassert>

#include "Indexer.h"

int main(int argc, char* argv[]) {
  {
    std::vector<std::string> v = {"a", "b", "c", "d", "a"};
    Indexer indexer(v);
    assert(indexer.hasDuplication());
    assert(indexer["a"] == 0);
    assert(indexer["d"] == 3);    
  }
  
  {
    std::vector<std::string> v = {"a", "b", "c"};
    Indexer indexer(v);
    assert(!indexer.hasDuplication());
    assert(indexer["a"] == 0);
    assert(indexer["d"] == -1);    
  }

  
  return 0;
}

