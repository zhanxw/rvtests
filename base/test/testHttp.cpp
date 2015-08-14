#include "Http.h"
#include <cassert>

int main() {
  std::vector<std::string> res;

  {
    Http http("www.google.com", 80);
    assert(http.read(&res) >= 0);
    assert(res.size() > 0);
    fprintf(stderr, "Google.com responded %zu lines\n", res.size());
  }

  {
    Http http("http://zhanxw.com/rvtests/version", 80);
    assert(http.read(&res) >= 0);
    assert(res.size() > 0);
    for(size_t i = 0; i < res.size(); ++i) {
      printf("'%s'\n", res[i].c_str());
    }
  }
  
  {
    Http s("zhanxw.com", 12345);
    assert(s.read(&res) < 0);
    
  }

  {
    Http s("unreachable.xxx", 12345);
    assert(s.read(&res) != 0);
  }

  return 0;
}
