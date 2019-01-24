#include <cassert>
#include "Http.h"

int main() {
  HttpResponse res;

  {
    Http http("http://www.google.com");
    assert(http.read(&res) >= 0);
    assert(res.size() > 0);
    fprintf(stderr, "Google.com responded %zu lines\n", res.size());
  }

  {
    Http http("http://zhanxw.com:80/rvtests/version");
    assert(http.read(&res) >= 0);
    assert(res.size() > 0);
    for (size_t i = 0; i < res.size(); ++i) {
      printf("'%s'\n", res.getBody()[i].c_str());
    }
  }

  {
    Http s("zhanxw.com");
    assert(s.read(&res) < 0);
  }

  {
    Http s("http://unreachable.xxx");
    assert(s.read(&res) != 0);
  }

  return 0;
}
