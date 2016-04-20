#include <vector>
#include <cassert>

#include "CommonFunction.h"

void assertSame(const std::vector<double>& in, double* out) {
  for (size_t i = 0; i != in.size(); ++i) {
    if (in[i] != out[i]) {
      fprintf(stderr, "rank does not look good, %g != %g\n", in[i], out[i]);
      assert(false);
    };
  }
}

int main(int argc, char* argv[]) {
  {
    // Test rank ////////////////////////////////////////////////
    std::vector<double> rank;
    {
      int a[] = {1, 2, 3};
      double correct[] = {0, 1, 2};
      std::vector<double> in(a, a + sizeof(a) / sizeof(a[0]));
      calculateRank(in, &rank);
      assertSame(rank, correct);
    }

    {
      int a[] = {1, 1, 2};
      double correct[] = {0.5, 0.5, 2};
      std::vector<double> in(a, a + sizeof(a) / sizeof(a[0]));
      calculateRank(in, &rank);
      assertSame(rank, correct);
    }

    {
      int a[] = {1, 1, 1, 1};
      double correct[] = {1.5, 1.5, 1.5, 1.5};
      std::vector<double> in(a, a + sizeof(a) / sizeof(a[0]));
      calculateRank(in, &rank);
      assertSame(rank, correct);
    }

    {
      int a[] = {3, 1, 2, 4};
      double correct[] = {2, 0, 1, 3};
      std::vector<double> in(a, a + sizeof(a) / sizeof(a[0]));
      calculateRank(in, &rank);
      assertSame(rank, correct);
    }
  }

  {
    std::vector<int> a = {1, 2, 3};
    std::vector<int> b = {2, 4};
    remove(b, &a);
    assert(a.size() == 2);
    assert(a[0] == 1);
    assert(a[1] == 3);
  }

  {
    std::vector<int> a = {1, 2, 3};
    std::vector<int> b = {2, 4};
    remove(a, &b);
    assert(b.size() == 1);
    assert(b[0] == 4);
  }
  
  return 0;
}
