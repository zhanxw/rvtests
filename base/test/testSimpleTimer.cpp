#include "SimpleTimer.h"
#include <stdio.h>
#include <cassert>

int main(int argc, char *argv[])
{
  {
    SimpleTimer m;
    int a = 1;
    for (int i = 1; i < 10000; i++) {
      for (int j = 1; j < 10000; j++) {
        a *= i;
        a *= j;
        a /= j;
        a /= i;
      }
    }
    fprintf(stderr, "elapsed %g seconds, and a = %d\n", m.stop(), a);
  }
  {
    AccurateTimer m;
    int a = 1;
    for (int i = 1; i < 10000; i++) {
      for (int j = 1; j < 10000; j++) {
        a *= i;
        a *= j;
        a /= j;
        a /= i;
      }
    }
    fprintf(stderr, "elapsed %g seconds, and a = %d\n", m.stop(), a);
  }

  return 0;
}
