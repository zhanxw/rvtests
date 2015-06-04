#include <stdio.h>
#include <stdlib.h>

#include "Profiler.h"

void f(int n) {
  PROFILE_FUNCTION() ;
      
  double d = 0;
  for (int i = 0; i < n; ++i) {
    d += i;
  }
  return;
}

void g(int n) {
  PROFILE_FUNCTION() ;
      
  for (int i = 0; i < 10; ++i) {
    f(1000000);
  }
  return;
}

int main(int argc, char *argv[])
{

  g(10000);

  for (int i = 0 ; i < 15; ++ i) {
    PROFILE_SCOPE("loop i");
    PROFILE_NAME_START("loop 2 i");
    for (int j = 1; j < 5000000; ++j) {
      i *= j;
      i /= j;
    }
    PROFILE_NAME_STOP("loop 2 i");
  }
  
  PROFILE_DUMP();
  
  return 0;
}
