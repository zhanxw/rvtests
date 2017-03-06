#include "testArgumentAux.h"

#include "base/Argument.h"

DECLARE_INT_PARAMETER(intParam);

int test() {
  fprintf(stderr, "int parameter = %d\n", FLAG_intParam);
  return 0;
}
