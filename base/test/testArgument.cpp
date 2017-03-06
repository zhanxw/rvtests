#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "base/Argument.h"
#include "testArgumentAux.h"

ADD_PARAMETER_GROUP("My parameter group");
ADD_INT_PARAMETER(intParam, "-i", "an integer parameter");

int main(int argc, char** argv) {
  PARSE_PARAMETER(argc, argv);

  fprintf(stderr,
          "\n------------------------- Help -------------------------\n");
  PARAMETER_HELP();
  fprintf(stderr,
          "\n------------------------- Status -------------------------\n");
  PARAMETER_STATUS();
  fprintf(stderr,
          "\n------------------------- Other -------------------------\n");

  // external argument

  test();

  return 0;
};
