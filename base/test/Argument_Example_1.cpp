#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#if 1
#include "Argument.h"

// BEGIN_PARAMETER_LIST();  // this line is optional
ADD_PARAMETER_GROUP("My parameter group");
ADD_BOOL_PARAMETER(isHelp, "-h", "provide help");
// END_PARAMETER_LIST();    // this line is optional
#endif

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

  if (FLAG_isHelp)
    fprintf(stdout, "-h is set\n");
  else
    fprintf(stdout, "-h is NOT set\n");

  return 0;
};
