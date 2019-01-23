#include <stdlib.h>
#include "third/tabix/bgzf.h"

int main(int argc, char** argv) {
  const char* fn = argv[1];
  int64_t a, b;
  FILE* fp = fopen(fn, "rb");
  while (1) {
    if (1 != fread(&a, sizeof(int64_t), 1, fp)) {
      break;
    }
    if (1 != fread(&b, sizeof(int64_t), 1, fp)) {
      break;
    }
    fprintf(stderr, "%ld\t%ld\n", a, b);
  }
  fclose(fp);
  return 0;
}
