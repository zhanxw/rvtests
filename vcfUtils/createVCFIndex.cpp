#include <stdlib.h>
#include "third/tabix/bgzf.h"

int main(int argc, char** argv) {
  const char* fn = argv[1];
  BGZF* fp = bgzf_open(fn, "rb");
  kstring_t* str;
  str = (kstring_t*)calloc(1, sizeof(kstring_t));
  kstring_t& s = *str;
  FILE* fIndex = fopen("test.idx", "wb");
  int ret;
  int64_t offset;
  int64_t pos;
  do {
    offset = bgzf_tell(fp);
    ret = bgzf_getline(fp, '\n', &s);

    if (ret <= 0) {
      break;
    }
    size_t beg = 0;
    for (; beg < s.l; ++beg) {
      if (s.s[beg] == '\t') {
        pos = strtol(s.s + beg + 1, NULL, 0);
        break;
      }
    }
    if (s.s[0] == '#') {
      if (s.s[1] == '#') {
        continue;
      } else if (s.s[1] == '#') {  // header line
        pos = 0;
      } else {
        fprintf(stderr, "Strange header line!\n");
      }
    }
    // printf("%ld %ld\n", pos, offset);
    fwrite(&pos, sizeof(int64_t), 1, fIndex);
    fwrite(&offset, sizeof(int64_t), 1, fIndex);
  } while (1);

  bgzf_close(fp);
  fclose(fIndex);
  return 0;
}
