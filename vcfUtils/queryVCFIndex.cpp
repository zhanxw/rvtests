#include <stdlib.h>
#include <algorithm>

#include "base/MmapFile.h"
#include "third/tabix/bgzf.h"

struct Record {
  int64_t pos;
  int64_t offset;
  bool operator<(Record& o) { return (pos < o.pos); }
};

bool comparator(const Record& a, const Record& b) { return a.pos < b.pos; }

int main(int argc, char** argv) {
  const char* fBG = argv[1];
  const char* fIndex = argv[2];
  int64_t pos = strtol(argv[3], NULL, 0);
  // const int Nrecord = 10;

  // read everything
  MmapFile mmapFile;
  mmapFile.open(fIndex);
  size_t Nrecord = mmapFile.getFileSize() / 16 - 1;
  Record* r = (Record*)mmapFile.data;

  // FILE* fp = fopen(fIndex, "rb");
  // if (Nrecord != fread(r, sizeof(Record), Nrecord, fp)) {
  //   fprintf(stderr, "Read error!\n");
  // }

  // binary search for file position
  int64_t offset = -1;
  Record query;
  query.pos = pos;
  // Comparator comparator;
  Record* lb = std::lower_bound(r, r + Nrecord, query,
                                comparator);  // r[lb].pos >= query.pos
  Record* ub = std::upper_bound(lb, r + Nrecord, query,
                                comparator);  // r[ub].pos > query.pos
  for (Record* pi = lb; pi != ub; ++pi) {
    printf("%ld %ld\n", pi->pos, pi->offset);
    offset = pi->offset;
    // (TODO) only store one virtual offset for now.
    break;
  }

  // int64_t offset = -1;
  // for (int i = 0; i < Nrecord; ++i) {
  //   if (r[i].pos == pos) {
  //     offset = r[i].offset;
  //     break;
  //   }
  // }
  if (offset < 0) {
    fprintf(stderr, "Cannot find position!\n");
  } else {
    printf("found: %ld %ld\n", pos, offset);
  }
  BGZF* fp2 = bgzf_open(fBG, "rb");
  if (bgzf_seek(fp2, offset, SEEK_SET)) {
    fprintf(stderr, "seek error!\n");
  }
  kstring_t* str;
  str = (kstring_t*)calloc(1, sizeof(kstring_t));
  kstring_t& s = *str;
  int ret = bgzf_getline(fp2, '\n', &s);
  if (ret <= 0) {
    fprintf(stderr, "getline error, ret = %d!\n", ret);
  }
  for (size_t i = 0; i < s.l; ++i) {
    if (i >= 50) break;
    printf("%c", s.s[i]);
  }
  printf("\n");

  free(str);
  bgzf_close(fp2);
  // fclose(fp);

  return 0;
}
