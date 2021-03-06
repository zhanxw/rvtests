#include <stdio.h>
#include <stdlib.h>

#include "base/IO.h"

uint64_t offset = 0;

void read(BufferedReader* br, char* ret, int len) {
  if (len != br->read(ret, len)) {
    fprintf(stderr, "Does not read enough bytes!\n");
  }
  offset += 4;
}
void read(BufferedReader* br, int32_t* ret) {
  char buf[4];
  int len = 4;
  if (len != br->read(buf, len)) {
    fprintf(stderr, "Does not read enough bytes!\n");
  }
  (*ret) = *(int32_t*)buf;
  offset += 4;
}
void read(BufferedReader* br, uint32_t* ret) {
  char buf[4];
  int len = 4;
  if (len != br->read(buf, len)) {
    fprintf(stderr, "Does not read enough bytes!\n");
  }
  (*ret) = *(uint32_t*)buf;
  offset += 4;
}
void read(BufferedReader* br, uint64_t* ret) {
  char buf[8];
  int len = 8;
  if (len != br->read(buf, len)) {
    fprintf(stderr, "Does not read enough bytes!\n");
  }
  (*ret) = *(uint64_t*)buf;
  offset += 4;
}

void printOffset(int len) { printf("[%8zu] ", offset + len); }

void dump(const char* name, const char* ret, int len) {
  printOffset(-len);
  printf("%s = ", name);
  for (int i = 0; i < len; ++i) {
    printf("%c", ret[i]);
  }
  printf(" ");
  for (int i = 0; i < len; ++i) {
    if (i) {
      printf("|");
    }
    printf("%d", ret[i]);
  }
  printf("\n");
}

void dump(const char* name, const char* ret, int len, int maxChar) {
  printOffset(-len);
  printf("%s = ", name);
  for (int i = 0;; ++i) {
    if (i == len) {
      break;
    }
    if (i == maxChar) {
      printf("...");
      break;
    }
    if (ret[i] < 0x21 || ret[i] > 0x7e) {
      printf("\\%d", ret[i]);
    } else {
      printf("%c", ret[i]);
    }
  }
  printf("\n");
}

void dump(const char* name, int32_t ret) {
  printOffset(-4);
  printf("%s = %d", name, ret);
  printf("\n");
}

void dump(const char* name, uint32_t ret) {
  printOffset(-4);
  printf("%s = %ld", name, (int64_t)ret);
  printf("\n");
}

void dump(const char* name, uint64_t ret) {
  printOffset(-8);
  printf("%s = %zu", name, ret);
  printf("\n");
}

void dumpVirtualOffset(const char* name, uint64_t ret) {
  printOffset(-8);
  printf("%s = %zu (coffset = %zu, uoffset = %zu)", name, ret, ret >> 16,
         ret & ((1 << 16) - 1));
  printf("\n");
}

// intentionally do not delete[] for the sake of writing less codes
#define READ_WRITE_CHAR(var_type, var_name, var_len, max_out_len) \
  var_type* var_name = new var_type[var_len];                     \
  read(&br, var_name, var_len);                                   \
  dump(#var_name, var_name, var_len, max_out_len);                \
// delete[] var_name;

#define READ_WRITE_NUMBER(var_type, var_name) \
  var_type var_name;                          \
  read(&br, &(var_name));                     \
  dump(#var_name, var_name);

#define READ_WRITE_VIRTUAL_OFFSET(var_type, var_name) \
  var_type var_name;                                  \
  read(&br, &(var_name));                             \
  dumpVirtualOffset(#var_name, var_name);

// this program is based on:
// https://samtools.github.io/hts-specs/tabix.pdf
int main(int argv, const char** argc) {
  BufferedReader br(argc[1], 1024);
  READ_WRITE_CHAR(char, magic, 4, 50);
  const char TBI_MAGIC[] = "TBI\1";
  if (strncmp(magic, TBI_MAGIC, 4)) {
    fprintf(
        stderr,
        "Tabix magic is not found -- the input file is not tabix tbi file!\n");
  }

  READ_WRITE_NUMBER(int32_t, n_ref);
  READ_WRITE_NUMBER(int32_t, format);
  printf("  NOTE: format = 0: generic, 1: sam, 2: vcf\n");
  READ_WRITE_NUMBER(int32_t, col_seq);
  READ_WRITE_NUMBER(int32_t, col_beg);
  READ_WRITE_NUMBER(int32_t, col_end);
  READ_WRITE_NUMBER(int32_t, meta);
  printf("  NOTE: meta character is %c\n", meta);
  READ_WRITE_NUMBER(int32_t, skip);
  READ_WRITE_NUMBER(int32_t, l_nm);

  READ_WRITE_CHAR(char, names, l_nm, 50);

  for (int32_t ref_idx = 0; ref_idx < n_ref; ++ref_idx) {
    // char* names = new char[l_nm];
    // read(&br, names, l_nm);
    // dump("names", names, l_nm, 50);

    READ_WRITE_NUMBER(int32_t, n_bin);
    for (int32_t i = 0; i < n_bin; ++i) {
      READ_WRITE_NUMBER(uint32_t, bin);
      READ_WRITE_NUMBER(int32_t, n_chunk);
      for (int32_t j = 0; j < n_chunk; ++j) {
        READ_WRITE_VIRTUAL_OFFSET(uint64_t, chk_beg);
        READ_WRITE_VIRTUAL_OFFSET(uint64_t, chk_end);
      }
    }
    READ_WRITE_NUMBER(int32_t, n_intv);
    for (int32_t i = 0; i < n_intv; ++i) {
      READ_WRITE_VIRTUAL_OFFSET(uint64_t, ioff);
    }
  }

  if (!br.isEof()) {
    READ_WRITE_NUMBER(uint64_t, n_no_coor);
  }

  if (br.isEof()) {
    printOffset(0);
    printf("--FILE END--\n");
  }
  return 0;
}
