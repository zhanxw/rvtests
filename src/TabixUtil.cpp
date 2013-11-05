#include "TabixUtil.h"

#include "third/tabix/tabix.h"

int tabixIndexFile(const std::string& fn,
                   int skip,
                   char metaChar,
                   int chrom,
                   int startPos,
                   int endPos
                   ){

  // typedef struct {
  //       int32_t preset;
  //       int32_t sc, bc, ec; // seq col., beg col. and end col.
  //       int32_t meta_char, line_skip;
  // } ti_conf_t;
  ti_conf_t meta_conf = {0, chrom, startPos, endPos, metaChar, skip};
  return ti_index_build(fn.c_str(), &meta_conf);
}

