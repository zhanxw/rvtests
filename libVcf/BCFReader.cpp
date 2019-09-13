#include "BCFReader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "third/htslib/include/htslib/bgzf.h"
#include "third/htslib/include/htslib/kseq.h"
#include "third/htslib/include/htslib/kstring.h"
#include "third/htslib/include/htslib/vcf.h"

htsFile *my_vcf_open(const char *fn, const char *mode) {
  if (strchr(mode, 'b')) return hts_open(fn, mode);

  return hts_open(fn, mode);
}

#define my_vcf_read vcf_read

// adopted from vcf_write() in vcf.c
int my_vcf_write(htsFile *bp, bcf_hdr_t *h, bcf1_t *b, std::string *line) {
  // copied from bcf_write (equivalent to bcf_write1) from vcf.c in htslib
  if (h->dirty) bcf_hdr_sync(h);
  if (bcf_hdr_nsamples(h) != b->n_sample) {
    hts_log_error(
        "Broken VCF record, the number of columns at %s:%d does not match the "
        "number of samples (%d vs %d)",
        bcf_seqname(h, b), b->pos + 1, b->n_sample, bcf_hdr_nsamples(h));
    return -1;
  }

  // now copy from vcf_write
  kstring_t ks = {0, 0, 0};
  if (vcf_format1(h, b, &ks) != 0) {
    return -1;
  }
  (*line) = ks.s;
  return 0;
}

int BCFReader::open(const std::string &fn) {
  // check file existance
  bp = my_vcf_open(fn.c_str(), "rb");  // vcf file handle
  if (!bp) {
    this->cannotOpen = true;
    return -1;
  }
  b = bcf_init();
  if (b == 0) {
    return -1;
  }
  // write header
  hin = bcf_hdr_read(bp);

  // adapted from htslib vcf.c vcf_hdr_write()
  kstring_t htxt = {0, 0, 0};
  bcf_hdr_format(hin, 0, &htxt);
  while (htxt.l && htxt.s[htxt.l - 1] == '\0') --htxt.l;  // kill trailing zeros
  header = htxt.s;
  free(htxt.s);

  // open index
  this->hasIndex = this->openIndex(fn);

  // set up range iterator
  resetRangeIterator();

  cannotOpen = false;
  readyToRead = true;
  return 0;
}

bool BCFReader::readLine(std::string *line) {
  // openOK?
  if (cannotOpen) return false;

  // check read mode
  if (range.empty()) {
    // read line by line
    // if (my_vcf_read(bp, hin, b) > 0) {
    if (bcf_read(bp, hin, b) == 0) {
      if (my_vcf_write(bp, hin, b, line) == 0) {
        return true;
      }
    }
    return false;
  }

  // read by region
  // check index
  assert(!range.empty());
  if (!hasIndex) {
    this->readyToRead = false;
    return false;
  }

  if (iter) {
    // read a record
    while (bcf_itr_next(bp, iter, b) >= 0) {
      if (my_vcf_write(bp, hin, b, line) == 0) {
        return true;
      }
    }
  }

  // seek to region
  for (; this->rangeIterator != this->rangeEnd; ++rangeIterator) {
    char rangeBuffer[128];
    snprintf(rangeBuffer, 128, "%s:%u-%u",
             this->rangeIterator.getChrom().c_str(),
             this->rangeIterator.getBegin(), this->rangeIterator.getEnd());
    rangeBuffer[127] = '\0';

    iter = bcf_itr_querys(idx, hin, rangeBuffer);
    if (!iter) {
      continue;
    }
    while (bcf_itr_next(bp, iter, b) >= 0) {
      ++rangeIterator;
      if (my_vcf_write(bp, hin, b, line) == 0) {
        return true;
      }
    }
  }
  return false;
};
