#ifndef _BCFREADER_H_
#define _BCFREADER_H_

#include "base/RangeList.h"
#include "third/samtools/bcftools/bcf.h"

// static void write_header(bcf_hdr_t *h);

class BCFReader {
 public:
  BCFReader(const std::string& fn)
      : cannotOpen(false),
        hasIndex(false),
        readyToRead(false),
        b(0),
        bout(0),
        off(0),
        str2id(0) {
    open(fn);
  };

  virtual ~BCFReader() { close(); };

 private:
  // don't copy
  BCFReader(const BCFReader&);
  BCFReader& operator=(const BCFReader&);

 public:
  bool good() const { return this->readyToRead; }

  bool readLine(std::string* line);

  /**
   * @return 0 if adding region is valid
   */
  int setRange(const std::string& r) {
    RangeList a;
    a.addRangeList(r);
    return this->setRange(a);
  }
  int setRange(const RangeList& r) {
    range.setRange(r);
    resetRangeIterator();
    return 0;
  }
  int addRange(const std::string& r) {
    RangeList a;
    a.addRangeList(r);
    return this->addRange(a);
  }
  int addRange(const RangeList& r) {
    range.addRange(r);
    resetRangeIterator();
    return 0;
  }

  /**
   * Some ranges may be overlapping, thus we merge those
   */
  void mergeRange() {
    range.sort();
    resetRangeIterator();
  };

  const std::string& getHeader() const { return this->header; }

  bool indexed() const { return this->hasIndex; }

 private:
  bool openIndex(const std::string& fn) {
    idx = bcf_idx_load(fn.c_str());
    if (idx) {
      this->hasIndex = true;
    } else {
      this->hasIndex = false;
    }
    return this->hasIndex;
  };

  void closeIndex() { bcf_idx_destroy(idx); };

  int open(const std::string& fn);

  void close() {
    // destroy range iterator
    // close index
    bcf_hdr_destroy(hin);
    bcf_destroy(b);   // bcf_destroy(blast);
    vcf_close(bp);    // close bcf handle for input
    vcf_close(bout);  // close bcf handle for output
    /* // resume stdout */
    /* stdout = fdopen(this->origStdout, "w"); */
    /* assert(stdout); */

    if (str2id) {
      bcf_str2id_destroy(str2id);
    }
    closeIndex();
  };
  void resetRangeIterator() {
    this->rangeBegin = this->range.begin();
    this->rangeEnd = this->range.end();
    this->rangeIterator = this->range.begin();
  }

 private:
  RangeList range;
  bool cannotOpen;
  bool hasIndex;
  bool readyToRead;

  // variable used for accessing by range
  RangeList::iterator rangeBegin;
  RangeList::iterator rangeEnd;
  RangeList::iterator rangeIterator;

  // bcftools part
  bcf_t* bp;
  bcf1_t* b;
  bcf_t* bout;
  bcf_hdr_t* hin;
  bcf_hdr_t* hout;
  bcf_idx_t* idx;
  int tid, begin, end;
  uint64_t off;
  void* str2id;

  /* BCF_t* BCFHandle; */
  /* ti_iter_t iter; */
  const char* line;
  int line_len;
  // int origStdout;
  std::string header;
};

#endif /* _BCFREADER_H_ */
