#ifndef _BCFREADER_H_
#define _BCFREADER_H_

#include "bcf.h"
#include "RangeList.h"

static void write_header(bcf_hdr_t *h);

class BCFReader {
 public:
  BCFReader(const std::string& fn)
      : cannotOpen(false),
        inReading(false),
        hasIndex(false),
        b(0),
        bout(0),
        off(0),
        str2id(0){
    open(fn);
  };

  virtual ~BCFReader() {
    close();
  };

  bool readLine(std::string* line) {
    // openOK?
    if (cannotOpen) return false;

    // read
    if (!inReading) {
      resetRangeIterator();
      inReading = true;
    };

    // check read mode
    if (range.empty()) {
      // read line by line
      if (vcf_read(bp, hin, b) > 0) {
        vcf_write(bout, hout, b, line);
        return true;
      }
      return false;
    }

    // read by region
    // check index
    assert(!range.empty());
    if (!hasIndex) return false;

    if (off != 0) {
      // read a record
      while (vcf_read(bp, hin, b) > 0) {
        if (tid >= 0) {
          int l = strlen(b->ref);
          l = b->pos + (l > 0? l : 1);
          if (b->tid != tid || b->pos >= end) break; // current record has passed prespecified region
          if (!(l > begin && end > b->pos)) continue; // not sure when this will happen

          vcf_write(bout, hout, b, line);
          return true;
        }
      }
    }

    // seek to region
    for (; this->rangeIterator != this->rangeEnd; ++ rangeIterator) {
      char rangeBuffer[128];
      snprintf(rangeBuffer, 128, "%s:%u-%u", this->rangeIterator.getChrom().c_str(),
               this->rangeIterator.getBegin(), this->rangeIterator.getEnd());
      rangeBuffer[127] = '\0';
      // int tid, beg, end, len;
      if (!str2id) {
        str2id = bcf_build_refhash(hout);
      }
      if (bcf_parse_region(str2id, rangeBuffer, &tid, &begin, &end) >= 0) {
        off = bcf_idx_query(idx, tid, begin);
        if (off == 0) {
          // fprintf(stderr, "[%s] no records in the query region.\n", __func__);
          // return 1; // FIXME: a lot of memory leaks...
          continue;
        }
        bgzf_seek(bp->fp, off, SEEK_SET);

        while (vcf_read(bp, hin, b) > 0) {
          if (tid >= 0) {
            int l = strlen(b->ref);
            l = b->pos + (l > 0? l : 1);
            if (b->tid != tid || b->pos >= end) break; // current record has passed prespecified region
            if (!(l > begin && end > b->pos)) continue; // not sure when this will happen

            ++rangeIterator;
            vcf_write(bout, hout, b, line);
            return true;
          }
        }
      }
    }
    return false;
  };

  /**
   * @return 0 if adding region is valid
   */
  int addRange(const std::string& r) {
    if (inReading) {
      // don't allow updating region when reading starts
      return -1;
    }
    range.addRangeList(r.c_str());
    resetRangeIterator();
    return 0;
  };

  /**
   * Some ranges may be overlapping, thus we merge those
   */
  void mergeRange() {
    range.sort();
    resetRangeIterator();
  };

 private:
  // don't copy
  BCFReader(const BCFReader& );
  BCFReader& operator=(const BCFReader& );

 private:
  /* bool seekByRegion(region) { */

  /*   void *str2id = bcf_build_refhash(hout); */
  /*   if (bcf_parse_region(str2id, argv[optind+1], &tid, &begin, &end) >= 0) { */
  /*     bcf_idx_t *idx; */
  /*     idx = bcf_idx_load(argv[optind]); */
  /*     if (idx) { */
  /*       uint64_t off; */
  /*       off = bcf_idx_query(idx, tid, begin); */
  /*       if (off == 0) { */
  /*         fprintf(stderr, "[%s] no records in the query region.\n", __func__); */
  /*         return 1; // FIXME: a lot of memory leaks... */
  /*       } */
  /*       bgzf_seek(bp->fp, off, SEEK_SET); */
  /*       bcf_idx_destroy(idx); */
  /*     } */
  /*   } */

  /* } */


  bool openIndex(const std::string& fn) {
    idx = bcf_idx_load(fn.c_str());
    if (idx) {
      this->hasIndex = true;
    } else {
      this->hasIndex = false;
    }
    return this->hasIndex;
  };

  void closeIndex(){
    bcf_idx_destroy(idx);
  };


  int open(const std::string& fn);

  void close() {
    // destroy range iterator
    // close index
    bcf_hdr_destroy(hin);
    bcf_destroy(b); // bcf_destroy(blast);
    vcf_close(bp); // close bcf handle for input
    vcf_close(bout); // close bcf handle for output

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

  int vcf_write(bcf_t *bp, bcf_hdr_t *h, bcf1_t *b, std::string* line);

 private:
  // don't copy
  BCFReader(BCFReader& t);
  BCFReader& operator=(BCFReader& t);
 private:
  RangeList range;
  bool cannotOpen;
  bool inReading; // indicate reading has already started
  bool hasIndex;

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
};


#endif /* _BCFREADER_H_ */
