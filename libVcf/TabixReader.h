#ifndef _TABIXREADER_H_
#define _TABIXREADER_H_

#include "base/RangeList.h"
#include "third/tabix/tabix.h"

class TabixReader {
 public:
  TabixReader(const std::string& fn)
      : cannotOpen(false),
        hasIndex(false),
        readyToRead(false),
        tabixHandle(0),
        ti_line(0) {
    open(fn);
  };

  virtual ~TabixReader() { close(); };

 private:
  // don't copy
  TabixReader(const TabixReader&);
  TabixReader& operator=(const TabixReader&);

 public:
  bool good() const { return this->readyToRead; }

  bool readLine(std::string* line) {
    // openOK?
    if (cannotOpen) return false;

    // read
    // check read mode
    if (range.empty()) {
      // read line by line
      if (!iter) {
        iter = ti_query(this->tabixHandle, 0, 0, 0);
        if (!iter) return false;
      }
      if (!this->firstLine.empty()) {
        (*line) = this->firstLine;
        this->firstLine.clear();
        return true;
      }
      while ((ti_line = ti_read(this->tabixHandle, iter, &ti_line_len)) != 0) {
        // need to skip header here
        if ((int)(*ti_line) == idxconf->meta_char) continue;
        (*line) = (ti_line);
        return true;
      }
      return false;
    }

    // read by region
    // check index
    assert(!range.empty());
    if (!hasIndex) {
      readyToRead = false;
      return false;
    }

    if (iter) {
      this->ti_line = ti_read(this->tabixHandle, iter, &ti_line_len);
      if (this->ti_line) {
        (*line) = ti_line;
        return true;
      }
    }

    // find valid iter
    for (; this->rangeIterator != this->rangeEnd; ++rangeIterator) {
      char rangeBuffer[128];
      snprintf(rangeBuffer, 128, "%s:%u-%u",
               this->rangeIterator.getChrom().c_str(),
               this->rangeIterator.getBegin(), this->rangeIterator.getEnd());
      rangeBuffer[127] = '\0';
      int tid, beg, end;  // , len;
      if (ti_parse_region(tabixHandle->idx, rangeBuffer, &tid, &beg, &end) !=
          0) {
        continue;
      }
      ti_iter_destroy(iter);
      iter = 0;
      this->iter = ti_queryi(this->tabixHandle, tid, beg, end);
      this->ti_line = ti_read(this->tabixHandle, this->iter, &ti_line_len);
      if (ti_line) {
        ++rangeIterator;
        (*line) = ti_line;
        return true;
      }
    }
    ti_iter_destroy(iter);
    iter = 0;

    return false;
  };

  /**
   * @return 0 if adding region is valid
   */
  int setRange(const std::string& r) {
    RangeList a;
    a.addRangeList(r);
    return this->setRange(a);
  }
  int setRange(const RangeList& r) {
    this->range.setRange(r);
    resetRangeIterator();
    if (this->iter) {
      ti_iter_destroy(iter);
      iter = 0;
    }
    return 0;
  }
  int addRange(const std::string& r) {
    RangeList a;
    a.addRangeList(r);
    return this->addRange(a);
  }
  int addRange(const RangeList& r) {
    this->range.addRange(r);
    resetRangeIterator();
    if (this->iter) {
      ti_iter_destroy(iter);
      iter = 0;
    }
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

 private:
  bool openIndex(const std::string& fn) {
    if (ti_lazy_index_load(this->tabixHandle) != 0) {
      // failed to open tabix index
      // fpritnf(stderr, "Cannot open index file for file [ %s ]!\n",
      // fn.c_str());
      this->hasIndex = false;
      return false;
    }

    this->hasIndex = true;
    return true;
  };

  void closeIndex() {
    // fpritnf(stderr, "close index...");
    if (!this->hasIndex) return;
    // fpritnf(stderr, "close index...");
    if (this->iter) {
      ti_iter_destroy(this->iter);
      this->iter = 0;
      // fpritnf(stderr, "close iter...");
    }

    /* fpritnf(stderr, "Close index\n"); */
    // fpritnf(stderr, "%x", this->tabixHandle);
    // fpritnf(stderr, "done. Close index\n");
  };

  int open(const std::string& fn) {
    ti_line = 0;

    // check file existance
    this->tabixHandle = ti_open(fn.c_str(), 0);
    if (!this->tabixHandle) {
      this->cannotOpen = true;
      return -1;
    }

    // open index
    this->hasIndex = this->openIndex(fn);

    // set up range iterator
    resetRangeIterator();

    // reset iterator
    this->iter = 0;

    // read header
    idxconf = ti_get_conf(this->tabixHandle->idx);

    // this command will let iter point to first record
    // (don't call this twice), or the internal fp pointer
    // will not re-shifted to its beginning position,
    // but will continue read.
    if (!this->hasIndex) {
      return -1;
    }

    this->iter = ti_query(this->tabixHandle, 0, 0, 0);
    while ((ti_line = ti_read(this->tabixHandle, this->iter,
                              &this->ti_line_len)) != 0) {
      if ((int)(*ti_line) != idxconf->meta_char) {
        this->firstLine = ti_line;
        break;
      }
      // fputs(ti_line, stdout); fputc('\n', stdout);
      this->header += ti_line;
      this->header += "\n";
    }

    cannotOpen = false;
    readyToRead = true;
    return 0;
  };

  void close() {
    // destroy range iterator
    // close index

    closeIndex();

    if (this->tabixHandle) {
      ti_close(this->tabixHandle);
      this->tabixHandle = 0;
      // fpritnf(stderr, "close handle...");
    }
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

  // tabix part
  tabix_t* tabixHandle;
  ti_iter_t iter;
  const char* ti_line;
  int ti_line_len;
  const ti_conf_t* idxconf;

  std::string header;
  std::string firstLine;
};

#endif /* _TABIXREADER_H_ */
