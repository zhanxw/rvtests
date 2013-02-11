#ifndef _TABIXREADER_H_
#define _TABIXREADER_H_

#include "third/tabix/tabix.h"
#include "RangeList.h"

class TabixReader{
public:
TabixReader(const std::string& fn) : inReading(false) {
    open(fn);
  };
  ~TabixReader() {
    close();
  };

  bool openIndex(const std::string& fn) {
    if (( this->tabixHandle = ti_open(fn.c_str(), 0)) == 0 ) {
      // failed to open tabix index
      this->hasIndex = false;
      return false;
    } else{
      if (ti_lazy_index_load(this->tabixHandle) != 0) {
        // failed to open tabix index
        this->hasIndex = false;
        return false;
      } else {
        this->hasIndex = true;
        return true;
      }
    }
    return true;
  };
  void closeIndex(){
    if (this->iter) {
      ti_iter_destroy(this->iter);
      this->iter = 0;
    }
    ti_close(this->tabixHandle);
    this->tabixHandle = 0;
  };

  bool readLine(std::string* line) {
    // check index
    if (!hasIndex) return false;

    // read
    if (!inReading) {
      resetRangeIterator();
      inReading = true;
    };

    while (this->rangeIterator != this->rangeEnd) {
      if (!this->ti_line) { // last time does not read a valid line
        // get range
        char rangeBuffer[128];
        snprintf(rangeBuffer, 128, "%s:%u-%u", this->rangeIterator.getChrom().c_str(),
                 this->rangeIterator.getBegin(), this->rangeIterator.getEnd());
        rangeBuffer[127] = '\0';
#if 0
        REprintf("Process range: %s\n", rangeBuffer);
        // this->range.dump();
#endif
        // parse range
        int tid, beg, end, len;
        if (ti_parse_region(tabixHandle->idx, rangeBuffer, &tid, &beg, &end) != 0){
#if 0
          REprintf("Maybe non-existing range: %s, pass....\n", rangeBuffer);
#endif
          // continue to next rangeIdx
          ti_iter_destroy(this->iter);
          this->iter = 0;
          ++ this->rangeIterator;
          continue;
          // FATAL("Cannot ti_parse_region");
        }
        this->iter =  ti_queryi(tabixHandle, tid, beg, end);
        this->ti_line = ti_read(this->tabixHandle, iter, &len);
        if (this->ti_line) { // s is valid
          (*line) = ti_line;
          return true;
        } else{
          // continue to next rangeIdx
          ti_iter_destroy(this->iter);
          this->iter = 0;
          ++ this->rangeIterator;
          continue;
        }
      } else {  // last time read a valid line
        int len;
        this->ti_line = ti_read(this->tabixHandle, iter, &len);
        if (!this->ti_line) {
          ++ this->rangeIterator;
          continue;
        } else {
          (*line) = ti_line;
          return true;
        }
      }
    } // end while
    return false;
  };

  /**
   * @return 0 if adding region is valid
   */
  int addRange(const std::string& r) {
    if (inReading)
      return -1;
    range.addRangeList(r.c_str());
    resetRangeIterator();
    return 0;
  };

  /**
   * Some ranges may be overlapping, thus we merge those
   */
  void mergeRange() {
    range.sort();
  };
  int open(const std::string& fn) {
    inReading = false;
    ti_line = 0;
    
    //open index
    this->tabixHandle = 0;
    this->iter = 0;
    this->hasIndex = this->openIndex(fn);

    // set up range iterator
    resetRangeIterator();

    return this->hasIndex ? 0 : -1;
  };
  
  void close() {
    // destroy range iterator
    // close index

    closeIndex();
  };
  void resetRangeIterator() {
    this->rangeBegin = this->range.begin();
    this->rangeEnd = this->range.end();
    this->rangeIterator = this->range.begin();
  }

private:
  // don't copy
  TabixReader(TabixReader& t);
  TabixReader& operator=(TabixReader& t);
private:
  RangeList range;
  bool inReading; // indicate reading has already started
  bool hasIndex;

  // variable used for accessing by range
  RangeList::iterator rangeBegin;
  RangeList::iterator rangeEnd;
  RangeList::iterator rangeIterator;

  // tabix part
  tabix_t * tabixHandle;
  ti_iter_t iter;
  const char* ti_line;

};

#endif /* _TABIXREADER_H_ */
