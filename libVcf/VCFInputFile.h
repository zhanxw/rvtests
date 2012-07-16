#ifndef _VCFINPUTFILE_H_
#define _VCFINPUTFILE_H_

#include "tabix.h"
#include "VCFRecord.h"
#include "VCFFilter.h"

/**
 * parse is equivalent to copy, meaning we will copy the content so speed up later reference by const char*
 */
class VCFInputFile{
public:
  typedef enum {
    LINE_MODE,           // read line by line
    RANGE_MODE          // read range by range
  } Mode;

public:
  VCFInputFile (const std::string& fn) {
    init(fn.c_str());
  }
  VCFInputFile (const char* fn) {
    init(fn);
  }
  void init(const char* fn) {
    this->fp = NULL;
    this->tabixHandle = NULL;
    this->mode = LINE_MODE;
    bool headerLoaded;

    this->fileName = fn;
    // open file
    this->fp = new LineReader(fn);

    // read header
    while (this->fp->readLine(&line)){
      if (line[0] == '#') {
        this->header.push_back(line);
        if (line.substr(0, 6) == "#CHROM") {
          this->record.createIndividual(line);
          headerLoaded = true;
          break;
        }
        continue;
      }
      if (line[0] != '#') {
        FATAL("Wrong VCF header");
      }
    }
    if (headerLoaded == false) {
      FATAL("VCF File does not have header!");
    }

    this->tabixHandle = 0;
    this->iter = 0;
    this->hasIndex = this->openIndex();

    this->clearRange();
    this->ti_line = 0;
  };
  ~VCFInputFile(){
    this->close();
  };
  void close() {
    closeIndex();
    this->record.deleteIndividual();
    if (this->fp) {
      delete this->fp;
      this->fp = NULL;
    }
  };
  VCFHeader* getVCFHeader() {
    this->rewriteVCFHeader();
    return &this->header;
  };
  bool openIndex() {
    return this->openIndex( (fileName + ".tbi").c_str());
  };
  bool openIndex(const char* fn) {
    if (( this->tabixHandle = ti_open(this->fileName.c_str(), 0)) == 0 ) {
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

  // use current subset of included people
  // to reconstruct a new VCF header line
  void rewriteVCFHeader();
  void setRangeMode();
  void clearRange();

  void setRangeFile(const char* fn) {
    if (!fn || strlen(fn) == 0)
      return;

    this->clearRange();
    this->range.addRangeFile(fn);
    this->setRangeMode();
  };
  // @param l is a string of range(s)
  void setRange(const char* chrom, int begin, int end) {
    this->clearRange();
    this->range.addRange(chrom, begin, end);
    this->setRangeMode();
  };
  void setRange(const RangeList& rl) {
    this->clearRange();
    this->range = rl;
    this->setRangeMode();
  };
  void setRangeList(const char* l){
    if (!l || strlen(l) == 0)
      return;

    this->clearRange();
    this->range.addRangeList(l);
    this->setRangeMode();
  };
  void setRangeList(RangeList& rl){
    if (rl.size() == 0) return;

    this->clearRange();
    this->range = rl;
    this->setRangeMode();
  };

  ///
  void reportReadError(const std::string& line) {
    if (line.size() > 50) {
      fprintf(stderr, "Error line [ %s ... ]", line.substr(0, 50).c_str());
    } else {
      fprintf(stderr, "Error line [ %s ]", line.c_str());
    }
  }
  /**
   * Check with VCFFilter to see if the current read line passed
   */
  virtual bool passFilter() {
    return true;
  };
  /**
   * @return true: a valid VCFRecord
   */
  bool readRecord(){
    int ret;
    // load contents
    if (this->mode == VCFInputFile::RANGE_MODE) {
      while (this->rangeIterator != this->rangeEnd) {
        if (!this->ti_line) { // last time does not read a valid line
          // get range
          char rangeBuffer[128];
          snprintf(rangeBuffer, 128, "%s:%u-%u", this->rangeIterator.getChrom().c_str(),
                   this->rangeIterator.getBegin(), this->rangeIterator.getEnd());
          rangeBuffer[127] = '\0';
#ifndef NDEBUG
          fprintf(stderr, "Process range: %s\n", rangeBuffer);
          // this->range.dump();
#endif        
          // parse range
          int tid, beg, end, len;
          if (ti_parse_region(tabixHandle->idx, rangeBuffer, &tid, &beg, &end) != 0){
#ifndef NDEBUG            
            fprintf(stderr, "Maybe non-existing range: %s, pass....\n", rangeBuffer);
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
            this->line = ti_line;
            ret = this->record.parse(this->line);
            if (ret) {
              reportReadError(this->line);
            }
            if (!this->passFilter())
              continue;
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
            this->line = ti_line;
            ret = this->record.parse(this->line);
            if (ret) {
              reportReadError(this->line);
            }
            if (!this->passFilter())
              continue;
            return true;
          }
        }
      } // end while
      return false;
    } else { // go line by line
      while (true) {
        if (this->fp->readLine(&this->line)){
          ret = this->record.parse(this->line);
          if (ret) {
            reportReadError(this->line);
          }
          if (!this->passFilter())
            continue;
          return true;
        } else {
          return false; // reach file end
        }
      }
    }
  };
  void includePeople(const char* s) {
    this->record.includePeople(s);
  };
  void includePeople(const std::vector<std::string>& v){
    this->record.includePeople(v);
  };
  void includePeopleFromFile(const char* fn) {
    this->record.includePeopleFromFile(fn);
  };
  void includeAllPeople(){
    this->record.includeAllPeople();
  };
  void excludePeople(const char* s) {
    this->record.excludePeople(s);
  };
  void excludePeople(const std::vector<std::string>& v){
    this->record.excludePeople(v);
  };
  void excludePeopleFromFile(const char* fn) {
    this->record.excludePeopleFromFile(fn);
  };
  void excludeAllPeople(){
    this->record.excludeAllPeople();
  };
  /**
   * @return number of ids have been changed.
   */
  int updateId(const char* fn);

  VCFRecord& getVCFRecord() {return this->record;};
  const char* getLine() const {return this->line.c_str();};
  const char* getFileName() const {return this->fileName.c_str();};
private:
  // disable copy-constructor
  VCFInputFile(const VCFInputFile& v);
  VCFInputFile& operator=(const VCFInputFile& v);

private:
  VCFHeader header;
  VCFRecord record;

  LineReader* fp;
  tabix_t * tabixHandle;
  bool hasIndex;
  RangeList range;

  std::string fileName;

  // variable used for accessing by rrange
  RangeList::iterator rangeBegin;
  RangeList::iterator rangeEnd;
  RangeList::iterator rangeIterator;
  ti_iter_t iter;
  const char* ti_line;

  Mode mode;
  std::string line;
};

#endif /* _VCFINPUTFILE_H_ */
