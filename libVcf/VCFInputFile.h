#ifndef _VCFINPUTFILE_H_
#define _VCFINPUTFILE_H_

// #include "tabix.h"
#include "VCFFilter.h"
#include "VCFRecord.h"

class TabixReader;
class BCFReader;

/**
 * parse is equivalent to copy, meaning we will copy the content to speed up
 * later reference using const char*
 * three modes to read a VCF files
 *  (1) BCF file mode, read by region or by line
 *  (2) VCF file, read by line
 *  (3) VCF file, read by region
 */
class VCFInputFile {
 public:
  typedef enum {
    BCF_MODE,
    VCF_LINE_MODE,  // read by line
    VCF_RANGE_MODE  // read by range
  } Mode;

 private:
  // disable copy-constructor
  VCFInputFile(const VCFInputFile&);
  VCFInputFile& operator=(const VCFInputFile&);

 public:
  VCFInputFile(const std::string& fn) { init(fn.c_str()); }

  VCFInputFile(const char* fn) { init(fn); }
  void init(const char* fn);

  virtual ~VCFInputFile() { this->close(); }
  void close();

  VCFHeader* getVCFHeader() {
    this->rewriteVCFHeader();
    return &this->header;
  };

  // use current subset of included people
  // to reconstruct a new VCF header line
  void rewriteVCFHeader();
  //  void setRangeMode();

  /**
   * Report a line that does not conform to VCF standard.
   */
  void reportReadError(const std::string& line) {
    if (line.size() > 50) {
      fprintf(stderr, "Error line [ %s ... ]\n", line.substr(0, 50).c_str());
    } else {
      fprintf(stderr, "Error line [ %s ]\n", line.c_str());
    }
  }
  /**
   * Check with VCFFilter to see if the current read line passed
   */
  virtual bool passFilter() { return true; };
  /**
   * Check if the current read VCF site is allowed
   */
  bool isAllowedSite() const {
    // no restriction on allowed sites
    if (this->allowedSite.empty()) return true;

    std::string chromPos = this->record.getChrom();
    chromPos += ":";
    chromPos += (this->record.getPosStr());
    if (this->allowedSite.count(chromPos)) return true;
    return false;
  }
  /**
   * @return true: a valid VCFRecord
   */
  bool readRecord();

  //////////////////////////////////////////////////
  // Sample inclusion/exclusion
  void includePeople(const char* s);
  void includePeople(const std::vector<std::string>& v);
  void includePeopleFromFile(const char* fn);
  void includeAllPeople();
  void excludePeople(const char* s);
  void excludePeople(const std::vector<std::string>& v);
  void excludePeopleFromFile(const char* fn);
  void excludeAllPeople();
  //////////////////////////////////////////////////
  // Adjust range collections
  void enableAutoMerge();
  void disableAutoMerge();
  // void clearRange();
  void setRangeFile(const char* fn);
  // @param l is a string of range(s)
  void setRange(const char* chrom, int begin, int end);
  void setRange(const RangeList& rl);
  void setRangeList(const std::string& l);
  // this function the entry point for all function add/change region list
  void setRangeList(const RangeList& rl);

  // which single-base chromosomal sites are allowed to read
  int setSiteFile(const std::string& fn);

  /**
   * @param fn load file and use first column as old id, second column as new id
   * @return number of ids have been changed.
   */
  int updateId(const char* fn);

  VCFRecord& getVCFRecord() { return this->record; };
  const char* getLine() const { return this->line.c_str(); };
  const char* getFileName() const { return this->fileName.c_str(); };
  void getIncludedPeopleName(std::vector<std::string>* p);

 private:
  void setRangeMode();

 private:
  VCFHeader header;
  VCFRecord record;

  // tabix_t * tabixHandle;
  // bool hasIndex;
  // RangeList range;

  std::string fileName;

  // variable used for accessing by rrange
  RangeList::iterator rangeBegin;
  RangeList::iterator rangeEnd;
  RangeList::iterator rangeIterator;
  // ti_iter_t iter;
  // const char* ti_line;
  bool autoMergeRange;

  Mode mode;
  std::string line;

  // readers
  LineReader* fp;
  TabixReader* tabixReader;
  BCFReader* bcfReader;

  // allow chromosomal sites
  std::set<std::string> allowedSite;
};

#endif /* _VCFINPUTFILE_H_ */
