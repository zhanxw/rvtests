#ifndef _BGENFILE_H_
#define _BGENFILE_H_

#include <stdint.h>  // uint32_t
#include <set>
#include <string>
#include <vector>

#include "third/zlib/zlib.h"
#include "third/zstd/lib/zstd.h"

#include "BGenIndex.h"
#include "BGenVariant.h"

// copied from libVcf/VCFConstant.h
#define MISSING_GENOTYPE -9
#define PLINK_MALE 1
#define PLINK_FEMALE 2

class RangeList;

class BGenFile {
 public:
  enum SNP_COMPRESSION { NO_COMPRESSION = 0, GZIP = 1, ZSTD = 2 };
  enum LAYOUT { UNSUPPORTED = 0, LAYOUT1 = 1, LAYOUT2 = 2 };
  enum SAMPLE_IDENTIFIER {
    NO_SAMPLE_IDENTIFIER = 0,
    HAS_SAMPLE_IDENTIFIER = 1
  };
  typedef enum {
    BGEN_LINE_MODE,  // read by line
    BGEN_RANGE_MODE  // read by range
  } Mode;

  BGenFile(const std::string& fn);
  /**
   * @return true: if a valid record is read
   */
  bool readRecord();

  //////////////////////////////////////////////////
  // Sample inclusion/exclusion
  void includePeople(const std::string& s);
  void includePeople(const std::vector<std::string>& v);
  void includePeopleFromFile(const char* fn);
  void includeAllPeople();
  void excludePeople(const std::string& s);
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
  // void setRangeMode();

  /**
   * Load sample file that accompanies bgen file
   * @param fn input sample file name.
   * Do nothing when @param fn is empty
   */
  int loadSampleFile(const std::string& fn);

 public:
  int getNumMarker() const { return M; }
  int getNumSample() const { return N; }
  int getNumEffectiveSample() const;
  const std::vector<std::string>& getSampleIdentifier() const {
    return sampleIdentifier;
  }
  void getIncludedSampleName(std::vector<std::string>* p) const;
  const BGenVariant& getVariant() const { return var; }
  int getEffectiveIndex(int idx) const;
  void printInfo();

 private:
  BGenFile(const BGenFile&);
  BGenFile& operator=(const BGenFile&);

 private:
  bool parseLayout1();
  bool parseLayout2();

  void parseString(FILE* fp, int lenByte, std::string* out);
  void parseUint32(FILE* fp, uint32_t* value);
  void parseUint16(FILE* fp, uint16_t* value);
  int choose(int n, int m);

  bool isFileEnd(FILE* fp);
  static long getFileSize(const std::string& fn);

  // sample inclusion/exclusion related
  void setPeopleMask(const std::string& s, bool b);
  void setPeopleMaskFromFile(const char* fn, bool b);
  void setRangeMode();
  // range list related
  void buildEffectiveIndex();

 private:
  std::string bgenFileName;
  FILE* fp;
  // first 4 bytes
  uint32_t offset;
  // header block
  uint32_t LH;
  uint32_t M;
  uint32_t N;
  uint8_t magic[4];
  std::vector<uint8_t> freeData;
  uint32_t flag;
  SNP_COMPRESSION snpCompression;
  LAYOUT layout;
  SAMPLE_IDENTIFIER flagSampleIdentifier;
  // sample identifier block
  std::vector<std::string> sampleIdentifier;

  long fileSize;
  std::vector<uint8_t> compressedBuf;
  uint32_t C;  // number of compressed bytes
  std::vector<uint8_t> buf;
  uint32_t D;  // number of bytes before compression

  BGenVariant var;
  BGenIndex index;

  bool autoMergeRange;
  Mode mode;                        /// read consecutively or read by index
  std::vector<bool> sampleMask;     // true means exclusion
  std::vector<int> effectiveIndex;  // index of unmasked samples
  // allow chromosomal sites
  std::set<std::string> allowedSite;
};  // class BGenFile

#endif /* _BGENFILE_H_ */
