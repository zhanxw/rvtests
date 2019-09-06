#include <stdlib.h>
#include <algorithm>

#include "base/MmapFile.h"
#include "base/RangeList.h"
#include "third/tabix/bgzf.h"

struct Record {
  int64_t pos;
  int64_t offset;
  bool operator<(Record& o) { return (pos < o.pos); }
};

bool comparator(const Record& a, const Record& b) { return a.pos < b.pos; }

class SingleChromosomeVCFIndex {
 public:
  SingleChromosomeVCFIndex(const std::string& vcfFile,
                           const std::string& indexFile);
  virtual ~SingleChromosomeVCFIndex();

  void close();
  void closeIndex();

  // @return 0 for success
  int createIndex();

  // open index
  int openIndex();

  // open index via mmap
  int mapIndex();

  // @return 1 if found, 0 not found, -1 if error
  int query(int chromPos, int64_t* voffset);

  // @return number of offsets found, -1 if error
  int query(int chromPosBeg, int chromPosEnd, int64_t* voffset);

  // @return length of @param line if success, or -1 if error
  int readLine(int64_t pos, std::string* line);

  int nextLine(std::string* line);

 private:
  std::string vcfFile_;  // must be bgzFile
  std::string indexFile_;
  void* data_;  // store indices
  MmapFile* mmapFile_;
  kstring_t* str;
  BGZF* fVcfFile_;
};

SingleChromosomeVCFIndex::SingleChromosomeVCFIndex(
    const std::string& vcfFile, const std::string& indexFile) {
  vcfFile_ = vcfFile;
  indexFile_ = indexFile;
  fVcfFile_ = bgzf_open(vcfFile_.c_str(), "rb");
  data_ = NULL;
  mmapFile_ = NULL;
  str = (kstring_t*)calloc(1, sizeof(kstring_t));
}

SingleChromosomeVCFIndex::~SingleChromosomeVCFIndex() { this->close(); }

void SingleChromosomeVCFIndex::close() {
  if (str) {
    free(str);
    str = NULL;
  }
  if (fVcfFile_) {
    bgzf_close(fVcfFile_);
    fVcfFile_ = NULL;
  }
  closeIndex();
}

void SingleChromosomeVCFIndex::closeIndex() {
  if (mmapFile_) {
    delete mmapFile_;
    data_ = NULL;
  }
  if (data_) {
    delete[](uint8_t*) data_;
    data_ = NULL;
  }
}

int SingleChromosomeVCFIndex::openIndex() {
  closeIndex();
  // read everything
  size_t fsize = getFileSize(indexFile_.c_str());
  printf("fsize = %ld\n", (long int)fsize);
  data_ = new uint8_t[fsize];
  FILE* fp = fopen(indexFile_.c_str(), "rb");
  if (fread(data_, sizeof(uint8_t), fsize, fp) != fsize) {
    printf("Read incomplete index\n");
    return -1;
  }

  // verify file integrity
  int64_t* d = (int64_t*)data_;
  if (fsize != sizeof(Record) * (2L + d[1])) {
    printf("Check file integrity!\n");
    printf("d = %ld %ld fsize = %ld\n", d[0], d[1], (long int)fsize);
    return -1;
  }
  return 0;
}

int SingleChromosomeVCFIndex::query(int chromPos, int64_t* pVirtualOffset) {
  return this->query(chromPos, chromPos, pVirtualOffset);
}

int SingleChromosomeVCFIndex::query(int chromPosBeg, int chromPosEnd,
                                    int64_t* voffset) {
  if (!data_) {
    printf("open index first!\n");
    return -1;
  }

  if (!voffset) {
    return -1;
  }
  printf("query [%d, %d]\n", chromPosBeg, chromPosEnd);

  Record* r = (Record*)data_;
  const int64_t Nrecord = r[0].offset;

  ++r;  // skip the first block, as first block is (#sample, #marker)

  // binary search for file position
  *voffset = -1;
  Record query;
  query.pos = chromPosBeg;
  // Comparator comparator;
  Record* lb =
      std::lower_bound(r, r + Nrecord + 1, query,
                       comparator);  // r[lb].pos >= query.pos = chromPosBeg
  query.pos = chromPosEnd;
  Record* ub =
      std::upper_bound(lb, r + Nrecord + 1, query,
                       comparator);  // r[ub].pos > query.pos = chromPosEnd
  printf("Found %d results\n", (int)(ub - lb));
  for (Record* pi = lb; pi != ub; ++pi) {
    // printf("%ld %ld\n", pi->pos, pi->offset);
    *voffset = lb->offset;
    break;
  }

  if (*voffset < 0) {
    printf("Cannot find position!\n");
    return -1;
  } else {
    printf("found %d position, e.g. %ld %ld\n", (int)(ub - lb), (*lb).pos,
           (*lb).offset);
    return ub - lb;
  }
}
int SingleChromosomeVCFIndex::readLine(int64_t offset, std::string* line) {
  if (bgzf_seek(fVcfFile_, offset, SEEK_SET)) {
    printf("seek error!\n");
  }
  kstring_t& s = *str;
  int ret = bgzf_getline(fVcfFile_, '\n', &s);
  if (ret <= 0) {
    printf("getline error, ret = %d!\n", ret);
  }
  // for (size_t i = 0; i < s.l; ++i) {
  //   if (i >= 50) break;
  //   printf("%c", s.s[i]);
  // }
  // printf("\n");

  *line = s.s;

  return s.l;
}

int SingleChromosomeVCFIndex::nextLine(std::string* line) {
  kstring_t& s = *str;
  int ret = bgzf_getline(fVcfFile_, '\n', &s);
  if (ret <= 0) {
    printf("getline error, ret = %d!\n", ret);
  }
  // for (size_t i = 0; i < s.l; ++i) {
  //   if (i >= 50) break;
  //   printf("%c", s.s[i]);
  // }
  // printf("\n");

  *line = s.s;
  return s.l;
}

int main(int argc, char** argv) {
  const char* fBG = argv[1];
  const char* fIndex = argv[2];
  // int64_t pos = strtol(argv[3], NULL, 0);
  const char* pos = argv[3];
  // const int Nrecord = 10;

  RangeList rl;
  rl.addRangeList(pos);

  SingleChromosomeVCFIndex sc(fBG, fIndex);
  sc.openIndex();
  int64_t offset;
  std::string ret;
  RangeList::iterator iter = rl.begin();
  for (; iter != rl.end(); ++iter) {
    printf("query: %d - %d \n", iter.getBegin(), iter.getEnd());
    int nFound = sc.query(iter.getBegin(), iter.getEnd(), &offset);
    if (nFound) {
      sc.readLine(offset, &ret);
      puts(ret.c_str());
    }
    nFound--;
    while (nFound > 0) {
      sc.nextLine(&ret);
      puts(ret.c_str());
      nFound--;
    }
  }
  return 0;
}
