#ifndef _BGENINDEX_H_
#define _BGENINDEX_H_

#include <string>

#include "third/sqlite/include/sqlite3.h"

#include "base/RangeList.h"

class BGenIndex {
  /**
   * definition of BGEN index
   * https://bitbucket.org/gavinband/bgen/wiki/The_bgenix_index_file_format
CREATE TABLE Variant (
  chromosome TEXT NOT NULL,
  position INT NOT NULL,
  rsid TEXT NOT NULL,
  number_of_alleles INT NOT NULL,
  allele1 TEXT NOT NULL,
  allele2 TEXT NULL,
  file_start_position INT NOT NULL,
  size_in_bytes INT NOT NULL,
  PRIMARY KEY (chromosome, position, rsid, allele1, allele2,
  file_start_position )
  ) WITHOUT ROWID;
   */
 public:
  BGenIndex() : db_(NULL), stmt_(NULL) {}
  BGenIndex(const std::string& fn) : db_(NULL), stmt_(NULL) { init(fn); }
  int init(const std::string& fn);
  /**
   * Contruct a SQL based on the specified range chrom:begin-end
   * @param int begin inclusive
   * @param int end exclusive
   */
  int setRange(const RangeList& r);
  /**
   * @return true if there is a genotype probability data block
   */
  bool next(int* file_start_position, int* size_in_bytes);
  ~BGenIndex();

 private:
  BGenIndex(const BGenIndex&);
  BGenIndex& operator=(const BGenIndex&);

 private:
  int queryRange(const std::string& chrom, int begin, int end);
  void resetRangeIterator();

 private:
  sqlite3* db_;
  sqlite3_stmt* stmt_;

 private:
  // variable used for accessing by range
  RangeList range;
  RangeList::iterator rangeBegin;
  RangeList::iterator rangeEnd;
  RangeList::iterator rangeIterator;
};

#endif /* _BGENINDEX_H_ */
