#include "BGenIndex.h"

#include <stdint.h>  // uint32_t
#include <cassert>
#include <string>
#include <vector>

int BGenIndex::init(const std::string& fn) {
  // printf("Open sqlite database: %s\n", fn.c_str());
  int rc = sqlite3_open(fn.c_str(), &db_);
  if (rc) {
    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_));
    // sqlite3_close(db_);
    return -1;
  }
  return 0;
}

int BGenIndex::queryRange(const std::string& chrom, int begin, int end) {
  char sql[512];
  snprintf(sql, 512,
           "SELECT file_start_position, size_in_bytes FROM Variant "
           "WHERE chromosome == '%s' AND position >= %d AND position < %d",
           chrom.c_str(), begin, end);
  // printf("SQL = \n%s\n", sql);
  int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt_, NULL);
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Can't prepare a SQL statement: %s\n", sqlite3_errmsg(db_));
    return (-1);
  }
  return 0;
}

int BGenIndex::setRange(const RangeList& r) {
  if (stmt_) {
    sqlite3_finalize(stmt_);
    stmt_ = NULL;
  }

  this->range.setRange(r);
  resetRangeIterator();

  this->queryRange(this->rangeIterator.getChrom(),
                   this->rangeIterator.getBegin(),
                   this->rangeIterator.getEnd());

  return 0;
}

void BGenIndex::resetRangeIterator() {
  this->rangeBegin = this->range.begin();
  this->rangeEnd = this->range.end();
  this->rangeIterator = this->range.begin();
}

/**
 *
 */
bool BGenIndex::next(int* file_start_position, int* size_in_bytes) {
  int rc = sqlite3_step(stmt_);

  while (rc != SQLITE_ROW) {
    // finish query
    if (rc == SQLITE_DONE) {
      sqlite3_finalize(stmt_);

      // find next range and query until sqlite_ROW
      ++rangeIterator;
      if (rangeIterator != this->rangeEnd) {
        this->queryRange(this->rangeIterator.getChrom(),
                         this->rangeIterator.getBegin(),
                         this->rangeIterator.getEnd());
        rc = sqlite3_step(stmt_);
      } else {
        return false;  // reach the end of range list
      }
    } else {
      // unsupported results
      fprintf(stderr, "Unhandled sqlite status [ %d ]: %s\n", rc,
              sqlite3_errmsg(db_));
      return false;
    }
  }

  assert(rc == SQLITE_ROW);  // has data to process
  *file_start_position = sqlite3_column_int(stmt_, 0);
  *size_in_bytes = sqlite3_column_int(stmt_, 1);
  return true;
}

BGenIndex::~BGenIndex() { sqlite3_close(db_); }
