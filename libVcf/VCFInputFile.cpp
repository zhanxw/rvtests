#include "VCFInputFile.h"

#include "base/IO.h"
#include "base/Utils.h"

#include "BCFReader.h"
#include "TabixReader.h"

// use current subset of included people
// to reconstruct a new VCF header line
void VCFInputFile::rewriteVCFHeader() {
  std::string s = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  VCFPeople& people = this->record.getPeople();
  for (unsigned int i = 0; i < people.size(); i++) {
    s += '\t';
    s += people[i]->getName();
  }
  this->header[this->header.size() - 1] = s;
}

void VCFInputFile::setRangeMode() {
  if (mode == VCF_LINE_MODE) {
    this->tabixReader = new TabixReader(this->fileName);
    if (!this->tabixReader->good()) {
      fprintf(stderr,
              "[ERROR] Cannot read VCF by range, please verify you have the "
              "index file"
              "(or create one using tabix).\nQuitting...");
      abort();
    } else {
      this->mode = VCFInputFile::VCF_RANGE_MODE;
    }
  } else if (mode == VCF_RANGE_MODE) {
    // Auto-merge should be handled by VCFInputFile, not in tabixReader
    // if (this->autoMergeRange) {
    //   this->tabixReader->enableAutoMerge();
    // }
  } else if (mode == BCF_MODE) {
    if (!this->bcfReader->good() || !this->bcfReader->indexed()) {
      fprintf(stderr,
              "[ERROR] Cannot read BCF by range, please verify you have the "
              "index file "
              "(or create one using bcftools).\nQuitting...");
      abort();
    }
    // Auto-merge should be handled by VCFInputFile, not in bcfReader
    // if (this->autoMergeRange) {
    //   this->bcfReader->enableAutoMerge();
    // }
  }

  // if (this->autoMergeRange) {
  //   this->range.sort();
  // }

  // this->rangeBegin = this->range.begin();
  // this->rangeEnd = this->range.end();
  // this->rangeIterator = this->range.begin();
}

// void VCFInputFile::clearRange() {
// #ifndef NDEBUG
//   if (this->range.size()) {
//     fprintf(stderr, "Clear existing %zu range.\n", this->range.size());
//   }
// #endif
//   if (mode == BCF_MODE) {
//     this->bcfReader->clearRange();
//   } else if (mode == VCF_RANGE_MODE) {
//     this->VCFRecord->clearRange();
//   }
//   // this->range.clear();
//   // this->ti_line = 0;
// };

/**
 * @param fn: the file contains two column: old_id new_id
 */
int VCFInputFile::updateId(const char* fn) {
  // load conversion table
  LineReader lr(fn);
  std::map<std::string, std::string> tbl;
  std::vector<std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")) {
    if (tbl.find(fd[0]) != tbl.end()) {
      fprintf(stderr,
              "Duplicated original ids: [ %s ], replace it to new id anyway.\n",
              fd[0].c_str());
    };
    if (fd.empty() || fd[0].empty() || fd.size() < 2) continue;
    tbl[fd[0]] = fd[1];
  }

  // rewrite each people's name
  std::string s = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  VCFPeople& people = this->record.getPeople();
  int n = 0;
  for (unsigned int i = 0; i < people.size(); i++) {
    if (tbl.find(people[i]->getName()) != tbl.end()) {
      ++n;
      people[i]->setName(tbl[people[i]->getName()]);
    }
  }
  this->rewriteVCFHeader();

  // return result
  return n;
}

void VCFInputFile::init(const char* fn) {
  this->fileName = fn;
  this->fp = NULL;
  this->tabixReader = NULL;
  this->bcfReader = NULL;
  this->autoMergeRange = false;

  // check whether file exists.
  FILE* fp = fopen(fn, "rb");
  if (!fp) {
    fprintf(stderr, "[ERROR] Failed to open file [ %s ]\n", fn);
    // quit to avoid segfault in the future
    exit(1);
  }
  fclose(fp);

  bool headerLoaded = false;
  // use file name to check file type
  if (endsWith(fn, ".bcf") || endsWith(fn, ".bcf.gz")) {
    this->mode = BCF_MODE;
    this->bcfReader = new BCFReader(fn);
    const std::string& h = this->bcfReader->getHeader();
    this->header.setHeader(h);
    this->record.createIndividual(this->header[this->header.size() - 1]);
    headerLoaded = true;
  } else {
    if (!endsWith(fn, ".vcf") && !endsWith(fn, "vcf.gz")) {
      fprintf(stderr, "[WARN] File name does not look like a VCF/BCF file.\n");
    }
    this->mode = VCF_LINE_MODE;
    this->fp = new LineReader(fn);

    // open file
    // read header
    while (this->fp->readLine(&line)) {
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
    // this->hasIndex = this->openIndex();
  }
  if (headerLoaded == false) {
    FATAL("VCF/BCF File does not have header!");
  }
  // this->clearRange();
}

void VCFInputFile::close() {
  // closeIndex();
  this->record.deleteIndividual();
  if (this->fp) {
    delete this->fp;
    this->fp = NULL;
  }
  if (this->tabixReader) {
    delete this->tabixReader;
    this->tabixReader = NULL;
  }
  if (this->bcfReader) {
    delete this->bcfReader;
    this->bcfReader = NULL;
  }
}

bool VCFInputFile::readRecord() {
  int nRead = 0;
  while (true) {
    if (this->mode == VCF_LINE_MODE) {
      nRead = this->fp->readLine(&this->line);
    } else if (this->mode == VCF_RANGE_MODE) {
      nRead = this->tabixReader->readLine(&this->line);
    } else if (this->mode == BCF_MODE) {
      nRead = this->bcfReader->readLine(&this->line);
    }
    if (!nRead) return false;

    // star parsing
    int ret;
    this->record.attach(&this->line);
    ret = this->record.parseSite();
    if (ret) {
      reportReadError(this->line);
    }
    if (!this->isAllowedSite()) continue;

    ret = this->record.parseIndividual();
    if (ret) {
      reportReadError(this->line);
    }
    if (!this->passFilter()) continue;

    // break;
    return true;
  }
}

//////////////////////////////////////////////////
// Sample inclusion/exclusion
void VCFInputFile::includePeople(const char* s) {
  this->record.includePeople(s);
}
void VCFInputFile::includePeople(const std::vector<std::string>& v) {
  this->record.includePeople(v);
}
void VCFInputFile::includePeopleFromFile(const char* fn) {
  this->record.includePeopleFromFile(fn);
}
void VCFInputFile::includeAllPeople() { this->record.includeAllPeople(); }
void VCFInputFile::excludePeople(const char* s) {
  this->record.excludePeople(s);
}
void VCFInputFile::excludePeople(const std::vector<std::string>& v) {
  this->record.excludePeople(v);
}
void VCFInputFile::excludePeopleFromFile(const char* fn) {
  this->record.excludePeopleFromFile(fn);
}
void VCFInputFile::excludeAllPeople() { this->record.excludeAllPeople(); }

//////////////////////////////////////////////////
// Adjust range collections
void VCFInputFile::enableAutoMerge() { this->autoMergeRange = true; }
void VCFInputFile::disableAutoMerge() { this->autoMergeRange = false; }
// void clearRange();
void VCFInputFile::setRangeFile(const char* fn) {
  if (!fn || strlen(fn) == 0) return;
  RangeList r;
  r.addRangeFile(fn);
  this->setRange(r);
}
// @param l is a string of range(s)
void VCFInputFile::setRange(const char* chrom, int begin, int end) {
  RangeList r;
  r.addRange(chrom, begin, end);
  this->setRange(r);
}
void VCFInputFile::setRange(const RangeList& rl) { this->setRangeList(rl); }
void VCFInputFile::setRangeList(const std::string& l) {
  if (l.empty()) return;

  RangeList r;
  r.addRangeList(l);
  this->setRange(r);
}
// this function the entry point for all function add/change region list
void VCFInputFile::setRangeList(const RangeList& rl) {
  if (rl.size() == 0) return;

  this->setRangeMode();

  RangeList l;
  l.setRange(rl);
  if (this->autoMergeRange) l.sort();

  if (mode == VCF_RANGE_MODE) {
    this->tabixReader->setRange(l);
  } else if (mode == BCF_MODE) {
    this->bcfReader->setRange(l);
  } else {
    fprintf(stderr, "[ERROR] invalid reading mode, quitting...\n");
    abort();
  }
}

int VCFInputFile::setSiteFile(const std::string& fn) {
  if (fn.empty()) return 0;

  std::vector<std::string> fd;
  LineReader lr(fn);
  int pos;
  std::string chromPos;
  while (lr.readLineBySep(&fd, "\t ")) {
    if (fd.empty()) continue;
    if (fd[0].find(':') != std::string::npos) {
      this->allowedSite.insert(fd[0]);
      continue;
    }
    if (fd.size() >= 2 && str2int(fd[1], &pos) && pos > 0) {
      chromPos = fd[0];
      chromPos += ':';
      chromPos += fd[1];
      this->allowedSite.insert(chromPos);
      continue;
    }
  }
  return 0;
}

void VCFInputFile::getIncludedPeopleName(std::vector<std::string>* p) {
  record.getIncludedPeopleName(p);
}
