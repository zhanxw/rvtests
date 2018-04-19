#include "KGGInputFile.h"

#include "base/Exception.h"
#include "base/IO.h"
#include "base/TypeConversion.h"

#ifndef UNUSED
#define UNUSED(x) (void)(x)
#endif

#define ROUND_UP_TO_4X(x) (((x) + 3) & ~0x03)

KGGInputFile::KGGInputFile(const std::string& fnPrefix) { init(fnPrefix, ""); }

KGGInputFile::KGGInputFile(const std::string& fnPrefix,
                           const std::string& fnSuffix) {
  if (!fnSuffix.empty() && fnSuffix[0] != '.') {
    std::string s = ".";
    s += fnSuffix;
    init(fnPrefix, s);
  } else {
    init(fnPrefix, fnSuffix);
  }
}

int KGGInputFile::init(const std::string& fnPrefix,
                       const std::string& fnSuffix) {
  this->prefix = fnPrefix;
  this->fpKed = new BufferedReader((prefix + ".ked" + fnSuffix).c_str(), 1024);
  this->fpKim = new BufferedReader((prefix + ".kim" + fnSuffix).c_str(), 1024);
  this->fpKam = new BufferedReader((prefix + ".kam" + fnSuffix).c_str(), 1024);
  if (!this->fpKed || !this->fpKim || !this->fpKam) {
    REPORT("Cannot open binary KGG file!");
    abort();
  }
  // read KED header
  int c;
  // magic number
  unsigned char magic[] = {0x9e, 0x82};
  c = this->fpKed->getc();
  if ((unsigned char)c != magic[0]) {
    fprintf(stderr, "Magic number of binary KGG file does not match!\n");
    abort();
  }
  c = this->fpKed->getc();
  if ((unsigned char)c != magic[1]) {
    fprintf(stderr, "Magic number of binary KGG file does not match!\n");
    abort();
  }

  // snp major mode
  const int PHASED_MODE = 0x01;  // 0b00000001;
  const int UNPHASED_MODE = 0x00;
  c = this->fpKed->getc();
  if (c == UNPHASED_MODE) {
    this->phased = false;
  } else if (c == PHASED_MODE) {
    this->phased = true;
  } else {
    fprintf(stderr, "Unrecognized major mode in binary KGG file.\n");
    exit(1);
  }

  // read bim
  LineReader* lr = new LineReader((this->prefix + ".kim" + fnSuffix).c_str());
  std::string chrPos;
  std::vector<std::string> fd;
  while (lr->readLineBySep(&fd, " \t")) {
    if (fd.size() != 7) {
      fprintf(stderr, "Wrong format in bim file.\n");
      exit(1);
    }

    // when rsid == ".", use "chr:pos" as dict key
    if (fd[1] == ".") {
      chrPos = fd[0];
      chrPos += ":";
      chrPos += fd[3];
    } else {
      chrPos = fd[1];
    }
    if (snp2Idx.find(chrPos) == snp2Idx.end()) {
      chrom.push_back(fd[0]);
      snp.push_back(fd[1]);
      snp2Idx[chrPos] = 0;
      const int val = snp2Idx.size() - 1;
      snp2Idx[chrPos] = val;
      mapDist.push_back(atof(fd[2].c_str()));
      pos.push_back(atoi(fd[3].c_str()));
      ref.push_back(fd[4]);
      alt.push_back(stringTokenize(fd[5], ","));
    } else {
      fprintf(stderr,
              "Error found: duplicated marker name or chromosomal position [ "
              "%s ]!\n",
              fd[1].c_str());
      exit(1);
    }
  }
  delete lr;

  // read kam
  lr = new LineReader((this->prefix + ".kam" + fnSuffix).c_str());
  while (lr->readLineBySep(&fd, " \t")) {
    if (fd.size() != 6) {
      fprintf(stderr, "Wrong format in kam file.\n");
      exit(1);
    }

    // will skip loading kam, fatherid, motherid
    if (pid2Idx.find(fd[1]) == pid2Idx.end()) {
      pid2Idx[fd[1]] = -1;
      const int val = pid2Idx.size() - 1;
      pid2Idx[fd[1]] = val;
      indv.push_back(fd[1]);
      sex.push_back(atoi(fd[4].c_str()));
      pheno.push_back(atof(fd[5].c_str()));
    } else {
      fprintf(stderr, "duplicated person id [ %s ], ignore!\n", fd[1].c_str());
      exit(1);
    }
  }
  delete lr;
  data.resize(getNumSample());

  variantIdx = 0;
  fprintf(stderr,
          "Finished loading %s.{kam,kim,ked}, %zu markers, %zu samples\n",
          fnPrefix.c_str(), snp2Idx.size(), indv.size());

  return 0;
}

KGGInputFile::~KGGInputFile() {
  delete (this->fpKed);
  delete (this->fpKim);
  delete (this->fpKam);
}

bool KGGInputFile::readRecord() {
  // read a record
  if (variantIdx >= getNumMarker()) {
    return false;
  }
  if (this->fpKed->isEof()) {
    return false;
  }

  const int numAllele = alt[variantIdx].size() + 1;
  const int sampleSize = getNumSample();
  int bytes;
  if (phased) {
    switch (numAllele) {
      case 2:
        bits = 3;
        break;
      case 3:
        bits = 4;
        break;
      case 4:
        bits = 5;
        break;
      default:
        bits = ceil(log(numAllele * numAllele + 1) / log(2));
        break;
    }
  } else {  // unphased, genotype
    switch (numAllele) {
      case 2:
        bits = 2;
        break;
      case 3:
        bits = 3;
        break;
      case 4:
        bits = 4;
        break;
      case 5:
        bits = ceil(log(((numAllele + 1) * numAllele / 2) + 1) / log(2));
        break;
    }
  }
  bytes = (bits * sampleSize + 8) / 8;
  if (phased) {
    buildPhasedTable(numAllele);
  } else {
    buildUnphasedTable(numAllele);
  }

  buffer.resize(bytes);
  this->fpKed->read(buffer.data(), bytes);

  for (int i = 0; i < sampleSize; ++i) {
    data[i] = 0;
  }
  int block = 0;
  unsigned char mask = 1 << 7;
  for (int j = 0; j < bits; ++j) {
    for (int i = 0; i < sampleSize; ++i) {
      data[i] <<= 1;
      data[i] |= (buffer[block] & mask) ? 1 : 0;
      mask >>= 1;
      if (mask == 0) {
        mask = 1 << 7;
        ++block;
        assert(block < bytes);
      }
    }
  }
  return true;
}

//////////////////////////////////////////////////
// Sample inclusion/exclusion
void KGGInputFile::setPeopleMask(const std::string& s, bool b) {
  for (size_t i = 0; i != indv.size(); ++i) {
    if (indv[i] == s) {
      sampleMask[i] = b;
    }
  }
  buildEffectiveIndex();
}
void KGGInputFile::includePeople(const std::string& s) {
  setPeopleMask(s, false);
}

void KGGInputFile::includePeople(const std::vector<std::string>& v) {
  for (size_t i = 0; i != v.size(); ++i) {
    includePeople(v[i].c_str());
  }
}
void KGGInputFile::setPeopleMaskFromFile(const char* fn, bool b) {
  if (!fn || strlen(fn) == 0) {
    return;
  }
  LineReader lr(fn);
  std::vector<std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")) {
    for (unsigned int i = 0; i < fd.size(); i++) {
      setPeopleMask(fd[i].c_str(), b);
    }
  }
  buildEffectiveIndex();
}
void KGGInputFile::includePeopleFromFile(const char* fn) {
  setPeopleMaskFromFile(fn, false);
}
void KGGInputFile::includeAllPeople() {
  std::fill(sampleMask.begin(), sampleMask.end(), false);
  buildEffectiveIndex();
}

void KGGInputFile::excludePeople(const std::string& s) {
  setPeopleMask(s, true);
}
void KGGInputFile::excludePeople(const std::vector<std::string>& v) {
  for (size_t i = 0; i != v.size(); ++i) {
    excludePeople(v[i]);
  }
}
void KGGInputFile::excludePeopleFromFile(const char* fn) {
  setPeopleMaskFromFile(fn, true);
}
void KGGInputFile::excludeAllPeople() {
  std::fill(sampleMask.begin(), sampleMask.end(), true);
  buildEffectiveIndex();
}

//////////////////////////////////////////////////
// Adjust range collections
#if 0
void KGGInputFile::enableAutoMerge() { warnUnsupported("enableAutoMerge"); }
void KGGInputFile::disableAutoMerge() { warnUnsupported("disableAutoMerge"); }
// void clearRange();
void KGGInputFile::setRangeFile(const char* fn) {
  warnUnsupported("setRangeFile");
}
// @param l is a string of range(s)
void KGGInputFile::setRange(const char* chrom, int begin, int end) {
  warnUnsupported("setRange");
}
void KGGInputFile::setRange(const RangeList& rl) {
  warnUnsupported("setRange");
}
void KGGInputFile::setRangeList(const std::string& l) {
  warnUnsupported("setRangeList");
}
// this function the entry point for all function add/change region list
void KGGInputFile::setRangeList(const RangeList& rl) {
  warnUnsupported("setRangeList");
}

void KGGInputFile::setRangeMode() { warnUnsupported("setRangeMode"); }
#endif
int KGGInputFile::setSiteFile(const std::string& fn) {
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

int KGGInputFile::getNumEffectiveSample() const {
  size_t ret = 0;
  for (size_t i = 0; i != sampleMask.size(); ++i) {
    if (sampleMask[i]) continue;
    ret++;
  }
  return ret;
}

void KGGInputFile::getIncludedSampleName(std::vector<std::string>* p) const {
  if (!p) return;
  p->clear();
  for (size_t i = 0; i != sampleMask.size(); ++i) {
    if (sampleMask[i]) continue;
    p->push_back(getSampleName()[i]);
  }
}

void KGGInputFile::buildEffectiveIndex() {
  effectiveIndex.resize(0);
  const size_t N = getNumSample();
  for (size_t i = 0; i != N; ++i) {
    if (sampleMask[i]) continue;
    effectiveIndex.push_back(i);
  }
}

int KGGInputFile::getEffectiveIndex(int idx) const {
  return this->effectiveIndex[idx];
}

int KGGInputFile::getGenotype(int indvIdx) {
  const int nAllele = alt[variantIdx].size() + 1;

  if (phased) {
    return phasedTable[nAllele][data[indvIdx]].x[0] +
           phasedTable[nAllele][data[indvIdx]].x[1];
  } else {
    return unphasedTable[nAllele][data[indvIdx]].x[0] +
           unphasedTable[nAllele][data[indvIdx]].x[1];
  }
}

void KGGInputFile::getAllele(int indvIdx, int* a1, int* a2) {
  const int nAllele = alt[variantIdx].size() + 1;

  if (phased) {
    *a1 = phasedTable[nAllele][data[indvIdx]].x[0];
    *a2 = phasedTable[nAllele][data[indvIdx]].x[1];
  } else {
    *a1 = unphasedTable[nAllele][data[indvIdx]].x[0];
    *a2 = unphasedTable[nAllele][data[indvIdx]].x[1];
  }
}

void KGGInputFile::buildUnphasedTable(int allele) {
  if (unphasedTable.count(allele)) {
    return;
  }
  std::map<char, TwoChar>& m = unphasedTable[allele];
  int val = 0;
  for (int i = 0; i < allele; ++i) {
    for (int j = 0; j <= i; ++j) {
      m[val].x[0] = j;
      m[val].x[1] = i;
      val++;
    }
  }
  m[val].x[0] = -9;
  m[val].x[1] = -9;
}

void KGGInputFile::buildPhasedTable(int allele) {
  if (phasedTable.count(allele)) {
    return;
  }
  std::map<char, TwoChar>& m = phasedTable[allele];
  int val = 0;
  for (int i = 0; i < allele; ++i) {
    for (int j = 0; j < allele; ++j) {
      m[val].x[0] = j;
      m[val].x[1] = i;
      val++;
    }
  }
  m[val].x[0] = -9;
  m[val].x[1] = -9;
}

void KGGInputFile::warnUnsupported(const char* tag) {
  fprintf(stderr, "Please remove unsupported features related to %s", tag);
}
