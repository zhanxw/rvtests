#include "libBgen/BGenFile.h"

#include <cassert>
#include "base/RangeList.h"
#include "base/TypeConversion.h"
#include "libBgen/BitReader.h"

BGenFile::BGenFile(const std::string& fn) : var(this->N) {
  this->bgenFileName = fn;
  this->mode = BGEN_LINE_MODE;

  int nRead;
  fp = fopen(fn.c_str(), "rb");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open file [ %s ]\n", fn.c_str());
    exit(1);
  }
  // first 4 bytes
  parseUint32(fp, &offset);
#ifdef DEBUG
  printf("nRead = %d, offset = %d\n", nRead, (int)offset);
#endif

  // header block
  parseUint32(fp, &LH);
  parseUint32(fp, &M);
  parseUint32(fp, &N);
#ifdef DEBUG
  printf("nRead = %d, LH = %d\n", nRead, (int)LH);
  printf("nRead = %d, M = %d\n", nRead, (int)M);
  printf("nRead = %d, N = %d\n", nRead, (int)N);
#endif
  sampleMask.resize(N);
  std::fill(sampleMask.begin(), sampleMask.end(), false);

  nRead = fread(&magic, sizeof(magic[0]), 4, fp);
#ifdef DEBUG
  printf("nRead = %d, magic = %c%c%c%c\n", nRead, magic[0], magic[1], magic[2],
         magic[3]);
#endif
  assert(nRead == 4);
  if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' &&
        magic[3] == 'n')) {
    // assert(false);
    fprintf(stderr,
            "bgen magic number [ %c%c%c%c ] does not match ['bgen']! \n",
            magic[0], magic[1], magic[2], magic[3]);
  }

  int freeDataLen = LH - 20;
  freeData.resize(freeDataLen);
  nRead = fread(freeData.data(), sizeof(freeData[0]), freeDataLen, fp);
#ifdef DEBUG
  printf("nRead = %d, freeData = ", nRead);
  for (int i = 0; i < freeDataLen; ++i) {
    printf("%c", freeData[i]);
  }
  printf("\n");
#endif

  parseUint32(fp, &flag);
#ifdef DEBUG
  printf("nRead = %d, flag = %zu\n", nRead, flag);
#endif

  snpCompression = static_cast<SNP_COMPRESSION>(flag & 3);
  layout = static_cast<LAYOUT>((flag >> 2) & 0xf);
  flagSampleIdentifier = static_cast<SAMPLE_IDENTIFIER>(flag >> 31);
#ifdef DEBUG
  printf("snpCompression = %d\t", (int)snpCompression);
  switch (snpCompression) {
    case NO_COMPRESSION:
      printf("No compression\n");
      break;
    case GZIP:
      printf("GZIP\n");
      break;
    case ZSTD:
      printf("ZSTD\n");
      break;
    default:
      printf("Wrong value!\n");
  }

  printf("layout= %d\n", (int)layout);
  printf("flagSampleIdentifier = %d %s\n", (int)flagSampleIdentifier,
         (int)flagSampleIdentifier == 1 ? "Have sample id"
                                        : "Do not have sample id");
#endif
  // sample identifier block
  if (flagSampleIdentifier == HAS_SAMPLE_IDENTIFIER) {
    uint32_t LSI;
    parseUint32(fp, &LSI);
#ifdef DEBUG
    printf("nRead = %d, LSI = %d\n", nRead, (int)LSI);
#endif
    if (!(LSI + LH <= offset)) {
      assert(false);
    }

    uint32_t N2;
    parseUint32(fp, &N2);
#ifdef DEBUG
    printf("nRead = %d, N2 = %d\n", nRead, (int)N2);
#endif
    assert(N2 == N);
    sampleIdentifier.resize(N);
    uint16_t siLen;
    for (uint32_t i = 0; i < N; ++i) {
      nRead = fread(&siLen, sizeof(siLen), 1, fp);
      assert(nRead == 1);

      sampleIdentifier[i].resize(siLen);
      nRead = fread(&sampleIdentifier[i][0], sizeof(sampleIdentifier[i][0]),
                    siLen, fp);
#ifdef DEBUG
      printf("nRead = %d, si[%d] = %d, %s\n", nRead, i, (int)siLen,
             sampleIdentifier[i].c_str());
#endif
    }
  }
  // advance fp to variant data blocks
  // NOTE: the assertion below may not work
  // sometimes, unspecified bytes may exists between header block and variant
  // data blocks.
  // long variantOffset = ftell(fp);
  // if (variantOffset != (long)offset + 4) {
  // }
  fseek(fp, offset + 4, SEEK_SET);

  fileSize = getFileSize(fn);
}

bool BGenFile::readRecord() {
  if (mode == BGEN_RANGE_MODE) {
    int file_pos, bytes;
    if (index.next(&file_pos, &bytes)) {
      fseek(fp, file_pos, SEEK_SET);
    } else {
      return false;
    }
  }

  switch (layout) {
    case LAYOUT1:
      return parseLayout1();
    case LAYOUT2:
      return parseLayout2();
    default:
      assert(false);
  }
  return false;
}

bool BGenFile::parseLayout1() {
  if (isFileEnd(fp)) {
    return false;
  }
  // variant identifying data
  uint32_t NinRow;
  int nRead = fread(&NinRow, sizeof(NinRow), 1, fp);
#ifdef DEBUG
  printf("nRead = %d, NinRow = %d\n", nRead, (int)NinRow);
#endif
  assert(nRead == 1);

  parseString(fp, 2, &var.varid);
  parseString(fp, 2, &var.rsid);
  parseString(fp, 2, &var.chrom);
  parseUint32(fp, &var.pos);

  //      const uint16_t K = 2;  // omitted in layout 1 by specification
  var.K = 2;
  var.alleles.resize(var.K);
  for (int i = 0; i < var.K; ++i) {
    parseString(fp, 4, &var.alleles[i]);
  }

  // genotype data block
  assert(snpCompression ==
         GZIP);  // do not deal with no-compression case for now
  nRead = fread(&C, sizeof(C), 1, fp);
  assert(nRead == 1);
#ifdef DEBUG
  printf("C = %zu\n", C);
#endif

  // std::vector<uint8_t> buf(NinRow * 6);
  D = NinRow * 6;
  buf.resize(D);
  // std::vector<uint8_t> bufCompress(C);
  compressedBuf.resize(C);
  nRead = fread(compressedBuf.data(), sizeof(uint8_t), C, fp);
  assert(nRead == (int)C);
  // unsigned long bufLen = NinRow * 6;
  unsigned long decompressedByte = NinRow * 6;
  int zlibStatus =
      uncompress(buf.data(), &decompressedByte, compressedBuf.data(), C);
  assert(zlibStatus == Z_OK);

  // parse probility
  var.missing.resize(N);
  var.ploidy.resize(N);
  var.isPhased = false;
  var.index.resize(N + 1);
  var.prob.resize(N * 3);
  float p[3];
  for (size_t i = 0; i < NinRow; ++i) {
    var.ploidy[i] = 2;
    var.index[i] = 3 * i;
    uint16_t* probPtr = (uint16_t*)(buf.data() + i * 6);
    p[0] = (float)(*probPtr) / 32768;
    ++probPtr;
    p[1] = (float)(*probPtr) / 32768;
    ++probPtr;
    p[2] = (float)(*probPtr) / 32768;
#ifdef DEBUG
    printf("prob = (%g, %g, %g)\n", i, p[0], p[1], p[2]);
#endif

    if (p[0] == 0 && p[1] == 0 && p[2] == 0) {
      var.missing[i] = true;
    } else {
      var.missing[i] = false;
    }
    var.prob[i * 3] = p[0];
    var.prob[i * 3 + 1] = p[1];
    var.prob[i * 3 + 2] = p[2];
  }
  var.index.push_back(3 * N);
#ifdef DEBUG
  printf("feof = %d\n", feof(fp));
#endif

  return true;
}

bool BGenFile::parseLayout2() {
  if (isFileEnd(fp)) {
    return false;
  }
  // variant identifying data
  parseString(fp, 2, &var.varid);
  parseString(fp, 2, &var.rsid);
  parseString(fp, 2, &var.chrom);
  parseUint32(fp, &var.pos);

  parseUint16(fp, &var.K);
#ifdef DEBUG
  printf("K = %zu\n", var.K);
#endif

  var.alleles.resize(var.K);
  for (int i = 0; i < var.K; ++i) {
    parseString(fp, 4, &var.alleles[i]);
  }

  uint32_t C;
  parseUint32(fp, &C);
#ifdef DEBUG
  printf("C = %zu\n", C);
#endif

  uint32_t D;
  if (snpCompression == NO_COMPRESSION) {
    D = C;
  } else {
    parseUint32(fp, &D);
#ifdef DEBUG
    printf("D = %zu\n", D);
#endif
  }
  // genotype data block
  std::vector<uint8_t> compressedBuf;
  if (snpCompression == NO_COMPRESSION) {
    compressedBuf.resize(C);
  } else {
    compressedBuf.resize(C - 4);
  }
  std::vector<uint8_t> buf(D);
  size_t nRead;
  if (snpCompression == GZIP) {
    nRead = fread(compressedBuf.data(), sizeof(uint8_t), C - 4, fp);
    assert(nRead == C - 4);
    unsigned long bufLen = D;
    nRead = uncompress(buf.data(), &bufLen, compressedBuf.data(), C);
    assert(nRead == Z_OK);
  } else if (snpCompression == ZSTD) {
    nRead = fread(compressedBuf.data(), sizeof(uint8_t), C - 4, fp);
    assert(nRead == C - 4);
    unsigned long bufLen = D;
    // TODO: create ZSTD context to save time
    size_t ret =
        ZSTD_decompress(buf.data(), bufLen, compressedBuf.data(), C - 4);
    if (ret > bufLen) {
      if (ZSTD_isError(ret)) {
#ifdef DEBUG
        printf("ZSTD error: %s\n", ZSTD_getErrorName(ret));
#endif
      }
    }
    assert(ret == bufLen);
  } else if (snpCompression == NO_COMPRESSION) {
    // nRead = fread(compressedBuf.data(), sizeof(uint8_t), C , fp);
    // assert(nRead == C);
    nRead = fread(buf.data(), sizeof(uint8_t), D, fp);
    assert(nRead == D);
  }

  const uint32_t nIndv = *(uint32_t*)buf.data();
  assert(nIndv == N);
  const uint16_t nAllele = *(uint16_t*)(buf.data() + 4);
  assert(nAllele >= 0);
  const uint8_t nMinAllele = *(uint8_t*)(buf.data() + 6);
  const uint8_t nMaxAllele = *(uint8_t*)(buf.data() + 7);
  assert(0 <= nMinAllele && nMinAllele <= 63);
  assert(0 <= nMaxAllele && nMaxAllele <= 63);
  const uint8_t* ploidityAndMissing = (uint8_t*)(buf.data() + 8);
  const uint8_t isPhased = *(uint8_t*)(buf.data() + 8 + N);
  var.isPhased = isPhased != 0;
  var.B = *(uint8_t*)(buf.data() + 8 + N + 1);  // bits
  assert(1 <= var.B && var.B <= 32);
#ifdef DEBUG
  int totalBit = 0;
  int cumBit = 0;
#endif

  var.missing.resize(N);
  var.ploidy.resize(N);
  var.index.reserve(N + 1);
  var.index.resize(0);
  var.prob.resize(0);
  BitReader br(buf.data() + 8 + N + 2, (D - 8 - N - 2), var.B);
  for (uint32_t i = 0; i < nIndv; ++i) {
    var.index.push_back(var.prob.size());
    const uint8_t ploidy = ploidityAndMissing[i] & 0x3f;
    const int Z = ploidy;
    const bool missing = (ploidityAndMissing[i] & 0x80) != 0;
    var.ploidy[i] = ploidy;
    var.missing[i] = missing;
#ifdef DEBUG
    printf("ploidy = %d, missing = %s\n", Z, missing ? "true" : "false");
#endif

    // ploidy = Z, allele = K
    if (var.isPhased) {
// total Z * (K- 1) * B bits bits used
#ifdef DEBUG
      totalBit = Z * (var.K - 1) * B;
      printf("Total phased bits = %d\n", totalBit);
#endif

      for (int j = 0; j < ploidy; ++j) {
        float remainProb = 1.0;
        for (int k = 0; k < var.K - 1; ++k) {
          float p = br.next();
          var.prob.push_back(p);
          remainProb -= p;
        }
        var.prob.push_back(remainProb);
      }
    } else {
      // total ( ( (Z+K-1) choose (K-1) ) -1 ) * B bits used
      const int nCombination = choose((Z + var.K - 1), (var.K - 1));
#ifdef DEBUG
      totalBit = (nCombination - 1) * B;
      printf("Total unphased bits = %d\n", totalBit);
#endif
      float remainProb = 1.0;
      for (int j = 0; j < nCombination - 1; ++j) {
        float p = br.next();
        var.prob.push_back(p);
        remainProb -= p;
      }
      var.prob.push_back(remainProb);
    }
#ifdef DEBUG
    cumBit += totalBit;
#endif
  }
  var.index.push_back(var.prob.size());

#ifdef DEBUG
  printf("CumBit = %d\n", cumBit);
  printf("Total chunk = %d\n", 10 + N + cumBit / 8);
  printf("feof = %d\n", feof(fp));
#endif

  return true;
}

void BGenFile::parseString(FILE* fp, int lenByte, std::string* out) {
  if (lenByte == 2) {
    uint16_t siLen;
    int nRead = fread(&siLen, sizeof(siLen), 1, fp);
    assert(nRead == 1);

    (*out).resize(siLen);
    // (*out)[siLen] = '\0';
    nRead = fread(&(*out)[0], sizeof((*out)[0]), siLen, fp);
#ifdef DEBUG
    printf("nRead = %d, read = %s\n", nRead, (*out).c_str());
#endif

  } else if (lenByte == 4) {
    uint32_t siLen;
    int nRead = fread(&siLen, sizeof(siLen), 1, fp);
    assert(nRead == 1);

    (*out).resize(siLen);
    nRead = fread(&(*out)[0], sizeof((*out)[0]), siLen, fp);
#ifdef DEBUG
    printf("nRead = %d, read = %s\n", nRead, (*out).c_str());
#endif

  } else {
    assert(false);
  }
}

void BGenFile::parseUint32(FILE* fp, uint32_t* value) {
  int nRead = fread(value, sizeof(value[0]), 1, fp);
  assert(nRead == 1);
}
void BGenFile::parseUint16(FILE* fp, uint16_t* value) {
  int nRead = fread(value, sizeof(value[0]), 1, fp);
  assert(nRead == 1);
}
// @return choose m out of n elements
int BGenFile::choose(int n, int m) {
  if (m == 1) {
    return n;
  }
  if (n == 1) {
    return 1;
  }
  int ret = 1;
  for (int i = 0; i < m; ++i) {
    ret *= (n - i);
  }
  for (int i = 0; i < m; ++i) {
    ret /= (i + 1);
  }
  return ret;
}

bool BGenFile::isFileEnd(FILE* fp) { return (fileSize - ftell(fp)) == 0; }
long BGenFile::getFileSize(const std::string& fn) {
  FILE* fp = fopen(fn.c_str(), "rb");
  fseek(fp, 0, SEEK_END);
  long ret = ftell(fp);
  fclose(fp);
  return ret;
}

//////////////////////////////////////////////////
// Sample inclusion/exclusion
void BGenFile::setPeopleMask(const std::string& s, bool b) {
  for (size_t i = 0; i != sampleIdentifier.size(); ++i) {
    if (sampleIdentifier[i] == s) {
      sampleMask[i] = b;
    }
  }
  buildEffectiveIndex();
}
void BGenFile::includePeople(const std::string& s) { setPeopleMask(s, false); }

void BGenFile::includePeople(const std::vector<std::string>& v) {
  for (size_t i = 0; i != v.size(); ++i) {
    includePeople(v[i].c_str());
  }
}
void BGenFile::setPeopleMaskFromFile(const char* fn, bool b) {
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
void BGenFile::includePeopleFromFile(const char* fn) {
  setPeopleMaskFromFile(fn, false);
}
void BGenFile::includeAllPeople() {
  std::fill(sampleMask.begin(), sampleMask.end(), false);
  buildEffectiveIndex();
}

void BGenFile::excludePeople(const std::string& s) { setPeopleMask(s, true); }
void BGenFile::excludePeople(const std::vector<std::string>& v) {
  for (size_t i = 0; i != v.size(); ++i) {
    excludePeople(v[i]);
  }
}
void BGenFile::excludePeopleFromFile(const char* fn) {
  setPeopleMaskFromFile(fn, true);
}
void BGenFile::excludeAllPeople() {
  std::fill(sampleMask.begin(), sampleMask.end(), true);
  buildEffectiveIndex();
}

//////////////////////////////////////////////////
// Adjust range collections
void BGenFile::enableAutoMerge() { this->autoMergeRange = true; }
void BGenFile::disableAutoMerge() { this->autoMergeRange = false; }
// void clearRange();
void BGenFile::setRangeFile(const char* fn) {
  if (!fn || strlen(fn) == 0) return;
  RangeList r;
  r.addRangeFile(fn);
  this->setRange(r);
}
// @param l is a string of range(s)
void BGenFile::setRange(const char* chrom, int begin, int end) {
  RangeList r;
  r.addRange(chrom, begin, end);
  this->setRange(r);
}
void BGenFile::setRange(const RangeList& rl) { this->setRangeList(rl); }
void BGenFile::setRangeList(const std::string& l) {
  if (l.empty()) return;

  RangeList r;
  r.addRangeList(l);
  this->setRange(r);
}
// this function the entry point for all function add/change region list
void BGenFile::setRangeList(const RangeList& rl) {
  if (rl.size() == 0) return;

  this->setRangeMode();

  RangeList l;
  l.setRange(rl);
  if (this->autoMergeRange) {
    l.sort();
  }

  if (mode == BGEN_RANGE_MODE) {
    index.setRange(rl);
    // this->tabixReader->setRange(l);
  } else {
    fprintf(stderr, "[ERROR] invalid reading mode, quitting...\n");
    abort();
  }
}

void BGenFile::setRangeMode() {
  if (this->index.init(bgenFileName + ".bgi")) {
    fprintf(stderr, "Cannot open BGEN index file [ %s ]!\n",
            (bgenFileName + ".bgi").c_str());
    abort();
  }
  this->mode = BGEN_RANGE_MODE;
}

int BGenFile::setSiteFile(const std::string& fn) {
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

int BGenFile::loadSampleFile(const std::string& fn) {
  if (fn.empty()) return 0;

  LineReader lr(fn);
  std::vector<std::string> fd;
  std::vector<std::string> sn;
  int lineNo = 0;
  while (lr.readLineBySep(&fd, " \t")) {
    ++lineNo;
    if (fd.size() < 3) {
      fprintf(stderr,
              "Error: Line [ %d ] of the input file [ %s ] has less than 3 "
              "columns!\n",
              lineNo, fn.c_str());
      return -1;
    }
    if (lineNo == 1) {
      if (fd[0] != "ID_1" || fd[1] != "ID_2" || fd[2] != "missing") {
        fprintf(stderr,
                "ERROR: The header line does not start with ID_1, ID_2 and "
                "missing!\n");
        return -1;
      }
    } else if (lineNo == 2) {
      if (fd[0] != "0" || fd[1] != "0" || fd[2] != "0") {
        fprintf(stderr,
                "ERROR: The second line (variable type line) does not start "
                "with 0, 0 and 0!\n");
        return -1;
      }
    } else {
      if (fd[0] != fd[1]) {
        fprintf(stderr,
                "WARN: Line [ %d ], ID_1 [ %s ] is different than ID_2 [ %s ], "
                "by default ID_2 column is used.\n",
                lineNo, fd[0].c_str(), fd[1].c_str());
      }
      sn.push_back(fd[1]);
    }
  }
  fprintf(stderr, "INFO: Loaded [ %d ] samples from .sample file [ %s ]\n",
          (int)sn.size(), fn.c_str());
  if (sn.size() != N) {
    fprintf(stderr,
            "ERROR: Sample file has [ %d ] samples, but BGEN file have [ %d ] "
            "samples, skipped loading this sample file\n",
            (int)sn.size(), N);
    return -1;
  } else {
    this->sampleIdentifier = sn;
    return 0;
  }
}

int BGenFile::getNumEffectiveSample() const {
  size_t ret = 0;
  for (size_t i = 0; i != sampleMask.size(); ++i) {
    if (sampleMask[i]) continue;
    ret++;
  }
  return ret;
}

void BGenFile::getIncludedSampleName(std::vector<std::string>* p) const {
  if (!p) return;
  p->clear();
  for (size_t i = 0; i != sampleMask.size(); ++i) {
    if (sampleMask[i]) continue;
    p->push_back(sampleIdentifier[i]);
  }
}

void BGenFile::buildEffectiveIndex() {
  effectiveIndex.resize(0);
  for (size_t i = 0; i != N; ++i) {
    if (sampleMask[i]) continue;
    effectiveIndex.push_back(i);
  }
}

int BGenFile::getEffectiveIndex(int idx) const {
  return this->effectiveIndex[idx];
}

void BGenFile::printInfo() {
  printf("--First 4 bytes--\n");
  printf("\toffset = %d\n", (int)offset);

  // header block
  printf("--Header block--\n");
  printf("\tLH = %d\n", (int)LH);
  printf("\tM = %d\n", (int)M);
  printf("\tN = %d\n", (int)N);

  if (freeData.empty()) {
    printf("\tfreeData = []\n");
  } else {
    printf("\tfreeData = [");
    for (size_t i = 0; i < freeData.size(); ++i) {
      printf("%c", freeData[i]);
    }
    printf("]\n");
  }
  printf("\tflag = %x\n", flag);
  printf("\tsnpCompression = %d ", (int)snpCompression);
  switch (snpCompression) {
    case NO_COMPRESSION:
      printf("(No compression)\n");
      break;
    case GZIP:
      printf("(GZIP)\n");
      break;
    case ZSTD:
      printf("(ZSTD)\n");
      break;
    default:
      printf("Wrong value!\n");
  }

  printf("\tlayout= %d\n", (int)layout);
  printf("\tflagSampleIdentifier = %d %s\n", (int)flagSampleIdentifier,
         (int)flagSampleIdentifier == 1 ? "(Have sample id)"
                                        : "(Do not have sample id)");

  // sample identifier block
  if (flagSampleIdentifier == HAS_SAMPLE_IDENTIFIER) {
    printf("--Sample identifier block--\n");
    // printf("LSI = %d\n", (int)LSI);
    // printf("N2 = %d\n", (int)N2);
    // assert(N2 == N);
    for (uint32_t i = 0; i < N; ++i) {
      printf("\tsi[%d] = %s\n", i, sampleIdentifier[i].c_str());
    }
  }

  // variant data blocks
  printf("--Variant data block--\n");
  for (size_t i = 0; i < M; ++i) {
    if (!readRecord()) {
      printf("\tNo variants presented, file truncated?\n");
      break;
    }

    printf("\t[Variant %d]\n", (int)i);

    printf("\tChromosomal position: %s:%d\n", var.chrom.c_str(), var.pos);
    printf("\tRSID = %s\n", var.rsid.c_str());
    printf("\tVarID = %s\n", var.varid.c_str());
    printf("\tAlleles = %s ", var.alleles[0].c_str());
    for (size_t j = 1; j != var.alleles.size(); ++j) {
      printf("/ %s ", var.alleles[j].c_str());
    }
    printf("\n");

    if (layout == LAYOUT1) {
      printf("\tPolidy = 2\n");
      printf("\tUnphased\n");
      printf("\tAlleles = 2\n");
      printf("\tBitsPerGenotypeProbability = 16\n");  // 2 bytes per genotype
                                                      // probability
      int nMissing = 0;
      for (size_t j = 0; j != N; ++j) {
        if (var.prob[j * 3] == 0.0 && var.prob[j * 3 + 1] == 0.0 &&
            var.prob[j * 3 + 2] == 0.0)
          ++nMissing;
      }
      printf("Missing = %d\t", nMissing);
    } else if (layout == LAYOUT2) {
      int nMissing = 0;
      for (size_t i = 0; i < var.missing.size(); ++i) {
        if (var.missing[i]) ++nMissing;
      }
      std::set<uint8_t> s = makeSet(var.ploidy);
      std::string ss = toString(s, ",");
      printf("\tPolidy = %s\n", ss.c_str());
      printf("\t%s\n", var.isPhased ? "Phased" : "Unphased");
      printf("\tAlleles = %d\n", var.K);
      printf("\tBitsPerGenotypeProbability = %d\n", (int)var.B);
      printf("\tMissing = %d\n", nMissing);
    } else {
      assert(false);
    }
  }
}
