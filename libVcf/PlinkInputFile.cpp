#include "PlinkInputFile.h"
#include "base/SimpleMatrix.h"

#ifndef UNUSED
#define UNUSED(x) (void)(x)
#endif

#define ROUND_UP_TO_4X(x) (((x) + 3) & ~0x03)

// @param m: people by marker matrix
int PlinkInputFile::readIntoMatrix(SimpleMatrix* mat) const {
  assert(mat);

  // read bed
  int numPeople = this->getNumIndv();
  int numMarker = this->getNumMarker();
  fprintf(stderr, "%d people, %d marker\n", numPeople, numMarker);
  if (snpMajorMode) {
    // unsigned char mask[] = { 0x3, 0xc, 0x30, 0xc0 }; //0b11, 0b1100,
    // 0b110000, 0b11000000
    unsigned char mask = 0x3;  // 0b0000011
    unsigned char c;
    (*mat).resize(numPeople, numMarker);
    for (int m = 0; m < numMarker; m++) {
      for (int p = 0; p < numPeople; p++) {
        int offset = p & (4 - 1);
        if (offset == 0) {
          int ret = fread(&c, sizeof(unsigned char), 1, fpBed);
          assert(ret == 1);
        }
        unsigned char geno = (c >> (offset << 1)) & mask;
        switch (geno) {
          case HOM_REF:
            (*mat)[p][m] = 0;
            break;
          case HET:
            (*mat)[p][m] = 1;
            break;
          case HOM_ALT:
            (*mat)[p][m] = 2;
            break;
          case MISSING:
            (*mat)[p][m] = -9;
            break;
          default:
            REPORT("Read PLINK genotype error!\n");
            break;
        };
        //
        // unsigned char temp =  (c >> (offset << 1) ) & 3 ;
        // printf("m=%d, p=%d, offset=%d, c=%d, geno=%d, geno = %d, temp=%d\n",
        // m,p,offset,c, geno, (int)(*mat)[p][m],temp);
      }
    }
  } else {  // Indv_Major_Mode
    unsigned char mask[] = {0x3, 0xc, 0x30,
                            0xc0};  // 0b11, 0b1100, 0b110000, 0b11000000
    unsigned char c;
    (*mat).resize(numPeople, numMarker);
    for (int p = 0; p < numPeople; p++) {
      for (int m = 0; m < numMarker; m++) {
        int offset = m & (4 - 1);
        if (offset == 0) {
          int ret = fread(&c, sizeof(unsigned char), 1, fpBed);
          assert(ret == 1);
        }
        unsigned char geno = (c & mask[offset]) >> (offset << 1);
        switch (geno) {
          case HOM_REF:
            (*mat)[m][p] = 0;
            break;
          case HET:
            (*mat)[m][p] = 1;
            break;
          case HOM_ALT:
            (*mat)[m][p] = 2;
            break;
          case MISSING:
            (*mat)[m][p] = -9;
            break;
          default:
            REPORT("Read PLINK genotype error!\n");
            break;
        }
      }
    }
  }
  return this->getNumMarker() * this->getNumIndv();
};

int PlinkInputFile::readIntoMatrix(
    SimpleMatrix* mat, std::vector<std::string>* peopleNames,
    std::vector<std::string>* markerNames) const {
  assert(mat);

  // read bed
  int numPeople = getNumIndv();
  int numMarker = getNumMarker();

  mat->resize(peopleNames == NULL ? numPeople : peopleNames->size(),
              markerNames == NULL ? numMarker : markerNames->size());

  // get people index
  std::vector<int> peopleIdx;
  if (peopleNames == NULL || peopleNames->size() == 0) {
    // all peoples
    peopleIdx.resize(pid2Idx.size());
    for (unsigned int i = 0; i < pid2Idx.size(); i++) peopleIdx[i] = (i);
  } else {
    for (unsigned int i = 0; i < peopleNames->size(); i++) {
      std::map<std::string, int>::const_iterator iter =
          pid2Idx.find((*peopleNames)[i]);
      if (iter != pid2Idx.end()) {
        peopleIdx.push_back(iter->second);
      } else {
        fprintf(stderr, "Skip reading the non-exist sample [ %s ]\n",
                (*peopleNames)[i].c_str());
      }
    }
  }

  // get marker index
  std::vector<int> markerIdx;
  if (markerNames == NULL || markerNames->size() == 0) {
    // all markers
    markerIdx.resize(snp.size());
    for (unsigned int i = 0; i < snp.size(); i++) markerIdx[i] = (i);
  } else {
    for (unsigned int i = 0; i < markerNames->size(); i++) {
      std::map<std::string, int>::const_iterator iter =
          snp2Idx.find((*markerNames)[i]);
      if (iter != snp2Idx.end()) {
        markerIdx.push_back(iter->second);
      } else {
        fprintf(stderr, "Skip reading the non-exist marker [ %s ]\n",
                (*markerNames)[i].c_str());
      }
    }
  }

  int peopleToRead = peopleIdx.size();
  int markerToRead = markerIdx.size();

  unsigned char mask[] = {0x3, 0xc, 0x30,
                          0xc0};  // 0b11, 0b1100, 0b110000, 0b11000000
  (*mat).resize(peopleToRead, markerToRead);
  if (snpMajorMode) {
    for (int p = 0; p < peopleToRead; p++) {
      for (int m = 0; m < markerToRead; m++) {
        // get file position
        int pos = 3 + (numPeople / 4 + 1) * markerIdx[m] + peopleIdx[p] / 4;
        int offset = peopleIdx[p] % 4;
        unsigned char c;
        int tmp = fseek(this->fpBed, pos, SEEK_SET);
        UNUSED(tmp);
        int nRead = fread(&c, sizeof(unsigned char), 1, fpBed);
        if (nRead != 1) return -1;
        unsigned char geno = (c & mask[offset]) >> (offset << 1);
        switch (geno) {
          case HOM_REF:
            (*mat)[p][m] = 0;
            break;
          case HET:
            (*mat)[p][m] = 1;
            break;
          case HOM_ALT:
            (*mat)[p][m] = 2;
            break;
          case MISSING:
            (*mat)[p][m] = -9;
            break;
          default:
            REPORT("Read PLINK genotype error!\n");
            break;
        };
      }
    }
  } else {  // Indv_Major_Mode
    for (int p = 0; p < numPeople; p++) {
      for (int m = 0; m < numMarker; m++) {
        // get file position
        int pos = 3 + (numMarker / 4 + 1) * peopleIdx[p] + markerIdx[m] / 4;
        int offset = markerIdx[m] % 4;
        unsigned char c;
        fseek(this->fpBed, pos, SEEK_SET);
        int nRead = fread(&c, sizeof(unsigned char), 1, fpBed);
        if (nRead != 1) return -1;

        unsigned char geno = (c & mask[offset]) >> (offset << 1);
        switch (geno) {
          case HOM_REF:
            (*mat)[m][p] = 0;
            break;
          case HET:
            (*mat)[m][p] = 1;
            break;
          case HOM_ALT:
            (*mat)[m][p] = 2;
            break;
          case MISSING:
            (*mat)[m][p] = -9;
            break;
          default:
            REPORT("Read PLINK genotype error!\n");
            break;
        };
      }
    }
  }
  return this->getNumMarker() * this->getNumIndv();
}

int PlinkInputFile::calculateMAF(std::vector<double>* maf) {
  assert(this->fpBed);

  int N = getNumIndv();
  int M = getNumMarker();
  maf->resize(M);

  long fpPrevPosition = ftell(this->fpBed);
  fseek(this->fpBed, 3, SEEK_SET);
  if (snpMajorMode) {
    int stride = ROUND_UP_TO_4X(N) / 4;
    unsigned char* buf = new unsigned char[stride];
    for (int i = 0; i < M; ++i) {
      fread(buf, sizeof(unsigned char), stride, this->fpBed);
      int numAllele = 0;
      int numMiss = 0;
      for (int j = 0; j < N; ++j) {
        switch (extract2Bit(buf[j >> 2], j & 0x03)) {
          case HET:
            numAllele++;
            break;
          case HOM_ALT:
            numAllele += 2;
            break;
          case MISSING:
            numMiss++;
        }
      }
      if (N == numMiss) {
        (*maf)[i] = 0;
      } else {
        (*maf)[i] = 0.5 * (numAllele) / (N - numMiss);
      }
    }

    delete[] buf;
  } else {  // individual-major mode
    int stride = ROUND_UP_TO_4X(M) / 4;
    unsigned char* buf = new unsigned char[stride];
    std::vector<int> numAllele(N);
    std::vector<int> numMiss(N);
    for (int i = 0; i < N; ++i) {
      fread(buf, sizeof(unsigned char), stride, this->fpBed);
      for (int j = 0; j < M; ++j) {
        switch (extract2Bit(buf[j >> 2], j & 0x03)) {
          case HET:
            numAllele[j] += 1;
            break;
          case HOM_ALT:
            numAllele[j] += 2;
            break;
          case MISSING:
            numMiss[i] += 1;
        }
      }
    }
    for (int j = 0; j < M; ++j) {
      if (N == numMiss[j]) {
        (*maf)[j] = 0.;
      } else {
        (*maf)[j] = 0.5 * numAllele[j] / (N - numMiss[j]);
      }
    }
    delete[] buf;
  }
  fseek(this->fpBed, fpPrevPosition, SEEK_SET);
  return 0;
}

int PlinkInputFile::calculateMissing(std::vector<double>* imiss,
                                     std::vector<double>* lmiss) {
  assert(this->fpBed);

  int N = getNumIndv();
  int M = getNumMarker();
  imiss->resize(N);
  lmiss->resize(M);

  long fpPrevPosition = ftell(this->fpBed);
  fseek(this->fpBed, 3, SEEK_SET);
  if (snpMajorMode) {
    int stride = ROUND_UP_TO_4X(N) / 4;
    unsigned char* buf = new unsigned char[stride];
    for (int i = 0; i < M; ++i) {
      fread(buf, sizeof(unsigned char), stride, this->fpBed);
      for (int j = 0; j < N; ++j) {
        if (extract2Bit(buf[j >> 2], (j & 0x03)) == MISSING) {
          (*imiss)[j]++;
          (*lmiss)[i]++;
        }
      }
    }

    delete[] buf;
  } else {  // individual-major mode
    int stride = ROUND_UP_TO_4X(M) / 4;
    unsigned char* buf = new unsigned char[stride];
    std::vector<int> numAllele(N);
    std::vector<int> numMiss(N);
    for (int i = 0; i < N; ++i) {
      fread(buf, sizeof(unsigned char), stride, this->fpBed);
      for (int j = 0; j < M; ++j) {
        if ((buf[j >> 2] >> ((j & ~0x03) >> 1)) == MISSING) {
          (*imiss)[i]++;
          (*lmiss)[j]++;
        }
      }
    }
    delete[] buf;
  }
  fseek(this->fpBed, fpPrevPosition, SEEK_SET);

  for (int i = 0; i < N; ++i) {
    (*imiss)[i] /= M;
  }
  for (int i = 0; i < M; ++i) {
    (*lmiss)[i] /= N;
  }
  return 0;
}

int PlinkInputFile::readBED(unsigned char* buf, size_t n) {
  size_t nRead = 0;
  while (nRead < n) {
    nRead += fread(buf + nRead, sizeof(unsigned char), n, this->fpBed);
  }
  return nRead;
}

unsigned char PlinkInputFile::get2BitGenotype(int sample, int marker) {
  assert(0 <= sample && sample < getNumSample());
  assert(0 <= marker && marker < getNumMarker());

  const int M = getNumMarker();
  const int N = getNumSample();

  unsigned char c;
  if (snpMajorMode) {
    int stride = ROUND_UP_TO_4X(N) / 4;
    fseek(this->fpBed, 3 + stride * marker + (sample >> 2), SEEK_SET);
    fread(&c, sizeof(unsigned char), 1, this->fpBed);
    return extract2Bit(c, sample & 0x03);
  } else {
    int stride = ROUND_UP_TO_4X(M) / 4;
    fseek(this->fpBed, 3 + stride * sample + (marker >> 2), SEEK_SET);
    fread(&c, sizeof(unsigned char), 1, this->fpBed);
    return extract2Bit(c, sample & 0x03);
  }
}

unsigned char PlinkInputFile::extract2Bit(unsigned char g, int i) {
  assert(0 <= i && i <= 3);
  return (g >> (i << 1)) & 3;
}
