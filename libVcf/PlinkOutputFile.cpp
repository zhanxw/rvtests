#include "libVcf/PlinkOutputFile.h"

#include "base/SimpleMatrix.h"
#include "libVcf/PlinkInputFile.h"
#include "libVcf/VCFHeader.h"
#include "libVcf/VCFRecord.h"

void PlinkOutputFile::init(const char* fnPrefix) {
  std::string prefix = fnPrefix;
  this->fpBed = fopen((prefix + ".bed").c_str(), "wb");
  this->fpBim = fopen((prefix + ".bim").c_str(), "wt");
  this->fpFam = fopen((prefix + ".fam").c_str(), "wt");
  if (!this->fpBed || !this->fpBim || !this->fpFam) {
    REPORT("Cannot create binary PLINK file!");
    abort();
  }
  // write Bed header
  char c;
  // magic number
  c = 0x6c;  // 0b01101100;
  fwrite(&c, sizeof(char), 1, this->fpBed);
  c = 0x1b;  // 0b00011011;
  fwrite(&c, sizeof(char), 1, this->fpBed);
  // snp major mode
  c = 0x01;  // 0b00000001;
  fwrite(&c, sizeof(char), 1, this->fpBed);
}

void PlinkOutputFile::writeHeader(const VCFHeader* h) {
  std::vector<std::string> people;
  h->getPeopleName(&people);
  this->writeFAM(people);
}

int PlinkOutputFile::isMultiAllelic(const char* r) {
  if (strchr(r, ',') == NULL) return 0;

  return -1;
}

int PlinkOutputFile::writeRecord(VCFRecord* r) {
  // write BIM
  if (isMultiAllelic(r->getRef()) || isMultiAllelic(r->getAlt())) {
    fprintf(stdout, "%s:%d Skip with ref = [ %s ] and alt= [ %s ]\n", __FILE__,
            __LINE__, r->getRef(), r->getAlt());
    return -1;
  }

  this->writeBIM(r->getChrom(), r->getID(), 0, r->getPos(), r->getRef(),
                 r->getAlt());

  // write BED
  int GTidx = r->getFormatIndex("GT");
  VCFPeople& people = r->getPeople();
  unsigned char c = 0;
  VCFIndividual* indv;
  int offset;
  for (unsigned int i = 0; i < people.size(); i++) {
    indv = people[i];
    offset = i & (4 - 1);
    if (indv->justGet(GTidx).isHaploid()) {  // 0: index of GT
      int a1 = indv->justGet(GTidx).getAllele1();
      if (a1 == 0)
        setGenotype(&c, offset, HOM_REF);
      else if (a1 == 1)
        setGenotype(&c, offset, HET);
      else
        setGenotype(&c, offset, MISSING);
    } else {
      int a1 = indv->justGet(GTidx).getAllele1();
      int a2 = indv->justGet(GTidx).getAllele2();
      if (a1 == 0) {
        if (a2 == 0) {
          // homo ref: 0b00
        } else if (a2 == 1) {
          setGenotype(&c, offset, HET);  // het: 0b01
        } else {
          setGenotype(&c, offset, MISSING);  // missing 0b10
        }
      } else if (a1 == 1) {
        if (a2 == 0) {
          setGenotype(&c, offset, HET);  // het: 0b01
        } else if (a2 == 1) {
          setGenotype(&c, offset, HOM_ALT);  // hom alt: 0b11
        } else {
          setGenotype(&c, offset, MISSING);  // missing
        }
      } else {
        // NOTE: Plink does not support tri-allelic
        // so have to set genotype as missing.
        setGenotype(&c, offset, MISSING);  // missing
      }
    }
    if (offset == 3) {  // 3: 4 - 1, so every 4 genotype we will flush
      fwrite(&c, sizeof(char), 1, this->fpBed);
      c = 0;
    }
  }
  if (people.size() % 4 != 0) {  // remaining some bits
    fwrite(&c, sizeof(char), 1, this->fpBed);
  }
  return 0;
}

int PlinkOutputFile::writeRecordWithFilter(VCFRecord* r, const double minGD,
                                           const double minGQ) {
  // write BIM
  if (isMultiAllelic(r->getRef()) || isMultiAllelic(r->getAlt())) {
    fprintf(stdout, "%s:%d Skip with ref = [ %s ] and alt= [ %s ]\n", __FILE__,
            __LINE__, r->getRef(), r->getAlt());
    return -1;
  }

  this->writeBIM(r->getChrom(), r->getID(), 0, r->getPos(), r->getRef(),
                 r->getAlt());

  // write BED
  int GTidx = r->getFormatIndex("GT");
  int GDidx = r->getFormatIndex("GD");
  int GQidx = r->getFormatIndex("GQ");
  bool missing = false;

  VCFPeople& people = r->getPeople();
  unsigned char c = 0;
  VCFIndividual* indv;
  int offset;
  for (unsigned int i = 0; i < people.size(); i++) {
    indv = people[i];
    offset = i & (4 - 1);

    missing = false;
    if (minGD > 0 &&
        (GDidx < 0 || (GDidx > 0 && indv->justGet(GDidx).toDouble() < minGD))) {
      missing = true;
    }
    if (missing && minGQ > 0 &&
        (GQidx < 0 || (GQidx > 0 && indv->justGet(GQidx).toDouble() < minGQ))) {
      missing = true;
    }
    if (!missing) {
      if (indv->justGet(GTidx).isHaploid()) {  // 0: index of GT
        int a1 = indv->justGet(GTidx).getAllele1();
        if (a1 == 0)
          setGenotype(&c, offset, HOM_REF);
        else if (a1 == 1)
          setGenotype(&c, offset, HET);
        else
          setGenotype(&c, offset, MISSING);
      } else {
        int a1 = indv->justGet(GTidx).getAllele1();
        int a2 = indv->justGet(GTidx).getAllele2();
        if (a1 == 0) {
          if (a2 == 0) {
            // homo ref: 0b00
          } else if (a2 == 1) {
            setGenotype(&c, offset, HET);  // het: 0b01
          } else {
            setGenotype(&c, offset, MISSING);  // missing 0b10
          }
        } else if (a1 == 1) {
          if (a2 == 0) {
            setGenotype(&c, offset, HET);  // het: 0b01
          } else if (a2 == 1) {
            setGenotype(&c, offset, HOM_ALT);  // hom alt: 0b11
          } else {
            setGenotype(&c, offset, MISSING);  // missing
          }
        } else {
          // NOTE: Plink does not support tri-allelic
          // so have to set genotype as missing.
          setGenotype(&c, offset, MISSING);  // missing
        }
      }
    } else {                             // lower GD or GT
      setGenotype(&c, offset, MISSING);  // missing
    }
    if (offset == 3) {  // 3: 4 - 1, so every 4 genotype we will flush
      fwrite(&c, sizeof(char), 1, this->fpBed);
      c = 0;
    }
  }
  if (people.size() % 4 != 0)  // remaining some bits
    fwrite(&c, sizeof(char), 1, this->fpBed);

  return 0;
}

int PlinkOutputFile::writeBIM(const char* chr, const char* id, double mapDist,
                              int pos, const char* ref, const char* alt) {
  // printf("In writeBIM(), id = %s and its address is id = %p \n", id, id);

  std::string chrom = chr;
  if (atoi(chr) > 0) {
    fputs(chr, this->fpBim);
    fputc('\t', this->fpBim);
  } else if (chrom == "X")
    fputs("23\t", this->fpBim);
  else if (chrom == "Y")
    fputs("24\t", this->fpBim);
  else if (chrom == "MT")
    fputs("25\t", this->fpBim);
  else {
    fprintf(stdout, "skip chrom %s\n", chr);
    return -1;
  }
  if (id && id[0] != '.')
    fprintf(this->fpBim, "%s\t", id);
  else
    fprintf(this->fpBim, "%s:%d\t", chrom.c_str(), pos);

  fprintf(this->fpBim, "0\t");
  fprintf(this->fpBim, "%d\t", pos);
  fprintf(this->fpBim, "%s\t", ref);
  fprintf(this->fpBim, "%s\n", alt);
  return 0;
}

int PlinkOutputFile::writeBIM(const std::vector<std::string>& chr,
                              const std::vector<std::string>& id,
                              const std::vector<double>& mapDist,
                              const std::vector<int>& pos,
                              const std::vector<std::string>& ref,
                              const std::vector<std::string>& alt) {
  for (size_t i = 0; i < chr.size(); ++i) {
    this->writeBIM(chr[i].c_str(), id[i].c_str(), mapDist[i], pos[i],
                   ref[i].c_str(), alt[i].c_str());
  }
  return 0;
}

void PlinkOutputFile::writeFAM(const std::string& people) {
  fprintf(this->fpFam, "%s\t%s\t0\t0\t0\t-9\n", people.c_str(), people.c_str());
}

void PlinkOutputFile::writeFAM(const std::vector<std::string>& people) {
  for (size_t i = 0; i < people.size(); i++) {
    this->writeFAM(people[i]);
  }
}

void PlinkOutputFile::writeFAM(const std::vector<std::string>& fid,
                               const std::vector<std::string>& iid,
                               std::vector<double>& pheno) {
  assert(fid.size() == iid.size() && fid.size() == pheno.size());
  for (unsigned int i = 0; i < fid.size(); i++) {
    fprintf(this->fpFam, "%s\t%s\t0\t0\t0\t%g\n", fid[i].c_str(),
            iid[i].c_str(), pheno[i]);
  }
}

void PlinkOutputFile::writeFAM(const PlinkInputFile& pin, int i) {
  fprintf(this->fpFam, "%s\t%s\t0\t0\t%d\t%g\n", pin.getSampleName()[i].c_str(),
          pin.getSampleName()[i].c_str(), pin.getSex()[i], pin.getPheno()[i]);
}

// NOTE: m should be: marker x people
void PlinkOutputFile::writeBED(SimpleMatrix* mat, int nPeople, int nMarker) {
  /* int nPeople = mat->cols; */
  /* int nMarker = mat->rows; */
  unsigned char c = 0;
  int offset = 0;
  for (int m = 0; m < nMarker; m++) {
    for (int i = 0; i < nPeople; i++) {
      offset = i & (4 - 1);
      int geno = (int)((*mat)[m][i]);
      switch (geno) {
        case 0:
          setGenotype(&c, offset, HOM_REF);  // het: 0b01
          break;
        case 1:
          setGenotype(&c, offset, HET);  // het: 0b01
          break;
        case 2:
          setGenotype(&c, offset, HOM_ALT);  // hom alt: 0b11
          break;
        default:
          setGenotype(&c, offset, MISSING);  // missing
          break;
      }
      if (offset == 3) {  // 3: 4 - 1, so every 4 genotype we will flush
        fwrite(&c, sizeof(char), 1, this->fpBed);
        c = 0;
      }
    }                        // end for i
    if (nPeople % 4 != 0) {  // remaining some bits
      fwrite(&c, sizeof(char), 1, this->fpBed);
      c = 0;
    }
  }  // end for m
}

/**
 * This is barely minimal extraction with low efficiency.
 */
int PlinkOutputFile::extract(const std::string& prefix,
                             const std::vector<int>& sampleIdx,
                             const std::vector<int>& snpIdx) {
  PlinkInputFile pin(prefix);
  return extract(pin, sampleIdx, snpIdx);
}

int PlinkOutputFile::extract(PlinkInputFile& pin,
                             const std::vector<int>& sampleIdx,
                             const std::vector<int>& snpIdx) {
  extractFAM(pin, sampleIdx);
  extractBIM(pin, snpIdx);
  extractBED(pin, sampleIdx, snpIdx);
  return 0;
}

int PlinkOutputFile::extractFAM(PlinkInputFile& pin,
                                const std::vector<int>& sampleIdx) {
  for (size_t i = 0; i != sampleIdx.size(); ++i) {
    this->writeFAM(pin, i);  /// TODO: should also output family id
  }
  return 0;
}
int PlinkOutputFile::extractFAMWithPhenotype(PlinkInputFile& pin,
                                             const std::vector<int>& sampleIdx,
                                             const SimpleMatrix& pheno) {
  // assert((int)sampleIdx.size() == pheno.nrow());
  for (size_t i = 0; i != sampleIdx.size(); ++i) {
    /// TODO: should also output family id
    fprintf(this->fpFam, "%s\t%s\t0\t0\t%d\t%g\n",
            pin.getSampleName()[i].c_str(), pin.getSampleName()[i].c_str(),
            pin.getSex()[i], pheno[i][0]);
  }
  return 0;
}
int PlinkOutputFile::extractBIM(PlinkInputFile& pin,
                                const std::vector<int>& snpIdx) {
  for (size_t i = 0; i != snpIdx.size(); ++i) {
    this->writeBIM(pin.getChrom()[i].c_str(), pin.getMarkerName()[i].c_str(),
                   pin.getMapDist()[i], pin.getPosition()[i],
                   pin.getRef()[i].c_str(), pin.getAlt()[i].c_str());
  }
  return 0;
}
int PlinkOutputFile::extractBED(PlinkInputFile& pin,
                                const std::vector<int>& sampleIdx,
                                const std::vector<int>& snpIdx) {
  const int M = snpIdx.size();
  const int N = sampleIdx.size();
  unsigned char c = 0;
  int offset = 0;
  for (int i = 0; i < M; ++i) {  // assume SNP-major
    for (int j = 0; j < N; ++j) {
      offset = j & (4 - 1);
      setGenotype(&c, offset, pin.get2BitGenotype(sampleIdx[j], snpIdx[i]));
      if (offset == 3) {  // 3: 4 - 1, so every 4 genotype we will flush
        fwrite(&c, sizeof(char), 1, this->fpBed);
        c = 0;
      }
    }                  // end for j
    if (N % 4 != 0) {  // remaining some bits
      fwrite(&c, sizeof(char), 1, this->fpBed);
      c = 0;
    }
  }  // end for i
  return 0;
}
