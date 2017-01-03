#include "libVcf/PlinkOutputFile.h"

#include "base/SimpleMatrix.h"
#include "libVcf/PlinkInputFile.h"

int PlinkOutputFile::writeRecord(VCFRecord* r) {
  int ret;
  // write BIM
  // printf("id= %s and its address id = %p\n", r->getID(), r->getID());
  ret = this->writeBIM(r->getChrom(), r->getID(), 0, r->getPos(), r->getRef(),
                       r->getAlt());
  if (ret) return ret;  // unsuccess
  // printf("id= %s and its address id = %p\n", r->getID(), r->getID());

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
  if (people.size() % 4 != 0)  // remaining some bits
    fwrite(&c, sizeof(char), 1, this->fpBed);

  return 0;
}

int PlinkOutputFile::writeRecordWithFilter(VCFRecord* r, const double minGD,
                                           const double minGQ) {
  int ret;
  // write BIM
  // printf("id= %s and its address id = %p\n", r->getID(), r->getID());
  ret = this->writeBIM(r->getChrom(), r->getID(), 0, r->getPos(), r->getRef(),
                       r->getAlt());
  if (ret) return ret;  // unsuccess
  // printf("id= %s and its address id = %p\n", r->getID(), r->getID());

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
  int refLen = strlen(ref);
  int altLen = strlen(alt);
  if (refLen > 1) {
    if (ref[1] != ',') {
      fprintf(stdout, "skip with ref = %s and alt = %s\n", ref, alt);
      return -1;
    }
  }
  if (altLen > 1) {
    if (alt[1] != ',') {
      fprintf(stdout, "skip with ref = %s and alt= %s\n", ref, alt);
      return -1;
    }
  }

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
                             const std::vector<int>& snpIdx,
                             const std::vector<int>& sampleIdx) {
  PlinkInputFile pin(prefix);

  for (size_t i = 0; i != sampleIdx.size(); ++i) {
    this->writeFAM(pin.getIID()[i]);
  }

  for (size_t i = 0; i != snpIdx.size(); ++i) {
    this->writeBIM(pin.getChrom()[i].c_str(), pin.getMarkerName()[i].c_str(),
                   pin.getMapDist()[i], pin.getPosition()[i],
                   pin.getRef()[i].c_str(), pin.getAlt()[i].c_str());
  }

  const int M = snpIdx.size();
  const int N = sampleIdx.size();
  unsigned char c = 0;
  int offset = 0;
  for (int i = 0; i < M; ++i) {
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
