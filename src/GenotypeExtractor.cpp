#include "GenotypeExtractor.h"

#include "GenotypeCounter.h"
#include "Result.h"
#include "base/Logger.h"

#include "libVcf/VCFUtil.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/MathVector.h"

extern Logger* logger;

GenotypeExtractor::GenotypeExtractor(const std::string& fn)
    : vin(NULL),
      freqMin(-1),
      freqMax(-1),
      GDmin(-1),
      GDmax(-1),
      needGD(false),
      GQmin(-1),
      GQmax(-1),
      needGQ(false),
      parRegion(NULL),
      sex(NULL) {
  this->vin = new VCFExtractor(fn.c_str());
}

GenotypeExtractor::~GenotypeExtractor() {
  if (this->vin) {
    delete this->vin;
    this->vin = NULL;
  }
}

int GenotypeExtractor::extractMultipleGenotype(Matrix* g) {
  // static Matrix m;  // make it static to reduce memory allocation
  int row = 0;
  // std::vector<std::string> colNames;
  // std::string name;
  //   GenotypeCounter genoCounter;
  this->genotype.clear();
  this->altAlleleToParse = 0;
  while (true) {
    if (altAlleleToParse == 0) {
      if (!this->vin->readRecord()) {
        // finished reading
        break;
      }
    }
    // parse alt allele one at a time
    VCFRecord& r = this->vin->getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    countAltAllele(r.getAlt());

    this->sampleSize = people.size();
    // m.Dimension(row + 1, people.size());
    row++;
    this->variantName.resize(row);
    this->counter.resize(row);
    this->hemiRegion.resize(row);

    // get GT index and cannot assume it is a constant
    int genoIdx;
    const bool useDosage = (!this->dosageTag.empty());
    if (useDosage) {
      genoIdx = r.getFormatIndex(dosageTag.c_str());
    } else {
      genoIdx = r.getFormatIndex("GT");
    }
    if (useDosage && altAlleleToParse > 1) {
      logger->error(
          "Unsupported scenario: multiple alleles and use dosage as "
          "genotypes! - we will only use dosage");
      altAlleleToParse = 1;
    }
    int GDidx = r.getFormatIndex("GD");
    int GQidx = r.getFormatIndex("GQ");
    assert(this->parRegion);
    bool isHemiRegion = this->parRegion->isHemiRegion(r.getChrom(), r.getPos());
    // e.g.: Loop each (selected) people in the same order as in the VCF
    for (int i = 0; i < sampleSize; i++) {
      indv = people[i];
      const double geno = getGenotype(*indv, useDosage, isHemiRegion, (*sex)[i],
                                      genoIdx, GDidx, GQidx);
      genotype.push_back(geno);
      counter.back().add(geno);

#if 0
      if (genoIdx >= 0) {
        if (useDosage) {
          if (!isHemiRegion) {
            // m[row][i] = indv->justGet(genoIdx).toDouble();
            genotype.push_back(indv->justGet(genoIdx).toDouble());
          } else {
            // for male hemi region, imputated dosage is usually between 0 and 1
            // need to multiply by 2.0
            if ((*sex)[i] == PLINK_MALE) {
              //m[row][i] = indv->justGet(genoIdx).toDouble() * 2.0;
              genotype.push_back(indv->justGet(genoIdx).toDouble() * 2.0);
            }
          }
        } else { // use hard-coded genotypes
          if (!isHemiRegion) {
            //m[row][i] = indv->justGet(genoIdx).getGenotype();
            genotype.push_back(indv->justGet(genoIdx).getGenotype());
          } else {
            if ((*sex)[i] == PLINK_MALE) {
              // m[row][i] = indv->justGet(genoIdx).getMaleNonParGenotype02();
              genotype.push_back(indv->justGet(genoIdx).getMaleNonParGenotype02());
            } else if ((*sex)[i] == PLINK_FEMALE) {
              // m[row][i] = indv->justGet(genoIdx).getGenotype();
              genotype.push_back(indv->justGet(genoIdx).getGenotype());
            } else {
              // m[row][i] = MISSING_GENOTYPE;
              genotype.push_back(MISSING_GENOTYPE);
            }
          }
        }
        if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
          // m[row][i] = MISSING_GENOTYPE;
          genotype.back() = MISSING_GENOTYPE;
        }
        //genoCounter.add(m[row][i]);
        genoCounter.add(genotype.back());
      } else {
        logger->error("Cannot find %s field!",
                      this->dosageTag.empty() ? "GT" : dosageTag.c_str());
        return -1;
      }
#endif
    }  // end for i

    // check frequency cutoffs
    const double maf = counter.back().getMAF();
    if (this->freqMin > 0. && this->freqMin > maf) continue;
    if (this->freqMax > 0. && this->freqMax < maf) continue;

    this->variantName.back() = r.getChrom();
    this->variantName.back() += ":";
    this->variantName.back() += r.getPosStr();
    this->hemiRegion.back() = (isHemiRegion);

#if 0
    // store genotype results
    name = r.getChrom();
    name += ":";
    name += r.getPosStr();
    colNames.push_back(name);
    ++row;

    assert(this->parRegion);
    if (this->parRegion &&
        this->parRegion->isHemiRegion(r.getChrom(), r.getPos())) {
      this->hemiRegion.push_back(true);
    } else {
      this->hemiRegion.push_back(false);
    }
    this->counter.push_back(genoCounter);
#endif
    this->altAlleleToParse--;
  }  // end while (this->vin->readRecord())

  // // delete rows (ugly code here, as we may allocate extra row in previous
  // // loop)
  // m.Dimension(row, m.cols);
  // todo: store to g

  // now transpose (marker by people -> people by marker)
  // g->Transpose(m);
  assert((int)genotype.size() == this->sampleSize * row);
  assign(this->genotype, sampleSize, row, g);
  for (int i = 0; i < row; ++i) {
    g->SetColumnLabel(i, variantName[i].c_str());
  }
  return SUCCEED;
}  // end GenotypeExtractor

int GenotypeExtractor::extractSingleGenotype(Matrix* g, Result* b) {
  // Matrix& genotype = *g;
  this->genotype.clear();
  Result& buf = *b;

  bool hasRead = this->vin->readRecord();
  if (!hasRead) return FILE_END;

  VCFRecord& r = this->vin->getVCFRecord();
  VCFPeople& people = r.getPeople();
  VCFIndividual* indv;

  buf.updateValue("CHROM", r.getChrom());
  buf.updateValue("POS", r.getPosStr());
  buf.updateValue("REF", r.getRef());
  buf.updateValue("ALT", r.getAlt());

  // genotype.Dimension(people.size(), 1);
  this->sampleSize = people.size();
  this->variantName.resize(1);
  this->counter.resize(1);
  this->hemiRegion.resize(1);

  // get GT index. if you are sure the index will not change, call this
  // function only once!
  const bool useDosage = (!this->dosageTag.empty());
  int genoIdx;
  if (useDosage) {
    genoIdx = r.getFormatIndex(dosageTag.c_str());
  } else {
    genoIdx = r.getFormatIndex("GT");
  }
  int GDidx = r.getFormatIndex("GD");
  int GQidx = r.getFormatIndex("GQ");

  bool isHemiRegion = this->parRegion->isHemiRegion(r.getChrom(), r.getPos());
  // e.g.: Loop each (selected) people in the same order as in the VCF
  for (int i = 0; i < sampleSize; i++) {
    indv = people[i];
    const double geno = getGenotype(*indv, useDosage, isHemiRegion, (*sex)[i],
                                    genoIdx, GDidx, GQidx);
    genotype.push_back(geno);
    counter.back().add(geno);

#if 0
    if (genoIdx >= 0) {
      if (useDosage) {
        if (!hemiRegion) {
          genotype[i][0] = indv->justGet(genoIdx).toDouble();
        } else {
          // for male hemi region, imputated dosage is usually between 0 and 1
          // need to multiply by 2.0
          if ((*sex)[i] == PLINK_MALE) {
            genotype[i][0] = indv->justGet(genoIdx).toDouble() * 2.0;
          }
        }
      } else {
        if (!hemiRegion) {
          genotype[i][0] = indv->justGet(genoIdx).getGenotype();
        } else {  // use hard-coded genotypes
          if ((*sex)[i] == PLINK_MALE) {
            genotype[i][0] = indv->justGet(genoIdx).getMaleNonParGenotype02();
          } else if ((*sex)[i] == PLINK_FEMALE) {
            genotype[i][0] = indv->justGet(genoIdx).getGenotype();
          } else {
            genotype[i][0] = MISSING_GENOTYPE;
          }
        }
      }
      if (!checkGD(*indv, GDidx) || !checkGQ(*indv, GQidx)) {
        genotype[i][0] = MISSING_GENOTYPE;
      }
      counter[0].add(genotype[i][0]);
      // logger->info("%d ", int(genotype[i][0]));
    } else {
      std::string s;
      indv->toStr(&s);
      logger->error(
          "Cannot find [ %s ] field when read individual information [ %s ]!",
          this->dosageTag.empty() ? "GT" : this->dosageTag.c_str(), s.c_str());
      return ERROR;
    }
#endif
  }

  // check frequency cutoffs
  const double maf = counter[0].getMAF();
  if (this->freqMin > 0. && this->freqMin > maf) return FAIL_FILTER;
  if (this->freqMax > 0. && this->freqMax < maf) return FAIL_FILTER;

  variantName.back() = r.getChrom();
  variantName.back() += ':';
  variantName.back() += r.getPosStr();
  hemiRegion.back() = isHemiRegion;

  assert((int)genotype.size() == sampleSize);
  assign(genotype, sampleSize, 1, g);
  g->SetColumnLabel(0, variantName.back().c_str());

  // this->hemiRegion.resize(1);
  // assert(this->parRegion);
  // if (this->parRegion &&
  //     this->parRegion->isHemiRegion(r.getChrom(), r.getPos())) {
  //   this->hemiRegion[0] = true;
  // } else {
  //   this->hemiRegion[0] = false;
  // }
  return SUCCEED;
}  // end extractSingleGenotype()

int loadMarkerFromVCF(const std::string& fileName, const std::string& marker,
                      std::vector<std::string>* rowLabel, Matrix* genotype) {
  if (!rowLabel || !genotype) {
    // invalid parameter
    return -1;
  }
  Matrix& m = *genotype;
  int col = 0;

  VCFInputFile vin(fileName);
  vin.setRangeList(marker);

  while (vin.readRecord()) {
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    m.Dimension(people.size(), col + 1);

    int GTidx = r.getFormatIndex("GT");
    for (int i = 0; i < (int)people.size(); i++) {
      indv = people[i];
      // get GT index. if you are sure the index will not change,
      // call this function only once!
      if (GTidx >= 0) {
        // printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first
        // field of each individual
        m[i][col] = indv->justGet(GTidx).getGenotype();
      } else {
        logger->error("Cannot find GT field!");
        return -1;
      }
    }
    if (col == 0) {
      // set-up names
      rowLabel->resize(people.size());
      for (size_t i = 0; i < people.size(); ++i) {
        (*rowLabel)[i] = people[i]->getName();
      }
    }
    std::string colLabel = r.getChrom();
    colLabel += ":";
    colLabel += r.getPosStr();
    m.SetColumnLabel(col, colLabel.c_str());
    ++col;
  }

  return 0;
}

bool GenotypeExtractor::setSiteFreqMin(const double f) {
  if (f < 0.0 || f > 1.0) {
    return false;
  }
  this->freqMin = f - 1e-10;  // allow rounding error
  return true;
}
bool GenotypeExtractor::setSiteFreqMax(const double f) {
  if (f < 0.0 || f > 1.0) {
    return false;
  }
  this->freqMax = f + 1e-10;  // allow rounding error
  return true;
}

void GenotypeExtractor::setSiteDepthMin(int d) {
  this->vin->setSiteDepthMin(d);
}
void GenotypeExtractor::setSiteDepthMax(int d) {
  this->vin->setSiteDepthMax(d);
}

// @return true if GD is valid
// if GD is missing, we will take GD = 0
bool GenotypeExtractor::checkGD(VCFIndividual& indv, int gdIdx) {
  if (!needGD) return true;
  int gd = indv.justGet(gdIdx).toInt();
  if (this->GDmin > 0 && gd < this->GDmin) return false;
  if (this->GDmax > 0 && gd > this->GDmax) return false;
  return true;
}
bool GenotypeExtractor::checkGQ(VCFIndividual& indv, int gqIdx) {
  if (!needGQ) return true;
  int gq = indv.justGet(gqIdx).toInt();
  if (this->GQmin > 0 && gq < this->GQmin) return false;
  if (this->GQmax > 0 && gq > this->GQmax) return false;
  return true;
}
void GenotypeExtractor::setGDmin(int m) {
  this->needGD = true;
  this->GDmin = m;
}
void GenotypeExtractor::setGDmax(int m) {
  this->needGD = true;
  this->GDmax = m;
}
void GenotypeExtractor::setGQmin(int m) {
  this->needGQ = true;
  this->GQmin = m;
}
void GenotypeExtractor::setGQmax(int m) {
  this->needGQ = true;
  this->GQmax = m;
}

void GenotypeExtractor::setSiteQualMin(int q) { this->vin->setSiteQualMin(q); }
void GenotypeExtractor::setSiteMACMin(int n) { this->vin->setSiteMACMin(n); }
int GenotypeExtractor::setAnnoType(const std::string& s) {
  return this->vin->setAnnoType(s.c_str());
}

void GenotypeExtractor::setRange(const RangeList& l) { this->vin->setRange(l); }
void GenotypeExtractor::setRangeList(const std::string& l) {
  this->vin->setRangeList(l);
}
void GenotypeExtractor::setRangeFile(const std::string& fn) {
  this->vin->setRangeFile(fn.c_str());
}
void GenotypeExtractor::includePeople(const std::string& v) {
  this->vin->includePeople(v.c_str());
}
void GenotypeExtractor::includePeople(const std::vector<std::string>& v) {
  this->vin->includePeople(v);
}
void GenotypeExtractor::includePeopleFromFile(const std::string& fn) {
  this->vin->includePeopleFromFile(fn.c_str());
}
void GenotypeExtractor::excludePeople(const std::string& v) {
  this->vin->excludePeople(v.c_str());
}
void GenotypeExtractor::excludePeople(const std::vector<std::string>& v) {
  this->vin->excludePeople(v);
}
void GenotypeExtractor::excludePeopleFromFile(const std::string& fn) {
  this->vin->excludePeopleFromFile(fn.c_str());
}
void GenotypeExtractor::excludePeople(const std::vector<std::string>& sample,
                                      const std::vector<int>& index) {
  for (size_t i = 0; i < index.size(); ++i) {
    this->excludePeople(sample[i]);
  }
}

void GenotypeExtractor::excludeAllPeople() { this->vin->excludeAllPeople(); }
void GenotypeExtractor::enableAutoMerge() { this->vin->enableAutoMerge(); }
void GenotypeExtractor::getPeopleName(std::vector<std::string>* p) {
  return this->vin->getVCFHeader()->getPeopleName(p);
}

void GenotypeExtractor::getIncludedPeopleName(
    std::vector<std::string>* p) const {
  this->vin->getIncludedPeopleName(p);
  return;
}

void GenotypeExtractor::countAltAllele(const char* s) {
  stringTokenize(s, ",", &altAllele);
  altAlleleToParse = altAllele.size();
}

double GenotypeExtractor::getGenotype(VCFIndividual& indv, const bool useDosage,
                                      const bool hemiRegion, const int sex,
                                      const int genoIdx, const int GDidx,
                                      const int GQidx) {
  double ret;
  if (genoIdx >= 0) {
    if (useDosage) {
      if (!hemiRegion) {
        // m[row][i] = indv->justGet(genoIdx).toDouble();
        ret = (indv.justGet(genoIdx).toDouble());
      } else {
        // for male hemi region, imputated dosage is usually between 0 and 1
        // need to multiply by 2.0
        if (sex == PLINK_MALE) {
          // m[row][i] = indv.justGet(genoIdx).toDouble() * 2.0;
          ret = (indv.justGet(genoIdx).toDouble() * 2.0);
        } else {
          ret = (indv.justGet(genoIdx).toDouble());
        }
      }
    } else {  // use hard-coded genotypes
      if (!hemiRegion) {
        // m[row][i] = indv.justGet(genoIdx).getGenotype();
        ret = (indv.justGet(genoIdx).getGenotype());
      } else {
        if (sex == PLINK_MALE) {
          // m[row][i] = indv.justGet(genoIdx).getMaleNonParGenotype02();
          ret = (indv.justGet(genoIdx).getMaleNonParGenotype02());
        } else if (sex == PLINK_FEMALE) {
          // m[row][i] = indv.justGet(genoIdx).getGenotype();
          ret = (indv.justGet(genoIdx).getGenotype());
        } else {
          // m[row][i] = MISSING_GENOTYPE;
          ret = (MISSING_GENOTYPE);
        }
      }
    }
    if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
      // m[row][i] = MISSING_GENOTYPE;
      return MISSING_GENOTYPE;
    }
    // genoCounter.add(m[row][i]);
    // genoCounter.add(genotype.back());
  } else {
    logger->error("Cannot find %s field!",
                  this->dosageTag.empty() ? "GT" : dosageTag.c_str());
    // return -1;
    return MISSING_GENOTYPE;
  }
}

void GenotypeExtractor::assign(const std::vector<double>& from, int nrow,
                               int ncol, Matrix* to) {
  assert(to);
  Matrix& out = *to;
  out.Dimension(nrow, ncol);
  assert((int)from.size() == nrow * ncol);
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      out[i][j] = from[nrow * j + i];
    }
  }
}
