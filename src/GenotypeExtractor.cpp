#include "GenotypeExtractor.h"

#include "Result.h"
#include "base/Logger.h"

#include "libVcf/VCFUtil.h"
#include "libsrc/MathVector.h"
#include "libsrc/MathMatrix.h"

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
      sex(NULL),
      claytonCoding(true) {
  this->vin = new VCFExtractor(fn.c_str());
}

GenotypeExtractor::~GenotypeExtractor() {
  if (this->vin) {
    delete this->vin;
    this->vin = NULL;
  }
}

int GenotypeExtractor::extractMultipleGenotype(Matrix* g) {
  static Matrix m;  // make it static to reduce memory allocation
  int row = 0;
  std::vector<std::string> colNames;
  std::string name;
  this->hemiRegion.clear();
  while (this->vin->readRecord()) {
    VCFRecord& r = this->vin->getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    m.Dimension(row + 1, people.size());

    int genoIdx;
    const bool useDosage = (!this->dosageTag.empty());
    if (useDosage) {
      genoIdx = r.getFormatIndex(dosageTag.c_str());
    } else {
      genoIdx = r.getFormatIndex("GT");
    }
    int GDidx = r.getFormatIndex("GD");
    int GQidx = r.getFormatIndex("GQ");
    assert(this->parRegion);
    bool hemiRegion = this->parRegion->isHemiRegion(r.getChrom(), r.getPos());
    // e.g.: Loop each (selected) people in the same order as in the VCF
    const int numPeople = (int)people.size();
    for (int i = 0; i < numPeople; i++) {
      indv = people[i];
      // get GT index. if you are sure the index will not change, call this
      // function only once!
      if (genoIdx >= 0) {
        // printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first
        // field of each individual
        if (useDosage) {
          if (!hemiRegion) {
            m[row][i] = indv->justGet(genoIdx).toDouble();
          } else {
            // for male hemi region, imputated dosage is usually between 0 and 1
            // need to multiply by 2.0
            if ((*sex)[i] == PLINK_MALE) {
              m[row][i] = indv->justGet(genoIdx).toDouble() * 2.0;
            }
          }
        } else {
          if (!hemiRegion) {
            m[row][i] = indv->justGet(genoIdx).getGenotype();
          } else {
            if ((*sex)[i] == PLINK_MALE) {
              m[row][i] = indv->justGet(genoIdx).getMaleNonParGenotype02();
            } else if ((*sex)[i] == PLINK_FEMALE) {
              m[row][i] = indv->justGet(genoIdx).getGenotype();
            } else {
              m[row][i] = MISSING_GENOTYPE;
            }
          }
        }
        if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
          m[row][i] = MISSING_GENOTYPE;
          continue;
        }
      } else {
        logger->error("Cannot find %s field!",
                      this->dosageTag.empty() ? "GT" : dosageTag.c_str());
        return -1;
      }
    }

    // check frequency cutoffs
    int numNonMissingPeople = 0;
    double maf = 0.;
    for (int i = 0; i < numPeople; ++i) {
      if (m[row][i] < 0) continue;
      maf += m[row][i];
      ++numNonMissingPeople;
    }
    if (numNonMissingPeople) {
      maf = maf / (2. * numNonMissingPeople);
    } else {
      maf = 0.0;
    }
    if (maf > .5) {
      maf = 1.0 - maf;
    }
    if (this->freqMin > 0. && this->freqMin > maf) continue;
    if (this->freqMax > 0. && this->freqMax < maf) continue;

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
  }

  // delete rows (ugly code here, as we may allocate extra row in previous
  // loop)
  m.Dimension(row, m.cols);

  // now transpose (marker by people -> people by marker)
  g->Transpose(m);
  for (int i = 0; i < row; ++i) {
    g->SetColumnLabel(i, colNames[i].c_str());
  }
  return SUCCEED;
} // end GenotypeExtractor

int GenotypeExtractor::extractSingleGenotype(Matrix* g, Result* b) {
  Matrix& genotype = *g;
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

  genotype.Dimension(people.size(), 1);

  // get GT index. if you are sure the index will not change, call this
  // function only once!
  const bool useDosage = (!this->dosageTag.empty());
  int genoIdx;
  if (useDosage) {
    genoIdx = r.getFormatIndex(dosageTag.c_str());
  } else {
    genoIdx = r.getFormatIndex("GT");
  }
  // int GTidx = r.getFormatIndex("GT");
  int GDidx = r.getFormatIndex("GD");
  int GQidx = r.getFormatIndex("GQ");

  bool hemiRegion = this->parRegion->isHemiRegion(r.getChrom(), r.getPos());
  // e.g.: Loop each (selected) people in the same order as in the VCF
  const int numPeople = (int)people.size();
  for (int i = 0; i < numPeople; i++) {
    indv = people[i];

    if (genoIdx >= 0) {
      // printf("%s ", indv->justGet(0).toStr());  // [0] meaning the first
      // field of each individual
      if (useDosage) {
        genotype[i][0] = indv->justGet(genoIdx).toDouble();
      } else {
        if (!hemiRegion) {
          genotype[i][0] = indv->justGet(genoIdx).getGenotype();
        } else {
          if ((*sex)[i] == PLINK_MALE) {
            genotype[i][0] = indv->justGet(genoIdx).getMaleNonParGenotype02();
          } else if ((*sex)[i] == PLINK_FEMALE) {
            genotype[i][0] = indv->justGet(genoIdx).getGenotype();
          } else {
            genotype[i][0] = MISSING_GENOTYPE;
          }
        }
      }
      if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
        genotype[i][0] = MISSING_GENOTYPE;
        continue;
      }
      // logger->info("%d ", int(genotype[i][0]));
    } else {
      logger->error(
          "Cannot find [ %s ] field when read individual information [ %s ]!",
          this->dosageTag.empty() ? "GT" : this->dosageTag.c_str(),
          indv->getSelf().toStr());
      return ERROR;
    }
  }

  // check frequency cutoffs
  double maf = 0.;
  for (int i = 0; i < numPeople; ++i) {
    maf += genotype[i][0];
  }
  maf = maf / (2. * numPeople);
  if (maf > .5) {
    maf = 1.0 - maf;
  }
  if (this->freqMin > 0. && this->freqMin > maf) return FAIL_FILTER;
  if (this->freqMax > 0. && this->freqMax < maf) return FAIL_FILTER;

  std::string label = r.getChrom();
  label += ':';
  label += r.getPosStr();
  genotype.SetColumnLabel(0, label.c_str());

  this->hemiRegion.resize(1);
  assert(this->parRegion);
  if (this->parRegion &&
      this->parRegion->isHemiRegion(r.getChrom(), r.getPos())) {
    this->hemiRegion[0] = true;
  } else {
    this->hemiRegion[0] = false;
  }
  return SUCCEED;
} // end extractSingleGenotype()

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
bool GenotypeExtractor::checkGD(VCFIndividual* indv, int gdIdx) {
  if (!needGD) return true;
  int gd = indv->justGet(gdIdx).toInt();
  if (this->GDmin > 0 && gd < this->GDmin) return false;
  if (this->GDmax > 0 && gd > this->GDmax) return false;
  return true;
}
bool GenotypeExtractor::checkGQ(VCFIndividual* indv, int gqIdx) {
  if (!needGQ) return true;
  int gq = indv->justGet(gqIdx).toInt();
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
