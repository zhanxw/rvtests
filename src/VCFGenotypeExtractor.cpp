#include "VCFGenotypeExtractor.h"

#include "GenotypeCounter.h"
#include "Result.h"

#include "base/Argument.h"
#include "base/Logger.h"
#include "libVcf/VCFUtil.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/MathVector.h"

DECLARE_BOOL_PARAMETER(outputID);

extern Logger* logger;

VCFGenotypeExtractor::VCFGenotypeExtractor(const std::string& fn)
    : GenotypeExtractor(fn), vin(NULL), altAlleleToParse(-1) {
  this->vin = new VCFExtractor(fn.c_str());
}

VCFGenotypeExtractor::~VCFGenotypeExtractor() {
  if (this->vin) {
    delete this->vin;
    this->vin = NULL;
  }
}

int VCFGenotypeExtractor::extractMultipleGenotype(Matrix* g) {
  assert(g);
  assert(g->rows >= 0 && g->cols >= 0);
  g->Dimension(0, 0);
  int row = 0;
  this->genotype.clear();
  this->altAlleleToParse = -1;

  while (true) {
    if (this->altAlleleToParse <= 0) {
      if (!this->vin->readRecord()) {
        // reached the end
        break;
      } else {
        const char* alt = this->vin->getVCFRecord().getAlt();
        if (multiAllelicMode) {
          parseAltAllele(alt);
          this->altAlleleToParse = altAllele.size();
        } else {
          this->altAllele.resize(1);
          this->altAllele[0] = alt;
          this->altAlleleToParse = 1;
        }
      }
    }
    assert(this->altAlleleToParse > 0);
    // parse alt allele one at a time
    VCFRecord& r = this->vin->getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    this->sampleSize = people.size();
    row++;
    this->variantName.resize(row);
    this->counter.resize(row);
    this->counter.back().reset();
    this->hemiRegion.resize(row);

    // get GT index and cannot assume it is a constant across variants
    int genoIdx;
    const bool useDosage = (!this->dosageTag.empty());
    if (useDosage) {
      genoIdx = r.getFormatIndex(dosageTag.c_str());
    } else {
      genoIdx = r.getFormatIndex("GT");
    }
    if (useDosage && multiAllelicMode) {
      logger->error(
          "Unsupported scenario: multiple mode and use dosage as "
          "genotypes! - we will only use dosage");
      this->altAlleleToParse = 1;
    }
    const int GDidx = r.getFormatIndex("GD");
    const int GQidx = r.getFormatIndex("GQ");
    assert(this->parRegion);
    const bool isHemiRegion =
        this->parRegion->isHemiRegion(r.getChrom(), r.getPos());
    // e.g.: Loop each (selected) people in the same order as in the VCF
    double geno;
    const int altAlleleGT = this->altAllele.size() - this->altAlleleToParse + 1;
    for (int i = 0; i < sampleSize; i++) {
      indv = people[i];
      if (multiAllelicMode) {
        geno =
            getGenotypeForAltAllele(*indv, useDosage, isHemiRegion, (*sex)[i],
                                    genoIdx, GDidx, GQidx, altAlleleGT);
      } else {
        geno = getGenotype(*indv, useDosage, isHemiRegion, (*sex)[i], genoIdx,
                           GDidx, GQidx);
      }
      this->genotype.push_back(geno);
      this->counter.back().add(geno);
    }  // end for i

    // check frequency cutoffs
    const double maf = counter.back().getMAF();
    if ((this->freqMin > 0. && this->freqMin > maf) ||
        (this->freqMax > 0. && this->freqMax < maf)) {
      // undo loaded contents
      row--;
      this->variantName.resize(row);
      this->counter.resize(row);
      this->hemiRegion.resize(row);
      this->genotype.resize(this->genotype.size() - this->sampleSize);
      this->altAlleleToParse--;
      continue;
    }

    this->variantName.back() = r.getChrom();
    this->variantName.back() += ":";
    this->variantName.back() += r.getPosStr();
    if (multiAllelicMode) {
      this->variantName.back() += r.getRef();
      this->variantName.back() += "/";
      this->variantName.back() +=
          altAllele[altAllele.size() - this->altAlleleToParse];
    }
    this->hemiRegion.back() = (isHemiRegion);

    this->altAlleleToParse--;
  }  // end while (this->vin->readRecord())

  // now transpose (marker by people -> people by marker)
  if (row > 0) {
    assert((int)genotype.size() == this->sampleSize * row);
    assign(this->genotype, sampleSize, row, g);
    for (int i = 0; i < row; ++i) {
      g->SetColumnLabel(i, variantName[i].c_str());
    }
  }
  return SUCCEED;
}  // end VCFGenotypeExtractor

int VCFGenotypeExtractor::extractSingleGenotype(Matrix* g, Result* b) {
  this->genotype.clear();
  Result& buf = *b;

  if (this->altAlleleToParse <= 0) {
    if (!this->vin->readRecord()) {
      return FILE_END;
    } else {
      const char* alt = this->vin->getVCFRecord().getAlt();
      if (multiAllelicMode) {
        parseAltAllele(alt);
        this->altAlleleToParse = altAllele.size();
      } else {
        this->altAllele.resize(1);
        this->altAllele[0] = alt;
        this->altAlleleToParse = 1;
      }
    }
  }
  assert(this->altAlleleToParse >= 0);
  VCFRecord& r = this->vin->getVCFRecord();
  VCFPeople& people = r.getPeople();
  VCFIndividual* indv;

  buf.updateValue("CHROM", r.getChrom());
  buf.updateValue("POS", r.getPosStr());
  if (FLAG_outputID) {
    buf.updateValue("ID", r.getID());
  }
  buf.updateValue("REF", r.getRef());
  buf.updateValue("ALT", altAllele[altAllele.size() - this->altAlleleToParse]);

  // genotype.Dimension(people.size(), 1);
  this->sampleSize = people.size();
  this->variantName.resize(1);
  this->counter.resize(1);
  this->counter.back().reset();
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
  const int GDidx = r.getFormatIndex("GD");
  const int GQidx = r.getFormatIndex("GQ");
  bool isHemiRegion = this->parRegion->isHemiRegion(r.getChrom(), r.getPos());
  // e.g.: Loop each (selected) people in the same order as in the VCF
  double geno;
  const int altAlleleGT = this->altAllele.size() - this->altAlleleToParse + 1;
  for (int i = 0; i < sampleSize; i++) {
    indv = people[i];
    if (multiAllelicMode) {
      geno = getGenotypeForAltAllele(*indv, useDosage, isHemiRegion, (*sex)[i],
                                     genoIdx, GDidx, GQidx, altAlleleGT);
    } else {
      geno = getGenotype(*indv, useDosage, isHemiRegion, (*sex)[i], genoIdx,
                         GDidx, GQidx);
    }
    genotype.push_back(geno);
    counter.back().add(geno);
  }

  // check frequency cutoffs
  const double maf = counter[0].getMAF();
  if ((this->freqMin > 0. && this->freqMin > maf) ||
      (this->freqMax > 0. && this->freqMax < maf)) {
    --this->altAlleleToParse;
    return FAIL_FILTER;
  }

  variantName.back() = r.getChrom();
  variantName.back() += ':';
  variantName.back() += r.getPosStr();
  hemiRegion.back() = isHemiRegion;

  assert((int)genotype.size() == sampleSize);
  assign(genotype, sampleSize, 1, g);
  g->SetColumnLabel(0, variantName.back().c_str());

  --this->altAlleleToParse;
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

bool VCFGenotypeExtractor::setSiteFreqMin(const double f) {
  if (f < 0.0 || f > 1.0) {
    return false;
  }
  this->freqMin = f - 1e-10;  // allow rounding error
  return true;
}
bool VCFGenotypeExtractor::setSiteFreqMax(const double f) {
  if (f < 0.0 || f > 1.0) {
    return false;
  }
  this->freqMax = f + 1e-10;  // allow rounding error
  return true;
}

void VCFGenotypeExtractor::setSiteDepthMin(int d) {
  this->vin->setSiteDepthMin(d);
}
void VCFGenotypeExtractor::setSiteDepthMax(int d) {
  this->vin->setSiteDepthMax(d);
}

// @return true if GD is valid
// if GD is missing, we will take GD = 0
bool VCFGenotypeExtractor::checkGD(VCFIndividual& indv, int gdIdx) {
  if (!needGD) return true;
  int gd = indv.justGet(gdIdx).toInt();
  if (this->GDmin > 0 && gd < this->GDmin) return false;
  if (this->GDmax > 0 && gd > this->GDmax) return false;
  return true;
}
bool VCFGenotypeExtractor::checkGQ(VCFIndividual& indv, int gqIdx) {
  if (!needGQ) return true;
  int gq = indv.justGet(gqIdx).toInt();
  if (this->GQmin > 0 && gq < this->GQmin) return false;
  if (this->GQmax > 0 && gq > this->GQmax) return false;
  return true;
}
void VCFGenotypeExtractor::setGDmin(int m) {
  this->needGD = true;
  this->GDmin = m;
}
void VCFGenotypeExtractor::setGDmax(int m) {
  this->needGD = true;
  this->GDmax = m;
}
void VCFGenotypeExtractor::setGQmin(int m) {
  this->needGQ = true;
  this->GQmin = m;
}
void VCFGenotypeExtractor::setGQmax(int m) {
  this->needGQ = true;
  this->GQmax = m;
}

void VCFGenotypeExtractor::setSiteFile(const std::string& fn) {
  this->vin->setSiteFile(fn);
}

void VCFGenotypeExtractor::setSiteQualMin(int q) {
  this->vin->setSiteQualMin(q);
}
void VCFGenotypeExtractor::setSiteMACMin(int n) { this->vin->setSiteMACMin(n); }
int VCFGenotypeExtractor::setAnnoType(const std::string& s) {
  return this->vin->setAnnoType(s.c_str());
}

void VCFGenotypeExtractor::setRange(const RangeList& l) {
  this->vin->setRange(l);
}
void VCFGenotypeExtractor::setRangeList(const std::string& l) {
  this->vin->setRangeList(l);
}
void VCFGenotypeExtractor::setRangeFile(const std::string& fn) {
  this->vin->setRangeFile(fn.c_str());
}
void VCFGenotypeExtractor::includePeople(const std::string& v) {
  this->vin->includePeople(v.c_str());
}
void VCFGenotypeExtractor::includePeople(const std::vector<std::string>& v) {
  this->vin->includePeople(v);
}
void VCFGenotypeExtractor::includePeopleFromFile(const std::string& fn) {
  this->vin->includePeopleFromFile(fn.c_str());
}
void VCFGenotypeExtractor::excludePeople(const std::string& v) {
  this->vin->excludePeople(v.c_str());
}
void VCFGenotypeExtractor::excludePeople(const std::vector<std::string>& v) {
  this->vin->excludePeople(v);
}
void VCFGenotypeExtractor::excludePeopleFromFile(const std::string& fn) {
  this->vin->excludePeopleFromFile(fn.c_str());
}
void VCFGenotypeExtractor::excludePeople(const std::vector<std::string>& sample,
                                         const std::vector<int>& index) {
  for (size_t i = 0; i < index.size(); ++i) {
    this->excludePeople(sample[i]);
  }
}

void VCFGenotypeExtractor::excludeAllPeople() { this->vin->excludeAllPeople(); }
void VCFGenotypeExtractor::enableAutoMerge() { this->vin->enableAutoMerge(); }
void VCFGenotypeExtractor::getPeopleName(std::vector<std::string>* p) {
  return this->vin->getVCFHeader()->getPeopleName(p);
}

void VCFGenotypeExtractor::getIncludedPeopleName(
    std::vector<std::string>* p) const {
  this->vin->getIncludedPeopleName(p);
  return;
}

void VCFGenotypeExtractor::parseAltAllele(const char* s) {
  stringTokenize(s, ",", &altAllele);
}

double VCFGenotypeExtractor::getGenotype(VCFIndividual& indv,
                                         const bool useDosage,
                                         const bool hemiRegion, const int sex,
                                         const int genoIdx, const int GDidx,
                                         const int GQidx) {
  double ret;
  if (genoIdx >= 0) {
    if (useDosage) {
      if (!hemiRegion) {
        ret = (indv.justGet(genoIdx).toDouble());
      } else {
        // for male hemi region, imputated dosage is usually between 0 and 1
        // need to multiply by 2.0
        if (sex == PLINK_MALE) {
          ret = (indv.justGet(genoIdx).toDouble() * 2.0);
        } else {
          ret = (indv.justGet(genoIdx).toDouble());
        }
      }
    } else {  // use hard-coded genotypes
      if (!hemiRegion) {
        ret = (indv.justGet(genoIdx).getGenotype());
      } else {
        if (sex == PLINK_MALE) {
          ret = (indv.justGet(genoIdx).getMaleNonParGenotype02());
        } else if (sex == PLINK_FEMALE) {
          ret = (indv.justGet(genoIdx).getGenotype());
        } else {
          ret = (MISSING_GENOTYPE);
        }
      }
    }
    if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
      return MISSING_GENOTYPE;
    } else {
      return ret;
    }
  } else {
    logger->error("Cannot find %s field!",
                  this->dosageTag.empty() ? "GT" : dosageTag.c_str());
    return MISSING_GENOTYPE;
  }
}

double VCFGenotypeExtractor::getGenotypeForAltAllele(
    VCFIndividual& indv, const bool useDosage, const bool hemiRegion,
    const int sex, const int genoIdx, const int GDidx, const int GQidx,
    const int alt) {
  double ret;
  if (genoIdx >= 0) {
    if (useDosage) {
      if (!hemiRegion) {
        ret = (indv.justGet(genoIdx).toDouble());
      } else {
        // for male hemi region, imputated dosage is usually between 0 and 1
        // need to multiply by 2.0
        if (sex == PLINK_MALE) {
          ret = (indv.justGet(genoIdx).toDouble() * 2.0);
        } else {
          ret = (indv.justGet(genoIdx).toDouble());
        }
      }
    } else {  // use hard-coded genotypes
      if (!hemiRegion) {
        ret = indv.justGet(genoIdx).countAltAllele(alt);
      } else {
        if (sex == PLINK_MALE) {
          // ret = (indv.justGet(genoIdx).getMaleNonParGenotype02());
          ret = (indv.justGet(genoIdx).countMaleNonParAltAllele2(alt));
        } else if (sex == PLINK_FEMALE) {
          // ret = (indv.justGet(genoIdx).getGenotype());
          ret = indv.justGet(genoIdx).countAltAllele(alt);
        } else {
          ret = (MISSING_GENOTYPE);
        }
      }
    }
    if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
      return MISSING_GENOTYPE;
    } else {
      return ret;
    }
  } else {
    logger->error("Cannot find %s field!",
                  this->dosageTag.empty() ? "GT" : dosageTag.c_str());
    return MISSING_GENOTYPE;
  }
}
