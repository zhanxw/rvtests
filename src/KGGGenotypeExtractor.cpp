#include "KGGGenotypeExtractor.h"

#include "GenotypeCounter.h"
#include "Result.h"

#include "base/Argument.h"
#include "base/Logger.h"
#include "libVcf/KGGInputFile.h"
#include "libVcf/VCFConstant.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/MathVector.h"

DECLARE_BOOL_PARAMETER(outputID);

extern Logger* logger;

KGGGenotypeExtractor::KGGGenotypeExtractor(const std::string& fn)
    : GenotypeExtractor(fn),
      kggIn(NULL),
      altAlleleToParse(-1),
      currentVariant(0) {
  this->kggIn = new KGGInputFile(fn.c_str(), ".gz");
}

KGGGenotypeExtractor::~KGGGenotypeExtractor() {
  if (this->kggIn) {
    delete this->kggIn;
    this->kggIn = NULL;
  }
}

int KGGGenotypeExtractor::extractMultipleGenotype(Matrix* g) {
  warnUnsupported("extractMultipleGenotype");
  return -1;
#if 0
  assert(g);
  assert(g->rows >= 0 && g->cols >= 0);
  g->Dimension(0, 0);
  int row = 0;
  this->genotype.clear();
  this->altAlleleToParse = -1;
  const bool useDosage = false;  // TODO: later may support genotype calls
  while (true) {
    if (this->altAlleleToParse <= 0) {
      if (!this->kggIn->readRecord()) {
        // reached the end
        break;
      } else {
        // const char* alt = this->kggIn->getKGGRecord().getAlt();
        if (multiAllelicMode) {
          // parseAltAllele(alt);
          // this->altAlleleToParse = altAllele.size();
          this->altAlleleToParse =
              var.alleles.size() -
              1;  // num of alternative alleles to be parsed/processed
        } else {
          // this->altAllele.resize(1);
          // this->altAllele[0] = alt;
          this->altAlleleToParse = 1;
        }
      }
    }
    assert(this->altAlleleToParse > 0);
    // parse alt allele one at a time

    this->sampleSize = this->kggIn->getNumEffectiveSample();
    row++;
    this->variantName.resize(row);
    this->counter.resize(row);
    this->counter.back().reset();
    this->hemiRegion.resize(row);

    assert(this->parRegion);
    const bool isHemiRegion =
        this->parRegion->isHemiRegion(var.chrom.c_str(), var.pos);
    // e.g.: Loop each (selected) people in the same order as in the KGG
    double geno;
    /// const int altAlleleGT = this->altAllele.size() - this->altAlleleToParse;
    const int altAlleleGT =
        var.alleles.size() -
        this->altAlleleToParse;  // when allele = K, alleleGT = 1, 2. ..., K-1
    for (int i = 0; i < sampleSize; i++) {
      // indv = people[i];
      const int indvIdx = this->kggIn->getEffectiveIndex(i);
      if (multiAllelicMode) {
        geno = getGenotypeForAltAllele(var, indvIdx, useDosage, isHemiRegion,
                                       (*sex)[i], altAlleleGT);
      } else {
        geno = getGenotype(var, indvIdx, useDosage, isHemiRegion, (*sex)[i]);
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

    this->variantName.back() = var.chrom;
    this->variantName.back() += ":";
    this->variantName.back() += var.pos;
    if (multiAllelicMode) {
      this->variantName.back() += var.alleles[0];
      this->variantName.back() += "/";
      this->variantName.back() +=
          var.alleles[altAllele.size() - this->altAlleleToParse];
    }
    this->hemiRegion.back() = (isHemiRegion);

    this->altAlleleToParse--;
  }  // end while (this->kggIn->readRecord())

  // now transpose (marker by people -> people by marker)
  if (row > 0) {
    assert((int)genotype.size() == this->sampleSize * row);
    assign(this->genotype, sampleSize, row, g);
    for (int i = 0; i < row; ++i) {
      g->SetColumnLabel(i, variantName[i].c_str());
    }
  }
  return SUCCEED;
#endif
}  // end extractMultipleGenotype(Matrix* g)

int KGGGenotypeExtractor::extractSingleGenotype(Matrix* g, Result* b) {
  this->genotype.clear();
  Result& buf = *b;
  const bool useDosage = true;
  int numAltAllele = this->kggIn->getAlt()[currentVariant].size();
  if (this->altAlleleToParse <= 0) {
    if (!this->kggIn->readRecord()) {
      return FILE_END;
    } else {
      if (multiAllelicMode) {
        // parseAltAllele(alt);
        this->altAlleleToParse = numAltAllele - 1;
      } else {
        // this->altAllele.resize(1);
        // this->altAllele[0] = alt;
        this->altAlleleToParse = 1;
      }
    }
  }
  assert(this->altAlleleToParse >= 0);

  buf.updateValue("CHROM", this->kggIn->getChrom()[currentVariant]);
  buf.updateValue("POS", this->kggIn->getPosition()[currentVariant]);
  if (FLAG_outputID) {
    buf.updateValue("ID", this->kggIn->getMarkerName()[currentVariant]);
  }
  buf.updateValue("REF", this->kggIn->getRef()[currentVariant]);
  buf.updateValue("ALT",
                  this->kggIn->getAlt()[currentVariant]
                                       [numAltAllele - this->altAlleleToParse]);

  // genotype.Dimension(people.size(), 1);
  this->sampleSize = this->kggIn->getNumEffectiveSample();
  this->variantName.resize(1);
  this->counter.resize(1);
  this->counter.back().reset();
  this->hemiRegion.resize(1);

  bool isHemiRegion = this->parRegion->isHemiRegion(
      this->kggIn->getChrom()[currentVariant].c_str(),
      this->kggIn->getPosition()[currentVariant]);
  // e.g.: Loop each (selected) people in the same order as in the KGG
  double geno;
  // const int altAlleleGT = this->altAllele.size() - this->altAlleleToParse +
  // 1;
  const int altAlleleGT = numAltAllele - this->altAlleleToParse;
  for (int i = 0; i < sampleSize; i++) {
    // indv = people[i];
    const int indvIdx = this->kggIn->getEffectiveIndex(i);
    if (multiAllelicMode) {
      geno = getGenotypeForAltAllele(indvIdx, useDosage, isHemiRegion,
                                     (*sex)[i], altAlleleGT);
    } else {
      geno = getGenotype(indvIdx, useDosage, isHemiRegion, (*sex)[i]);
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

  variantName.back() = this->kggIn->getChrom()[currentVariant];
  variantName.back() += ':';
  variantName.back() += toString(this->kggIn->getPosition()[currentVariant]);
  hemiRegion.back() = isHemiRegion;

  assert((int)genotype.size() == sampleSize);
  assign(genotype, sampleSize, 1, g);
  g->SetColumnLabel(0, variantName.back().c_str());

  --this->altAlleleToParse;
  if (this->altAlleleToParse == 0) {
    ++currentVariant;
  }
  return SUCCEED;
}  // end extractSingleGenotype()

#if 0
int loadMarkerFromKGG(const std::string& fileName, const std::string& marker,
                      std::vector<std::string>* rowLabel, Matrix* genotype) {
  if (!rowLabel || !genotype) {
    // invalid parameter
    return -1;
  }
  Matrix& m = *genotype;
  int col = 0;

  KGGInputFile kggIn(fileName);
  kggIn.setRangeList(marker);

  while (kggIn.readRecord()) {
    KGGRecord& r = kggIn.getKGGRecord();
    KGGPeople& people = r.getPeople();
    KGGIndividual* indv;

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
#endif
bool KGGGenotypeExtractor::setSiteFreqMin(const double f) {
  if (f < 0.0 || f > 1.0) {
    return false;
  }
  this->freqMin = f - 1e-10;  // allow rounding error
  return true;
}
bool KGGGenotypeExtractor::setSiteFreqMax(const double f) {
  if (f < 0.0 || f > 1.0) {
    return false;
  }
  this->freqMax = f + 1e-10;  // allow rounding error
  return true;
}

void KGGGenotypeExtractor::setSiteDepthMin(int d) {
  warnUnsupported("site depth");
}
void KGGGenotypeExtractor::setSiteDepthMax(int d) {
  warnUnsupported("site depth");
}

// // @return true if GD is valid
// // if GD is missing, we will take GD = 0
// bool KGGGenotypeExtractor::checkGD(KGGIndividual& indv, int gdIdx) {
//   warnUnsupported("GD");
//   return true;
// }
// bool KGGGenotypeExtractor::checkGQ(KGGIndividual& indv, int gqIdx) {
//   warnUnsupported("GQ");
//   return true;
// }
void KGGGenotypeExtractor::setGDmin(int m) { warnUnsupported("GD"); }
void KGGGenotypeExtractor::setGDmax(int m) { warnUnsupported("GD"); }
void KGGGenotypeExtractor::setGQmin(int m) { warnUnsupported("GQ"); }
void KGGGenotypeExtractor::setGQmax(int m) { warnUnsupported("GQ"); }

void KGGGenotypeExtractor::setSiteFile(const std::string& fn) {
  this->kggIn->setSiteFile(fn);
}

void KGGGenotypeExtractor::setSiteQualMin(int q) {
  warnUnsupported("SiteQual");
}
void KGGGenotypeExtractor::setSiteMACMin(int n) { warnUnsupported("SiteMAC"); }
int KGGGenotypeExtractor::setAnnoType(const std::string& s) {
  warnUnsupported("anno");
  return -1;
}

void KGGGenotypeExtractor::setRange(const RangeList& l) {
  warnUnsupported("setRange");
}
void KGGGenotypeExtractor::setRangeList(const std::string& l) {
  warnUnsupported("setRangeList");
}
void KGGGenotypeExtractor::setRangeFile(const std::string& fn) {
  warnUnsupported("setRangeFile");
}
void KGGGenotypeExtractor::includePeople(const std::string& v) {
  this->kggIn->includePeople(v.c_str());
}
void KGGGenotypeExtractor::includePeople(const std::vector<std::string>& v) {
  this->kggIn->includePeople(v);
}
void KGGGenotypeExtractor::includePeopleFromFile(const std::string& fn) {
  this->kggIn->includePeopleFromFile(fn.c_str());
}
void KGGGenotypeExtractor::excludePeople(const std::string& v) {
  this->kggIn->excludePeople(v.c_str());
}
void KGGGenotypeExtractor::excludePeople(const std::vector<std::string>& v) {
  this->kggIn->excludePeople(v);
}
void KGGGenotypeExtractor::excludePeopleFromFile(const std::string& fn) {
  this->kggIn->excludePeopleFromFile(fn.c_str());
}
void KGGGenotypeExtractor::excludePeople(const std::vector<std::string>& sample,
                                         const std::vector<int>& index) {
  for (size_t i = 0; i < index.size(); ++i) {
    this->excludePeople(sample[i]);
  }
}

void KGGGenotypeExtractor::excludeAllPeople() {
  this->kggIn->excludeAllPeople();
}
void KGGGenotypeExtractor::enableAutoMerge() {
  warnUnsupported("enableAutoMerge");
}
void KGGGenotypeExtractor::getPeopleName(std::vector<std::string>* p) {
  *p = this->kggIn->getSampleName();
  return;
}

void KGGGenotypeExtractor::getIncludedPeopleName(
    std::vector<std::string>* p) const {
  this->kggIn->getIncludedSampleName(p);
  return;
}

void KGGGenotypeExtractor::parseAltAllele(const char* s) {
  stringTokenize(s, ",", &altAllele);
}

double KGGGenotypeExtractor::getGenotype(int indvIdx, const bool useDosage,
                                         const bool hemiRegion, const int sex) {
  const int idx = this->kggIn->getEffectiveIndex(indvIdx);
  double ret;
  assert(!useDosage);
  int a1, a2;
  this->kggIn->getAllele(idx, &a1, &a2);
  if (a1 < 0 || a2 < 0) {
    return MISSING_GENOTYPE;
  }

  if (hemiRegion && sex == PLINK_MALE && a1 != a2) {
    logger->error(
        "In hemizygous region, male have different alleles, treated as missing "
        "genotypes");
    return MISSING_GENOTYPE;
  }

  if (a1 > 1) {
    a1 = 1;
  }
  if (a2 > 1) {
    a2 = 1;
  }
  return a1 + a2;
  return ret;
}

double KGGGenotypeExtractor::getGenotypeForAltAllele(int indvIdx,
                                                     const bool useDosage,
                                                     const bool hemiRegion,
                                                     const int sex,
                                                     const int alt) {
  const int idx = this->kggIn->getEffectiveIndex(indvIdx);
  assert(!useDosage);
  int a1, a2;
  this->kggIn->getAllele(idx, &a1, &a2);
  if (a1 < 0 || a2 < 0) {
    return MISSING_GENOTYPE;
  }

  if (hemiRegion && sex == PLINK_MALE && a1 != a2) {
    logger->error(
        "In hemizygous region, male have different alleles, treated as missing "
        "genotypes");
    return MISSING_GENOTYPE;
  }

  return (alt == a1 ? 1 : 0) + (alt == a2 ? 1 : 0);
}
#if 0
// TODO
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
#endif

void KGGGenotypeExtractor::warnUnsupported(const char* tag) {
  logger->warn("Please remove unsupported features related to %s", tag);
}
