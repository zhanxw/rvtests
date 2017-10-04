#include "BGenGenotypeExtractor.h"

#include "GenotypeCounter.h"
#include "Result.h"

#include "base/Argument.h"
#include "base/Logger.h"
#include "libBgen/BGenFile.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/MathVector.h"

DECLARE_BOOL_PARAMETER(outputID);

extern Logger* logger;

BGenGenotypeExtractor::BGenGenotypeExtractor(const std::string& fn)
    : GenotypeExtractor(fn), bgenIn(NULL), altAlleleToParse(-1) {
  this->bgenIn = new BGenFile(fn.c_str());
  if (this->bgenIn->getSampleIdentifier().empty()) {
    fprintf(stderr,
            "BGEN file does not include sample names, please specify sample "
            "file\n");
    exit(1);
  }
}

BGenGenotypeExtractor::BGenGenotypeExtractor(const std::string& fn,
                                             const std::string& sampleFile)
    : GenotypeExtractor(fn), bgenIn(NULL), altAlleleToParse(-1) {
  this->bgenIn = new BGenFile(fn.c_str());

  if (this->bgenIn->loadSampleFile(sampleFile)) {
    fprintf(stderr, "Cannot load sample file [ %s ]!\n", sampleFile.c_str());
    exit(1);
  }
}

BGenGenotypeExtractor::~BGenGenotypeExtractor() {
  if (this->bgenIn) {
    delete this->bgenIn;
    this->bgenIn = NULL;
  }
}

int BGenGenotypeExtractor::extractMultipleGenotype(Matrix* g) {
  assert(g);
  assert(g->rows >= 0 && g->cols >= 0);
  g->Dimension(0, 0);
  int row = 0;
  this->genotype.clear();
  this->altAlleleToParse = -1;
  const BGenVariant& var = this->bgenIn->getVariant();
  const bool useDosage = true;  // TODO: later may support genotype calls
  while (true) {
    if (this->altAlleleToParse <= 0) {
      if (!this->bgenIn->readRecord()) {
        // reached the end
        break;
      } else {
        // const char* alt = this->bgenIn->getBGENRecord().getAlt();
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
    // BGENRecord& r = this->bgenIn->getBGENRecord();
    // BGENPeople& people = r.getPeople();
    // BGENIndividual* indv;

    this->sampleSize = this->bgenIn->getNumEffectiveSample();
    row++;
    this->variantName.resize(row);
    this->counter.resize(row);
    this->counter.back().reset();
    this->hemiRegion.resize(row);

    // // get GT index and cannot assume it is a constant across variants
    // int genoIdx;
    // const bool useDosage = (!this->dosageTag.empty());
    // if (useDosage) {
    //   genoIdx = r.getFormatIndex(dosageTag.c_str());
    // } else {
    //   genoIdx = r.getFormatIndex("GT");
    // }
    // if (useDosage && multiAllelicMode) {
    //   logger->error(
    //       "Unsupported scenario: multiple mode and use dosage as "
    //       "genotypes! - we will only use dosage");
    //   this->altAlleleToParse = 1;
    // }
    // const int GDidx = r.getFormatIndex("GD");
    // const int GQidx = r.getFormatIndex("GQ");
    assert(this->parRegion);
    const bool isHemiRegion =
        this->parRegion->isHemiRegion(var.chrom.c_str(), var.pos);
    // e.g.: Loop each (selected) people in the same order as in the BGEN
    double geno;
    /// const int altAlleleGT = this->altAllele.size() - this->altAlleleToParse;
    const int altAlleleGT =
        var.alleles.size() -
        this->altAlleleToParse;  // when allele = K, alleleGT = 1, 2. ..., K-1
    for (int i = 0; i < sampleSize; i++) {
      // indv = people[i];
      const int indvIdx = this->bgenIn->getEffectiveIndex(i);
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
  }  // end while (this->bgenIn->readRecord())

  // now transpose (marker by people -> people by marker)
  if (row > 0) {
    assert((int)genotype.size() == this->sampleSize * row);
    assign(this->genotype, sampleSize, row, g);
    for (int i = 0; i < row; ++i) {
      g->SetColumnLabel(i, variantName[i].c_str());
    }
  }
  return SUCCEED;
}  // end extractMultipleGenotype(Matrix* g)

int BGenGenotypeExtractor::extractSingleGenotype(Matrix* g, Result* b) {
  this->genotype.clear();
  Result& buf = *b;
  const BGenVariant& var = this->bgenIn->getVariant();
  const bool useDosage = true;

  if (this->altAlleleToParse <= 0) {
    if (!this->bgenIn->readRecord()) {
      return FILE_END;
    } else {
      if (multiAllelicMode) {
        // parseAltAllele(alt);
        this->altAlleleToParse = var.alleles.size() - 1;
      } else {
        // this->altAllele.resize(1);
        // this->altAllele[0] = alt;
        this->altAlleleToParse = 1;
      }
    }
  }
  assert(this->altAlleleToParse >= 0);
  // BGENRecord& r = this->bgenIn->getBGENRecord();
  // BGENPeople& people = r.getPeople();
  // BGENIndividual* indv;

  buf.updateValue("CHROM", var.chrom);
  buf.updateValue("POS", (int)var.pos);
  if (FLAG_outputID) {
    buf.updateValue("ID", var.rsid);
  }
  buf.updateValue("REF", var.alleles[0]);
  buf.updateValue("ALT",
                  var.alleles[var.alleles.size() - this->altAlleleToParse]);

  // genotype.Dimension(people.size(), 1);
  this->sampleSize = this->bgenIn->getNumEffectiveSample();
  this->variantName.resize(1);
  this->counter.resize(1);
  this->counter.back().reset();
  this->hemiRegion.resize(1);

  // // get GT index. if you are sure the index will not change, call this
  // // function only once!
  // const bool useDosage = (!this->dosageTag.empty());
  // int genoIdx;
  // if (useDosage) {
  //   genoIdx = r.getFormatIndex(dosageTag.c_str());
  // } else {
  //   genoIdx = r.getFormatIndex("GT");
  // }
  // const int GDidx = r.getFormatIndex("GD");
  // const int GQidx = r.getFormatIndex("GQ");
  bool isHemiRegion = this->parRegion->isHemiRegion(var.chrom.c_str(), var.pos);
  // e.g.: Loop each (selected) people in the same order as in the BGEN
  double geno;
  // const int altAlleleGT = this->altAllele.size() - this->altAlleleToParse +
  // 1;
  const int altAlleleGT = var.alleles.size() - this->altAlleleToParse;
  for (int i = 0; i < sampleSize; i++) {
    // indv = people[i];
    const int indvIdx = this->bgenIn->getEffectiveIndex(i);
    if (multiAllelicMode) {
      geno = getGenotypeForAltAllele(var, indvIdx, useDosage, isHemiRegion,
                                     (*sex)[i], altAlleleGT);
    } else {
      geno = getGenotype(var, indvIdx, useDosage, isHemiRegion, (*sex)[i]);
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

  variantName.back() = var.chrom;
  variantName.back() += ':';
  variantName.back() += toString(var.pos);
  hemiRegion.back() = isHemiRegion;

  assert((int)genotype.size() == sampleSize);
  assign(genotype, sampleSize, 1, g);
  g->SetColumnLabel(0, variantName.back().c_str());

  --this->altAlleleToParse;
  return SUCCEED;
}  // end extractSingleGenotype()

#if 0
int loadMarkerFromBGEN(const std::string& fileName, const std::string& marker,
                      std::vector<std::string>* rowLabel, Matrix* genotype) {
  if (!rowLabel || !genotype) {
    // invalid parameter
    return -1;
  }
  Matrix& m = *genotype;
  int col = 0;

  BGENInputFile bgenIn(fileName);
  bgenIn.setRangeList(marker);

  while (bgenIn.readRecord()) {
    BGENRecord& r = bgenIn.getBGENRecord();
    BGENPeople& people = r.getPeople();
    BGENIndividual* indv;

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
bool BGenGenotypeExtractor::setSiteFreqMin(const double f) {
  if (f < 0.0 || f > 1.0) {
    return false;
  }
  this->freqMin = f - 1e-10;  // allow rounding error
  return true;
}
bool BGenGenotypeExtractor::setSiteFreqMax(const double f) {
  if (f < 0.0 || f > 1.0) {
    return false;
  }
  this->freqMax = f + 1e-10;  // allow rounding error
  return true;
}

void BGenGenotypeExtractor::setSiteDepthMin(int d) {
  warnUnsupported("site depth");
}
void BGenGenotypeExtractor::setSiteDepthMax(int d) {
  warnUnsupported("site depth");
}

// // @return true if GD is valid
// // if GD is missing, we will take GD = 0
// bool BGenGenotypeExtractor::checkGD(BGENIndividual& indv, int gdIdx) {
//   warnUnsupported("GD");
//   return true;
// }
// bool BGenGenotypeExtractor::checkGQ(BGENIndividual& indv, int gqIdx) {
//   warnUnsupported("GQ");
//   return true;
// }
void BGenGenotypeExtractor::setGDmin(int m) { warnUnsupported("GD"); }
void BGenGenotypeExtractor::setGDmax(int m) { warnUnsupported("GD"); }
void BGenGenotypeExtractor::setGQmin(int m) { warnUnsupported("GQ"); }
void BGenGenotypeExtractor::setGQmax(int m) { warnUnsupported("GQ"); }

void BGenGenotypeExtractor::setSiteFile(const std::string& fn) {
  this->bgenIn->setSiteFile(fn);
}

void BGenGenotypeExtractor::setSiteQualMin(int q) {
  warnUnsupported("SiteQual");
}
void BGenGenotypeExtractor::setSiteMACMin(int n) { warnUnsupported("SiteMAC"); }
int BGenGenotypeExtractor::setAnnoType(const std::string& s) {
  warnUnsupported("anno");
  return -1;
}

void BGenGenotypeExtractor::setRange(const RangeList& l) {
  this->bgenIn->setRange(l);
}
void BGenGenotypeExtractor::setRangeList(const std::string& l) {
  this->bgenIn->setRangeList(l);
}
void BGenGenotypeExtractor::setRangeFile(const std::string& fn) {
  this->bgenIn->setRangeFile(fn.c_str());
}
void BGenGenotypeExtractor::includePeople(const std::string& v) {
  this->bgenIn->includePeople(v.c_str());
}
void BGenGenotypeExtractor::includePeople(const std::vector<std::string>& v) {
  this->bgenIn->includePeople(v);
}
void BGenGenotypeExtractor::includePeopleFromFile(const std::string& fn) {
  this->bgenIn->includePeopleFromFile(fn.c_str());
}
void BGenGenotypeExtractor::excludePeople(const std::string& v) {
  this->bgenIn->excludePeople(v.c_str());
}
void BGenGenotypeExtractor::excludePeople(const std::vector<std::string>& v) {
  this->bgenIn->excludePeople(v);
}
void BGenGenotypeExtractor::excludePeopleFromFile(const std::string& fn) {
  this->bgenIn->excludePeopleFromFile(fn.c_str());
}
void BGenGenotypeExtractor::excludePeople(
    const std::vector<std::string>& sample, const std::vector<int>& index) {
  for (size_t i = 0; i < index.size(); ++i) {
    this->excludePeople(sample[i]);
  }
}

void BGenGenotypeExtractor::excludeAllPeople() {
  this->bgenIn->excludeAllPeople();
}
void BGenGenotypeExtractor::enableAutoMerge() {
  this->bgenIn->enableAutoMerge();
}
void BGenGenotypeExtractor::getPeopleName(std::vector<std::string>* p) {
  *p = this->bgenIn->getSampleIdentifier();
  return;
}

void BGenGenotypeExtractor::getIncludedPeopleName(
    std::vector<std::string>* p) const {
  this->bgenIn->getIncludedSampleName(p);
  return;
}

void BGenGenotypeExtractor::parseAltAllele(const char* s) {
  // TODO
  stringTokenize(s, ",", &altAllele);
}

double BGenGenotypeExtractor::getGenotype(const BGenVariant& var, int indvIdx,
                                          const bool useDosage,
                                          const bool hemiRegion,
                                          const int sex) {
  const int idx = this->bgenIn->getEffectiveIndex(indvIdx);
  double ret;
  if (var.missing[idx]) {
    return MISSING_GENOTYPE;
  }
  if (var.ploidy[idx] == 2) {
    if (hemiRegion && sex == PLINK_MALE) {
      logger->error(
          "Diploid variant detected, but it is in the hemizygous region from a "
          "male");
    }
    if (var.alleles.size() == 2) {
      ret = var.prob[var.index[idx] + 1] + var.prob[var.index[idx] + 2] * 2.0;
    } else if (var.alleles.size() == 1) {
      ret = 2;
    } else {
      double total = var.prob[var.index[idx]] + var.prob[var.index[idx] + 1] +
                     var.prob[var.index[idx] + 2];
      if (total > 0.) {
        ret = (var.prob[var.index[idx] + 1] +
               var.prob[var.index[idx] + 2] * 2.0) /
              total;
      } else {
        ret = MISSING_GENOTYPE;
      }
    }
  } else if (var.ploidy[idx] == 1) {
    if (!hemiRegion) {
      logger->error(
          "Haploid variant detected, but it is not in hemizygous region");
    }
    if (sex == PLINK_FEMALE) {
      logger->error("Haploid variant in detected, but it is from a female");
    }
    if (var.alleles.size() == 2) {
      ret = var.prob[var.index[idx] + 1] + var.prob[var.index[idx] + 2] * 2.0;
    } else if (var.alleles.size() == 1) {
      ret = 2;
    } else {
      double total = var.prob[var.index[idx]] + var.prob[var.index[idx] + 1] +
                     var.prob[var.index[idx] + 2];
      if (total > 0.) {
        ret = (var.prob[var.index[idx] + 1] +
               var.prob[var.index[idx] + 2] * 2.0) /
              total;
      } else {
        ret = MISSING_GENOTYPE;
      }
    }
  } else {
    ret = MISSING_GENOTYPE;
  }
  return ret;
}

double BGenGenotypeExtractor::getGenotypeForAltAllele(
    const BGenVariant& var, int indvIdx, const bool useDosage,
    const bool hemiRegion, const int sex, const int alt) {
  if (alt > 1) {
    // TODO, better support for multi-allelic sites
    return MISSING_GENOTYPE;
  }
  return getGenotype(var, indvIdx, useDosage, hemiRegion, sex);
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

void BGenGenotypeExtractor::warnUnsupported(const char* tag) {
  logger->warn("Please remove unsupported features related to %s", tag);
}
