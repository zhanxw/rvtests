#include "GenotypeExtractor.h"

#include "GenotypeCounter.h"
#include "Result.h"

#include "base/Argument.h"
#include "base/Logger.h"
#include "libVcf/VCFUtil.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/MathVector.h"

DECLARE_BOOL_PARAMETER(outputID);

extern Logger* logger;

GenotypeExtractor::GenotypeExtractor(const std::string& fn)
    : freqMin(-1),
      freqMax(-1),
      GDmin(-1),
      GDmax(-1),
      needGD(false),
      GQmin(-1),
      GQmax(-1),
      needGQ(false),
      parRegion(NULL),
      sex(NULL),
      sampleSize(-1),
      multiAllelicMode(false) {}

GenotypeExtractor::~GenotypeExtractor() {}

// bool GenotypeExtractor::setSiteFreqMin(const double f) {
//   if (f < 0.0 || f > 1.0) {
//     return false;
//   }
//   this->freqMin = f - 1e-10;  // allow rounding error
//   return true;
// }
// bool GenotypeExtractor::setSiteFreqMax(const double f) {
//   if (f < 0.0 || f > 1.0) {
//     return false;
//   }
//   this->freqMax = f + 1e-10;  // allow rounding error
//   return true;
// }

// void GenotypeExtractor::setSiteDepthMin(int d) {
//   this->vin->setSiteDepthMin(d);
// }
// void GenotypeExtractor::setSiteDepthMax(int d) {
//   this->vin->setSiteDepthMax(d);
// }

// // @return true if GD is valid
// // if GD is missing, we will take GD = 0
// bool GenotypeExtractor::checkGD(VCFIndividual& indv, int gdIdx) {
//   if (!needGD) return true;
//   int gd = indv.justGet(gdIdx).toInt();
//   if (this->GDmin > 0 && gd < this->GDmin) return false;
//   if (this->GDmax > 0 && gd > this->GDmax) return false;
//   return true;
// }
// bool GenotypeExtractor::checkGQ(VCFIndividual& indv, int gqIdx) {
//   if (!needGQ) return true;
//   int gq = indv.justGet(gqIdx).toInt();
//   if (this->GQmin > 0 && gq < this->GQmin) return false;
//   if (this->GQmax > 0 && gq > this->GQmax) return false;
//   return true;
// }
// void GenotypeExtractor::setGDmin(int m) {
//   this->needGD = true;
//   this->GDmin = m;
// }
// void GenotypeExtractor::setGDmax(int m) {
//   this->needGD = true;
//   this->GDmax = m;
// }
// void GenotypeExtractor::setGQmin(int m) {
//   this->needGQ = true;
//   this->GQmin = m;
// }
// void GenotypeExtractor::setGQmax(int m) {
//   this->needGQ = true;
//   this->GQmax = m;
// }

// void GenotypeExtractor::setSiteFile(const std::string& fn) {
//   this->vin->setSiteFile(fn);
// }

// void GenotypeExtractor::setSiteQualMin(int q) { this->vin->setSiteQualMin(q);
// }
// void GenotypeExtractor::setSiteMACMin(int n) { this->vin->setSiteMACMin(n); }
// int GenotypeExtractor::setAnnoType(const std::string& s) {
//   return this->vin->setAnnoType(s.c_str());
// }

// void GenotypeExtractor::setRange(const RangeList& l) {
// this->vin->setRange(l); }
// void GenotypeExtractor::setRangeList(const std::string& l) {
//   this->vin->setRangeList(l);
// }
// void GenotypeExtractor::setRangeFile(const std::string& fn) {
//   this->vin->setRangeFile(fn.c_str());
// }
// void GenotypeExtractor::includePeople(const std::string& v) {
//   this->vin->includePeople(v.c_str());
// }
// void GenotypeExtractor::includePeople(const std::vector<std::string>& v) {
//   this->vin->includePeople(v);
// }
// void GenotypeExtractor::includePeopleFromFile(const std::string& fn) {
//   this->vin->includePeopleFromFile(fn.c_str());
// }
// void GenotypeExtractor::excludePeople(const std::string& v) {
//   this->vin->excludePeople(v.c_str());
// }
// void GenotypeExtractor::excludePeople(const std::vector<std::string>& v) {
//   this->vin->excludePeople(v);
// }
// void GenotypeExtractor::excludePeopleFromFile(const std::string& fn) {
//   this->vin->excludePeopleFromFile(fn.c_str());
// }
// void GenotypeExtractor::excludePeople(const std::vector<std::string>& sample,
//                                       const std::vector<int>& index) {
//   for (size_t i = 0; i < index.size(); ++i) {
//     this->excludePeople(sample[i]);
//   }
// }

// void GenotypeExtractor::excludeAllPeople() { this->vin->excludeAllPeople(); }
// void GenotypeExtractor::enableAutoMerge() { this->vin->enableAutoMerge(); }
// void GenotypeExtractor::getPeopleName(std::vector<std::string>* p) {
//   return this->vin->getVCFHeader()->getPeopleName(p);
// }

// void GenotypeExtractor::getIncludedPeopleName(
//     std::vector<std::string>* p) const {
//   this->vin->getIncludedPeopleName(p);
//   return;
// }

// void GenotypeExtractor::parseAltAllele(const char* s) {
//   stringTokenize(s, ",", &altAllele);
// }

// double GenotypeExtractor::getGenotype(VCFIndividual& indv, const bool
// useDosage,
//                                       const bool hemiRegion, const int sex,
//                                       const int genoIdx, const int GDidx,
//                                       const int GQidx) {
//   double ret;
//   if (genoIdx >= 0) {
//     if (useDosage) {
//       if (!hemiRegion) {
//         ret = (indv.justGet(genoIdx).toDouble());
//       } else {
//         // for male hemi region, imputated dosage is usually between 0 and 1
//         // need to multiply by 2.0
//         if (sex == PLINK_MALE) {
//           ret = (indv.justGet(genoIdx).toDouble() * 2.0);
//         } else {
//           ret = (indv.justGet(genoIdx).toDouble());
//         }
//       }
//     } else {  // use hard-coded genotypes
//       if (!hemiRegion) {
//         ret = (indv.justGet(genoIdx).getGenotype());
//       } else {
//         if (sex == PLINK_MALE) {
//           ret = (indv.justGet(genoIdx).getMaleNonParGenotype02());
//         } else if (sex == PLINK_FEMALE) {
//           ret = (indv.justGet(genoIdx).getGenotype());
//         } else {
//           ret = (MISSING_GENOTYPE);
//         }
//       }
//     }
//     if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
//       return MISSING_GENOTYPE;
//     } else {
//       return ret;
//     }
//   } else {
//     logger->error("Cannot find %s field!",
//                   this->dosageTag.empty() ? "GT" : dosageTag.c_str());
//     return MISSING_GENOTYPE;
//   }
// }

// double GenotypeExtractor::getGenotypeForAltAllele(
//     VCFIndividual& indv, const bool useDosage, const bool hemiRegion,
//     const int sex, const int genoIdx, const int GDidx, const int GQidx,
//     const int alt) {
//   double ret;
//   if (genoIdx >= 0) {
//     if (useDosage) {
//       if (!hemiRegion) {
//         ret = (indv.justGet(genoIdx).toDouble());
//       } else {
//         // for male hemi region, imputated dosage is usually between 0 and 1
//         // need to multiply by 2.0
//         if (sex == PLINK_MALE) {
//           ret = (indv.justGet(genoIdx).toDouble() * 2.0);
//         } else {
//           ret = (indv.justGet(genoIdx).toDouble());
//         }
//       }
//     } else {  // use hard-coded genotypes
//       if (!hemiRegion) {
//         ret = indv.justGet(genoIdx).countAltAllele(alt);
//       } else {
//         if (sex == PLINK_MALE) {
//           // ret = (indv.justGet(genoIdx).getMaleNonParGenotype02());
//           ret = (indv.justGet(genoIdx).countMaleNonParAltAllele2(alt));
//         } else if (sex == PLINK_FEMALE) {
//           // ret = (indv.justGet(genoIdx).getGenotype());
//           ret = indv.justGet(genoIdx).countAltAllele(alt);
//         } else {
//           ret = (MISSING_GENOTYPE);
//         }
//       }
//     }
//     if (!checkGD(indv, GDidx) || !checkGQ(indv, GQidx)) {
//       return MISSING_GENOTYPE;
//     } else {
//       return ret;
//     }
//   } else {
//     logger->error("Cannot find %s field!",
//                   this->dosageTag.empty() ? "GT" : dosageTag.c_str());
//     return MISSING_GENOTYPE;
//   }
// }

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
