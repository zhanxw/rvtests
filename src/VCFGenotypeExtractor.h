#ifndef VCFGENOTYPEEXTRACTOR_H
#define VCFGENOTYPEEXTRACTOR_H

#include <string>
#include <vector>

#include "src/GenotypeExtractor.h"
// #include "base/ParRegion.h"
// #include "libsrc/MathMatrix.h"
// #include "libsrc/MathVector.h"

// class Matrix;
// class RangeList;
// class Result;
class VCFExtractor;
class VCFIndividual;
// class GenotypeCounter;

/**
 * Extract genotype from file @param fileName at the marker @param marker,
 * and store sample names in @param rowLabel and genotypes in @param genotype
 * (dimension is numSample x 1)
 * @return 0 if succeed
 */
int loadMarkerFromVCF(const std::string& fileName, const std::string& marker,
                      std::vector<std::string>* rowLabel, Matrix* genotype);

class VCFGenotypeExtractor : public GenotypeExtractor {
 public:
  explicit VCFGenotypeExtractor(const std::string& fn);
  virtual ~VCFGenotypeExtractor();

 private:
  VCFGenotypeExtractor(const VCFGenotypeExtractor&);
  VCFGenotypeExtractor& operator=(const VCFGenotypeExtractor&);

 public:
  /**
   * @param g, store people by marker matrix
   * @return 0 for success
   */
  int extractMultipleGenotype(Matrix* g);
  /**
   * @return 0 for success
   * @return -2 for reach end.
   * @param g: people by 1 matrix, where column name is like "chr:pos"
   * @param b: extract information, e.g. "1\t100\tA\tC"
   */
  int extractSingleGenotype(Matrix* g, Result* b);

  /* Site filters */
  bool setSiteFreqMin(const double f);
  bool setSiteFreqMax(const double f);
  void setSiteDepthMin(int d);
  void setSiteDepthMax(int d);
  // @return true if GD is valid
  // if GD is missing, we will take GD = 0
  bool checkGD(VCFIndividual& indv, int gdIdx);
  bool checkGQ(VCFIndividual& indv, int gqIdx);
  void setGDmin(int m);
  void setGDmax(int m);
  void setGQmin(int m);
  void setGQmax(int m);

  void setSiteFile(const std::string& fn);
  void setSiteQualMin(int q);
  void setSiteMACMin(int n);
  int setAnnoType(const std::string& s);

  void setRange(const RangeList& l);
  void setRangeList(const std::string& l);
  void setRangeFile(const std::string& fn);
  void includePeople(const std::string& v);
  void includePeople(const std::vector<std::string>& v);
  void includePeopleFromFile(const std::string& fn);
  void excludePeople(const std::string& v);
  void excludePeopleFromFile(const std::string& fn);
  void excludePeople(const std::vector<std::string>& sample);
  void excludePeople(const std::vector<std::string>& sample,
                     const std::vector<int>& index);
  void excludeAllPeople();
  void enableAutoMerge();
  void getPeopleName(std::vector<std::string>* p);
  void getIncludedPeopleName(std::vector<std::string>* p) const;

  /**
   * @return weigth, its length equals to # of markers
   */
  // std::vector<double>& getWeight() { return this->weight; };
  void setDosageTag(const std::string& tag) {
    if (tag.empty()) return;
    this->dosageTag = tag;
  }
  void unsetDosageTag() { this->dosageTag.clear(); }
  // bool isDosage() const { return !this->dosageTag.empty(); }
  // void setParRegion(ParRegion* p) { this->parRegion = p; }
  // //      Sex (1=male; 2=female; other=unknown)
  // void setSex(const std::vector<int>* sex) { this->sex = sex; }
  // coding male chromX as 0/2 instead of 0/1
  // similarly, for dosage, just multiply 2.0 from original dosage
  // void enableClaytonCoding() { this->claytonCoding = true; }
  // void disableClaytonCoding() { this->claytonCoding = false; }

  // check how many alt alleles at this site
  void parseAltAllele(const char* s);
  // extract genotype for @param indv
  inline double getGenotype(VCFIndividual& indv, const bool useDosage,
                            const bool hemiRegion, const int sex,
                            const int genoIdx, const int GDidx,
                            const int GQidx);
  double getGenotypeForAltAllele(VCFIndividual& indv, const bool useDosage,
                                 const bool hemiRegion, const int sex,
                                 const int genoIdx, const int GDidx,
                                 const int GQidx, const int alt);

  // assign extracted genotype @param from to a @param nrow by @param ncol
  // output matrix @param to
  // void assign(const std::vector<double>& from, int nrow, int ncol, Matrix*
  // to);
  // void enableMultiAllelicMode() { this->multiAllelicMode = true; }

 private:
  VCFExtractor* vin;
  std::vector<std::string> altAllele;  // store alt alleles
  int altAlleleToParse;                // number of alleles to parse
};                                     // class VCFGenotypeExtractor

#endif /* VCFGENOTYPEEXTRACTOR_H */
