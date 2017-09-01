#ifndef GENOTYPEEXTRACTOR_H
#define GENOTYPEEXTRACTOR_H

#include <string>
#include <vector>

#include "base/ParRegion.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/MathVector.h"

class Matrix;
class RangeList;
class Result;
// class VCFExtractor;
// class VCFIndividual;
class GenotypeCounter;

class GenotypeExtractor {
 public:
  explicit GenotypeExtractor(const std::string& fn);
  virtual ~GenotypeExtractor();

 private:
  GenotypeExtractor(const GenotypeExtractor&);
  GenotypeExtractor& operator=(const GenotypeExtractor&);

 public:
  /**
   * @param g, store people by marker matrix
   * @return 0 for success
   */
  virtual int extractMultipleGenotype(Matrix* g) = 0;
  /**
   * @return 0 for success
   * @return -2 for reach end.
   * @param g: people by 1 matrix, where column name is like "chr:pos"
   * @param b: extract information, e.g. "1\t100\tA\tC"
   */
  virtual int extractSingleGenotype(Matrix* g, Result* b) = 0;

  /* Site filters */
  virtual bool setSiteFreqMin(const double f) = 0;
  virtual bool setSiteFreqMax(const double f) = 0;
  virtual void setSiteDepthMin(int d) = 0;
  virtual void setSiteDepthMax(int d) = 0;
  // @return true if GD is valid
  // if GD is missing, we will take GD = 0
  virtual void setGDmin(int m) = 0;
  virtual void setGDmax(int m) = 0;
  virtual void setGQmin(int m) = 0;
  virtual void setGQmax(int m) = 0;

  virtual void setSiteFile(const std::string& fn) = 0;
  virtual void setSiteQualMin(int q) = 0;
  virtual void setSiteMACMin(int n) = 0;
  virtual int setAnnoType(const std::string& s) = 0;

  virtual void setRange(const RangeList& l) = 0;
  virtual void setRangeList(const std::string& l) = 0;
  virtual void setRangeFile(const std::string& fn) = 0;
  virtual void includePeople(const std::string& v) = 0;
  virtual void includePeople(const std::vector<std::string>& v) = 0;
  virtual void includePeopleFromFile(const std::string& fn) = 0;
  virtual void excludePeople(const std::string& v) = 0;
  virtual void excludePeopleFromFile(const std::string& fn) = 0;
  virtual void excludePeople(const std::vector<std::string>& sample) = 0;
  virtual void excludePeople(const std::vector<std::string>& sample,
                             const std::vector<int>& index) = 0;
  virtual void excludeAllPeople() = 0;
  virtual void enableAutoMerge() = 0;
  virtual void getPeopleName(std::vector<std::string>* p) = 0;
  virtual void getIncludedPeopleName(std::vector<std::string>* p) const = 0;

  const std::vector<GenotypeCounter>& getGenotypeCounter() const {
    return this->counter;
  }
  /**
   * @return weigth, its length equals to # of markers
   */
  // std::vector<double>& getWeight() { return this->weight; };
  void setDosageTag(const std::string& tag) {
    if (tag.empty()) return;
    this->dosageTag = tag;
  }
  void unsetDosageTag() { this->dosageTag.clear(); }
  bool isDosage() const { return !this->dosageTag.empty(); }
  void setParRegion(ParRegion* p) { this->parRegion = p; }
  //      Sex (1=male; 2=female; other=unknown)
  void setSex(const std::vector<int>* sex) { this->sex = sex; }
// coding male chromX as 0/2 instead of 0/1
// similarly, for dosage, just multiply 2.0 from original dosage
// void enableClaytonCoding() { this->claytonCoding = true; }
// void disableClaytonCoding() { this->claytonCoding = false; }

#if 0
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
#endif
  void assign(const std::vector<double>& from, int nrow, int ncol, Matrix* to);
  void enableMultiAllelicMode() { this->multiAllelicMode = true; }

 public:
  const static int SUCCEED = 0;
  const static int ERROR = -1;
  const static int FILE_END = -2;
  const static int FAIL_FILTER = -3;

 protected:
  // VCFExtractor* vin;
  double freqMin;
  double freqMax;
  int GDmin;
  int GDmax;
  bool needGD;
  int GQmin;
  int GQmax;
  bool needGQ;
  // std::vector<double> weight;   // per-variant weight
  // std::vector<double> af;       // per-variant alt allele freq
  // std::vector<int> numMissing;  // per-variant # of samples having missing
  // genotypes
  std::vector<GenotypeCounter> counter;
  std::string dosageTag;  // set if loading dosage instead of genotype

  // compensate sex chromosome
  ParRegion* parRegion;
  std::vector<bool>
      hemiRegion;               // true: if the extracted variant in hemi region
  const std::vector<int>* sex;  // external sex information
  // bool claytonCoding;  // code male hemi region genotype from 0/1 to 0/2
  std::vector<double> genotype;          // store extracted genotypes
  std::vector<std::string> variantName;  // store extracted variant names
  int sampleSize;                        // number of extracted vcf samples
  // for multiallelic
  bool multiAllelicMode;  // default is false
#if 0
  std::vector<std::string> altAllele;  // store alt alleles
  int altAlleleToParse;                // number of alleles to parse
#endif
};  // class GenotypeExtractor

#endif /* GENOTYPEEXTRACTOR_H */
