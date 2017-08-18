#ifndef _VCFFILTER_H_
#define _VCFFILTER_H_

#include "base/Regex.h"

class ParRegion;
/**
 * VCF library support the following filters (in additional to sampler filter
 *and range filter).
 *
 */
class VCFSiteFilter {
 public:
  // any region, PAR regions, and hemizygote regions
  enum ChromXExtraction { ANY = 0, PAR = 1, HEMI = 2 };

 public:
  VCFSiteFilter();
  virtual ~VCFSiteFilter();

  // setter function
  void setUseSiteDepthFromInfo() { this->siteDepthFromInfo = true; };
  void setSiteDepthMin(int d) { this->siteDepthMin = d; }
  void setSiteDepthMax(int d) { this->siteDepthMax = d; }
  void setSiteQualMin(int q) { this->siteQualMin = q; }
  void setUseSiteFreqFromInfo() { this->siteFreqFromInfo = true; };
  void setSiteFreqMin(double f) { this->siteFreqMin = f; };
  void setSiteFreqMax(double f) { this->siteFreqMax = f; };
  void setSiteMACMin(int n) {  // minar allele count
    this->siteMACMin = n;
  };
  int setAnnoType(const char* s) { return this->annoRegex.readPattern(s); };
  void setVariantSiteOnly() { this->onlyVariantSite = true; };
#if 0
  void setIndividualDepthMin(int d) {
    this->indvDepthMin = d;
  };
  void setIndividualDepthMax(int d) {
    this->indvDepthMax = d;
  };
  void setIndividualQualMin(int q) {
    this->indvQualMin = q;
  };
#endif

  // checking whether to check functions
  bool checkSiteDepth() const {
    return (this->siteDepthMin > 0 || this->siteDepthMax > 0);
  };
  bool checkSiteQual() const { return (this->siteQualMin > 0); };
  bool checkSiteFreq() const {
    return (this->siteFreqMin > 0 || this->siteFreqMax > 0);
  };
  bool checkSiteMAC() const { return (this->siteMACMin > 0); };
  // checking functions for each criterias
  bool useSiteDepthFromInfo() const { return this->siteDepthFromInfo; }
  bool siteDepthOK(int d) const {
    if (this->siteDepthMin > 0 && this->siteDepthMin > d) return false;
    if (this->siteDepthMax > 0 && this->siteDepthMax < d) return false;
    return true;
  };
  bool siteQualOK(int q) const {
    if (this->siteQualMin > 0 && this->siteQualMin > q) return false;
    return true;
  };
  bool useSiteFreqFromInfo() const { return this->siteFreqFromInfo; };
  bool siteFreqOK(double f) const {
    if (f > 0.5) {
      f = 1.0 - f;
    };
    if (this->siteFreqMin > 0 && this->siteFreqMin > f) return false;
    if (this->siteFreqMax > 0 && this->siteFreqMax < f) return false;
    return true;
  };
  bool siteMACOK(int n) const {
    if (this->siteMACMin > 0 && this->siteMACMin > n) return false;
    return true;
  };

  bool requiredAnnotation() const { return this->annoRegex.isInitialized(); }
  bool matchAnnotatoin(const char* s) { return this->annoRegex.match(s); }
  bool isVariantSiteOnly() const { return this->onlyVariantSite; };

  void setParRegion(ParRegion* p) { this->parRegion = p; }
  void setExtractChromXParRegion() { this->chromXExtraction = PAR; }

  void setExtractChromXHemiRegion() { this->chromXExtraction = HEMI; }

  void resetExtractChromXRegion() { this->chromXExtraction = ANY; }
#if 0
  bool individualDepthOK(int d) const {
    if (this->indvDepthMin > 0  && this->indvDepthMin > d) return false;
    if (this->indvDepthMax > 0  && this->indvDepthMax < d) return false;
    return true;
  };
  bool individualQualOK(double q) const {
    if ( this->indvQualMin > 0 && this->indvQualMin <= q) return false;
    return true;
  };
#endif
 private:
  // thresholds
  bool siteDepthFromInfo;
  int siteDepthMin;
  int siteDepthMax;
  int siteQualMin;

  bool siteFreqFromInfo;
  double siteFreqMin;
  double siteFreqMax;
  int siteMACMin;

  Regex annoRegex;       // for filter ANNO
  bool onlyVariantSite;  // only extract sites that are polymorphism

#if 0
  int    indvDepthMin;
  int    indvDepthMax;
  int    indvQualMin;
#endif

 protected:
  // allow X-chromosome analysis
  ParRegion* parRegion;
  enum ChromXExtraction chromXExtraction;
};

#endif /* _VCFFILTER_H_ */
