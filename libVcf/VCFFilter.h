#ifndef _VCFFILTER_H_
#define _VCFFILTER_H_

#include "Regex.h"

/**
 * Our VCF library support the following filters (in additional to sampler filter and range filter).
 *
 */
class VCFFilter{
public:
VCFFilter():
  // Site filter
  depthFromInfo(false),      // read depth from INFO field
      siteDepthMin(-1),
      siteDepthMax(-1),
      siteQualMin(-1),

      freqFromInfo(false),      // read AF from INFO field
      siteFreqMin(-1.0),
      siteFreqMax(-1.0),
      siteMACMin(-1),

      annoRegex(NULL),
      onlyVariantSite(false),

      // individual filter
      indvDepthMin(-1),
      indvDepthMax(-1),
      indvQualMin(-1) {
  };

  // setter function
  void setUseSiteDepthFromInfo(){
    this->siteDepthFromInfo = true;
  };
  void setSiteDepthMin(int d) {
    this->siteDepthMin = d;
  }
  void setSiteDepthMax(int d) {
    this->siteDepthMax = d;
  }
  void setSiteQualMin(int q) {
    this->siteQualMin = q;
  }
  void setUseFreqFromInfo() {
    this->freqFromInfo = true;
  };
  void setSiteFreqMin(double f) {
    this->siteFreqMin = f;
  };
  void setSiteFreqMax(double f) {
    this->siteFreqMax = f;
  };
  void setSiteMACMin(int n) { // minar allele count
    this->siteMACMin = n;
  };
  int setAnnoType(const char* s) {
    return this->annoRegex.readPattern(s);
  };
  void setVariantSiteOnly() {
    this->onlyVariantSite = true;
  };
  void setIndividualDepthMin(int d) {
    this->indvDepthMin = d;
  };
  void setIndividualDepthMax(int d) {
    this->indvDepthMax = d;
  };
  void setIndividualQualMin(int q) {
    this->indvQualMin = q;
  };

  // checking functions
  bool useSiteDepthFromInfo() const {
    return this->siteDepthFromInfo;
  }
  bool siteDepthOK(int d) const {
    if (this->siteDepthMin > 0  && this->siteDepthMin > d) return false;
    if (this->siteDepthMax > 0  && this->siteDepthMax < d) return false;
    return true;
  };
  bool siteQualOK(double q) const {
    if ( this->siteQualMin > 0 && this->siteQualMin > q) return false;
    return true;
  };
  bool useSiteFreqFromInfo() const {
    return this->siteFreqFromInfo;
  };
  bool siteFreqOK(double f) const {
    if (this->siteFreqMin > 0  && this->siteFreqMin > f) return false;
    if (this->siteFreqMax > 0  && this->siteFreqMax < f) return false;
    return true;
  };
  bool siteMACOK(int n) const {
    if ( this->siteMACMin > 0 && this->siteMACMin > n) return false;
    return true;
  };
  
  bool requiredAnnotation() const {
    return this->annoRegex.isInitialized();
  }
  bool matchAnnotatoin(const char* s) const {
    return this->annoRegex.match(s);
  }
  bool isVariantSiteOnly() const {
    return this->onlyVariantSite;
  };
  bool individualDepthOK(int d) const {
    if (this->individualDepthMin > 0  && this->individualDepthMin > d) return false;
    if (this->individualDepthMax > 0  && this->individualDepthMax < d) return false;
    return true;
  };
  bool individualQualOK(double q) const {
    if ( this->individualQualMin > 0 && this->individualQualMin <= q) return false;
    return true;
  };
  
private:
  // thresholds
  bool siteDepthFromInfo;
  int  siteDepthMin;
  int  siteDepthMax;
  int  siteQualMin;

  bool   siteFreqFromInfo;
  double siteFreqMin;
  double siteFreqMax;
  int    siteMACMin;

  Regex  annoRegex;          // for filter ANNO
  bool   onlyVariantSite;     // only extract sites that are polymorphism

  int    indvDepthMin;
  int    indvDepthMax;
  int    indvQualMin;
};

#endif /* _VCFFILTER_H_ */
