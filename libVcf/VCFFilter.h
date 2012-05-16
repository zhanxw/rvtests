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
  depthFromInfo(true),      // read depth from INFO field
      siteDepthMin(-1),
      siteDepthMax(-1),
      siteQualMin(-1),

      freqFromInfo(true),      // read AF from INFO field
      siteFreqMin(-1.0),
      siteFreqMax(-1.0),
      siteMACMin(-1),

      // annoRegex(NULL),
      onlyVariantSite(false),

      // individual filter
      indvDepthMin(-1),
      indvDepthMax(-1),
      indvQualMin(-1)
      {
      };

private:
  // thresholds
  bool depthFromInfo;
  int siteDepthMin;
  int siteDepthMax;
  int siteQualMin;

  bool freqFromInfo;
  double siteFreqMin;
  double siteFreqMax;
  int siteMACMin;

  Regex annoRegex;          // for filter ANNO
  bool onlyVariantSite;     // only extract sites that are polymorphism

  int indvDepthMin;
  int indvDepthMax;
  int indvQualMin;
};


#endif /* _VCFFILTER_H_ */
