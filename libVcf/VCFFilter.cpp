#include "VCFFilter.h"
#include "base/ParRegion.h"

VCFSiteFilter::VCFSiteFilter()
    :                            // Site filter
      siteDepthFromInfo(false),  // read depth from INFO field
      siteDepthMin(-1),
      siteDepthMax(-1),
      siteQualMin(-1),

      siteFreqFromInfo(false),  // read AF from INFO field
      siteFreqMin(-1.0),        // here freq means minor allele frequency
      siteFreqMax(-1.0),
      siteMACMin(-1),

      onlyVariantSite(false),
      parRegion(NULL),
      chromXExtraction(ANY) {
#if 0
  // individual filter
  indvDepthMin(-1),
      indvDepthMax(-1),
      indvQualMin(-1)
#endif
};
VCFSiteFilter::~VCFSiteFilter() {
  if (parRegion) {
    delete parRegion;
    parRegion = NULL;
  }
}
