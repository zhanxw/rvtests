#include "VCFInputFile.h"
#include "VCFFilter.h"

#include <string>

class VCFExtractor: public VCFInputFile, public VCFSiteFilter {
public:
VCFExtractor(const char* fn): VCFInputFile(fn){
  };
  bool passFilter() {
    VCFRecord& r = this->getVCFRecord();
    VCFPeople& people = r.getPeople();
    bool missing;

    // quick checks: depth, qual, af
    if (checkSiteDepth() && useSiteDepthFromInfo()) {
      const VCFValue& v = r.getInfoTag("DP", &missing);
      if (missing)
        return false;
      if (!siteDepthOK(v.toInt())){
        return false;
      }
    };
    if (checkSiteQual()) {
      int qual = r.getQualInt();
      if (!siteQualOK(qual)) {
        return false;
      }
    }
    if (checkSiteFreq() && useSiteFreqFromInfo()) {
      const VCFValue& v = r.getInfoTag("AF", &missing);
      if (missing)
        return false;
      if (!siteFreqOK(v.toDouble())){
        return false;
      }
    }

    // check annotation
    if (requiredAnnotation()) {
      const VCFValue& v = r.getInfoTag("ANNO", &missing);
      if (missing)
        return false;
      if (!matchAnnotatoin(v.toStr())){
        return false;
      }
    }

    // shall we loop each individuals?
    if (  (checkSiteDepth() && !useSiteDepthFromInfo()) ||
          (checkSiteFreq() && !useSiteFreqFromInfo()) ||
          (checkSiteMAC()) ||
          (isVariantSiteOnly())) {

      // loop each individual and check
      /* int dp = 0; */
      /* int qual = 0; */
      int mac = 0;
      int ac = 0;
      int an = 0;
      double af = -1.0;
      double nMiss = 0; // number of missing
      int gt = -9;

      for (unsigned int i = 0; i < people.size(); i++) {
        VCFIndividual* indv = people[i];
        int GTidx = r.getFormatIndex("GT");
        if (GTidx >= 0)
          gt = indv->get(0, &missing).getGenotype();  // [0] meaning the first field of each individual
        else
          missing = true;
        if (missing) {
          nMiss ++;
        } else{
          if (gt >= 0) {
            an += 2;
            ac += gt;
          }
        }
      };
      mac = (ac + ac > an) ? an - ac : ac;
      af = an == 0 ? 0.0 : 1.0 * ac / an;

      // check if it is variant site
      if (this->isVariantSiteOnly() && ac == 0) {
        return false;
      };

      // check site depth, freq, mac
      if (!siteDepthOK(ac)) {
        return false;
      }
      if (!siteMACOK(mac)) {
        return false;
      }
      if (!siteFreqOK(af)) {
        return false;
      }
    }
    return true;
  }; // end passFilter()
};
