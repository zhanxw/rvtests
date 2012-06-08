#include "VCFInputFile.h"
#include "VCFFilter.h"

#include <string>

class VCFExtractor: public VCFInputFile, public VCFFilter {
public:
VCFExtractor(const char* fn): VCFInputFile(fn){
  };
  bool passFilter() const{
    VCFRecord& record = vin->getVCFRecord();
    VCFPeople& people = record.getPeople();
    bool missing;

    // quick checks
    if (siteDepthOK(-1) && useSiteDepthFromInfo()) {
      VCFValue& v = r.getInfoTag("DP", &missing);
      if (missing) 
        return false;
      if (!siteDepthOK(v.toInt())){
        return false;
      }
    };
    if (siteFreqOK(-1.0) && useSiteFreqFromInfo()) {
       VCFValue& v = r.getInfoTag("AF", &missing);
      if (missing) 
        return false;
      if (!siteFreqOK(v.toDouble())){
        return false;
      }
    }     

    // check annotation
    if (requiredAnnotation()) {
      VCFValue& v  r.getInfoTag("ANNO", &missing);
      if (missing) 
        return false;
      if (!matchAnnotatoin(v.toStr())){
        return false;
      }
    }     

    // loop each individual and check
    int dp = 0;
    int qual = 0;
    int mac = 0;
    int ac = 0;
    int an = 0;
    double af;
    double nMiss; // number of missing

    for (int i = 0; i < people.sizse(); i++) {
      VCFIndividual* indv = people[i];
      int GTidx = r.getFormatIndex("GT");
      if (GTidx >= 0)
        int gt = indv->get(0, &missing).toGenotype();  // [0] meaning the first field of each individual
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
    mac = (ac + ac > an) ? an - ac; ac;
    af = 1.0 * (ac + 1) / (an + 2);
    
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
    return true;
  };
};
