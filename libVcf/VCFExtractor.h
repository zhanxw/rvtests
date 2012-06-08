#include "VCFUtil.h"
#include "Regex.h"
#include <string>

class VCFExtractor{
public:
VCFExtractor(const char* fn):
  indvDepthMin(-1),
      indvDepthMax(-1),
      indvQualMin(-1),
      depthFromInfo(true),      // read depth from INFO field
      siteDepthMin(-1),
      siteDepthMax(-1),
      siteQualMin(-1),
      freqFromInfo(true),      // read AF from INFO field
      siteFreqMin(-1.0),
      siteFreqMax(-1.0),
      siteMACMin(-1),
      annoRegex(NULL),
      onlyVariantSite(false) {
    super(fn);
    this->annoRegex = new Regex;
  };
  virtual ~VCFExtractor() {
    if (this->annoRegex) {
      delete this->annoRegex;
      this->annoRegex = NULL;
    }
  };

  void setIndvDepthMin(int d) {
    this->indvDepthMin = d;
  };
  void setIndvDepthMax(int d) {
    this->indvDepthMax = d;
  };
  void setIndvQualityMin(int d) {
    this->indvQualMin = d;
  };
  void setSiteDepthMin(int d) {
    this->siteDepthMin = d;
  };
  void setSiteDepthMax(int d) {
    this->siteDepthMax = d;
  };
  void setSiteQualityMin(int d) {
    this->siteQualMin = d;
  };
  void setSiteFreqMax(double d) {
    this->siteFreqMax = d;
  };
  void setSiteFreqMin(double d) {
    this->siteFreqMin = d;
  };
  void setSiteMACMin(int d) {
    this->siteMACMin = d;
  };
  void setOnlyVariantSite(bool b) {
    this->onlyVariantSite = b;
  };
  void setAnnotationFilter(const char* regex) {
    this->annoRegex->readPattern(regex);
  };

  int extract() {
    if (!this->vin) { return -1; }

    // read people id
    VCFHeader* h = this->vin->getVCFHeader();
    h->getPeopleName(&this->peopleId); // include all people
    int peopleSize = this->peopleId.size();

    // read each genotype
    int numRecordRead = 0;
    std::string snvId;
    bool checkGD = (this->indvDepthMax > -1) || (this->indvDepthMin > -1);
    bool checkGQ = (this->indvQualMin > -1);
    std::vector<int> genoVec;
    int qual = 0;
    int dp = 0;
    double freq = 0.0;
    while (vin->readRecord()){
      if (++numRecordRead % 10000 == 0) {
        fprintf(stdout, "Total %d record processed.\n", numRecordRead);
      }

      VCFRecord& record = vin->getVCFRecord();
      VCFPeople& people = record.getPeople();

      // apply site filter: QUAL
      if (this->siteQualMin >= 0) {
        qual = atoi(record.getQual());
        if (qual < siteQualMin)
          continue;
      }

      // apply site filter: ANNO
      if (this->annoRegex) {
        const char* v = record.getInfoTag("ANNO");
        if (! this->annoRegex->match(v) )
          continue;
      }
      // apply site filter: depth
      if (this->depthFromInfo) {
        dp = atoi(record.getInfoTag("DP"));
        if (this->siteDepthMax > -1) {
          if (dp > this->siteDepthMax)
            continue;
        }
        if (this->siteDepthMin > -1) {
          if (dp < this->siteDepthMin)
            continue;
        }
      }

      // apply site filter: frequency
      if (this->freqFromInfo) {
        freq = atof(record.getInfoTag("AF"));
        if (this->siteFreqMax > -1) {
          if (freq > this->siteFreqMax)
            continue;
        }
        if (this->siteFreqMin > -1) {
          if (freq < this->siteFreqMin)
            continue;
        }
      }

      // Loop into individual
      bool isVariantSite = false;
      int totalDepth = 0;
      int totalAN = 0;
      int totalAC = 0;
      int totalMissing = 0;
      genoVec.resize(0);
      VCFIndividual* indv;
      for (int i = 0; i < peopleSize; i++) {
        totalAN += 2;
        indv = people[i];

        int geno;
        // get GT index. if you are sure the index will not change, call this function only once!
        int GTidx = record.getFormatIndex("GT");
        if (GTidx < 0) {
          geno = MISSING_GENOTYPE;
          totalMissing += 2;
          continue;
        }
        geno = (*indv)[GTidx].getGenotype();

        // get GD index
        if (checkGD) {
          int GDidx = record.getFormatIndex("GD");
          int gd = 0;
          (*indv)[GDidx].toInt(&gd);
          if (!this->depthFromInfo) totalDepth += gd;
          if (this->indvDepthMin > 0 && this->indvDepthMin > gd) {
            geno = MISSING_GENOTYPE;
            totalMissing += 2;
            continue;
          }
          if (this->indvDepthMax > 0 && this->indvDepthMax < gd) {
            totalMissing += 2;
            geno = MISSING_GENOTYPE;
            continue;
          }
        }

        // get GQ index
        if (checkGQ) {
          int GQidx = record.getFormatIndex("GQ");
          int gq = 0;
          (*indv)[GQidx].toInt(&gq);
          if (this->indvQualMin > 0 && this->indvQualMin > gq) {
            totalMissing += 2;
            geno = MISSING_GENOTYPE;
            continue;
          }
        };
        // set variant tag
        if (geno != 0) { // 0 meaning 0/0
          isVariantSite = true;
        }
        // store genoType
        genoVec.push_back(geno);
      } // end each individual

      // apply site filter: Depth
      if (!this->depthFromInfo) {
        if (this->siteDepthMax > -1) {
          if (totalDepth > this->siteDepthMax)
            continue;
        }
        if (this->siteDepthMin > -1) {
          if (totalDepth < this->siteDepthMin)
            continue;
        }
      }

      // apply site filter: frequency
      if (!this->freqFromInfo) {
        if (totalAN <= 0){  // no allele here!
          fprintf(stderr, "totalAN == 0!\n");
          continue;
        }

        freq = (double)(totalAC) / totalAN;
        if (this->siteFreqMax > -1) {
          if (freq > this->siteFreqMax)
            continue;
        }
        if (this->siteFreqMin > -1) {
          if (freq < this->siteFreqMin)
            continue;
        }
      }

      // apply site filter: MAC
      if (this->siteMACMin > -1) {
        if (totalAC < this->siteMACMin)
          continue;
      }

      // apply site filter: output only variant site?
      if (this->onlyVariantSite && !isVariantSite) {
        continue;
      }

      // store site info
      snvId = record.getID();
      if (snvId == ".") {
        snvId = record.getChrom();
        snvId += ':';
        snvId += toString(record.getPos());
      }
      this->chrom.push_back(record.getChrom());
      this->pos.push_back(record.getPos());
      this->ref.push_back(record.getRef());
      this->alt.push_back(record.getAlt());
      this->id.push_back(snvId);

      this->ac.push_back(totalAC);
      this->an.push_back(totalAN);
      this->nMissing.push_back(totalMissing);
      this->af.push_back(freq);

      // store genotype
      int nr = this->geno.rows;
      int nc = this->geno.cols;
      this->geno.Dimension(nr + 1, nc);
      for (int i = 0; i < nc; i++) {
        this->geno[nr][i] = genoVec[i];
      }

    };  // for each VCFRecord
  };
#if 0
  int open(const char* fn) {
    this->vin = new VCFInputFile(fn);
    if (!this->vin){
      return 1;
    };
  };
  void close() {
    if (this->vin)
      delete this->vin;
    this->vin = NULL;
  };

  // code to extract data

  int extractGenoType(SimpleMatrix* m) {
    assert(m);
    int nr = this->geno.rows;
    int nc = this->geno.cols;
    m->Dimension(nr,nc);
    for (int r = 0; r < nr; r++){
      for (int c = 0; c < nc; c++ ){
        (*m)[nr][nc] = geno[nr][nc];
      }
    }
    return 0;
  };
#endif
#if 0
  const std::vector<int>& getAC() {return this->ac;};
  const std::vector<int>& getAN() {return this-> an;};
  const std::vector<int>& getMissing() {return this-> nMissing;}; // # allele missing
  const std::vector<double>& getAF() {return this-> af; };

  // get extracted variants chrom, pos, ref, alt, names
  const std::vector<std::string>& getChrom() {return this-> chrom;};
  const std::vector<int>& getPos() {return this-> pos;};
  const std::vector<std::string>& getRef() { return this-> ref;};
  const std::vector<std::string>& getAlt() { return this-> alt;};
  const std::vector<std::string>& getId() {return this-> id; };
#endif
private:
  // thresholds
  int indvDepthMin;
  int indvDepthMax;
  int indvQualMin;
  int siteDepthMin;
  int siteDepthMax;
  int siteQualMin;
  bool depthFromInfo;
  bool freqFromInfo;
  double siteFreqMin;
  double siteFreqMax;
  int siteMACMin;
  Regex* annoRegex;          // for filter ANNO
  bool onlyVariantSite;     // only extract sites that are polymorphism

  // store site stats.
  std::vector<int> ac;       // 1 allele count
  std::vector<int> an;       // total allele count
  std::vector<int> nMissing; // # allele missing
  std::vector<double> af;

  // get extracted people ids
  std::vector<std::string> peopleId;

  // get extracted variants chrom, pos, ref, alt, names
  std::vector<std::string> chrom;
  std::vector<int> pos;
  std::vector<std::string> ref;
  std::vector<std::string> alt;
  std::vector<std::string> id;

  // get extract genotypes
  SimpleMatrix geno;
};
