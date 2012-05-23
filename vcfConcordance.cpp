#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <cassert>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "Utils.h"
#include "VCFUtil.h"

#include "MathVector.h"
#include "MathMatrix.h"

#include "SiteSet.h"

class Value;

typedef std::unordered_map< std::string, Value > ConcordanceType;
typedef std::map<std::string, ConcordanceType> AllConcordanceType;
typedef std::vector< std::string > StringArray;

class Value{
public:
  const static int HOMREF = 1;
  const static int HET = 2;
  const static int HOMALT = 3;
  const static int MISSING = 4;
  const static int UNDEF = -1;
  const static int REFERENCE = 0;
  const static int COMPARISON = 1;
  int geno[2]; // VALUE
  int& operator[] (const int idx){
    assert (0<=idx && idx <= 1);
    return geno[idx];
  }
  void clearComparison() {
    geno[Value::COMPARISON] = -1;
  }
};

int loadGenotype(VCFInputFile& vin, AllConcordanceType* input, int idx) {
    AllConcordanceType& data = *input;
    std::string key;
  int lineNo = 0;
  while (vin.readRecord()){
    lineNo ++;
    key.clear();
    VCFRecord& r = vin.getVCFRecord();
    key += r.getChrom();
    key += ":";
    key += r.getPosStr();

    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    int GTidx = r.getFormatIndex("GT");
    if (GTidx < 0) continue;
    for (int i = 0; i < people.size(); ++i) {
        indv = people[i];
        const VCFValue& v = indv->justGet(GTidx);
        int a1 = v.getAllele1();
        int a2 = v.getAllele2();
        if (a1 == MISSING_GENOTYPE || a2 == MISSING_GENOTYPE) {
            data[indv->getName()][key][idx] = Value::MISSING;
        } else if (a1 == 0) {
            if (a2 == 0) {
                data[indv->getName()][key][idx] = Value::HOMREF;
            } else if (a2 == 1) {
                data[indv->getName()][key][idx] = Value::HET;
            }
        } else if (a1 == 1) {
            if (a2 == 0) {
                data[indv->getName()][key][idx] = Value::HET;
            } else if (a2 == 1) {
                data[indv->getName()][key][idx] = Value::HOMALT;
            }
        }
    }
  };
  fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
  return lineNo;
};

void clearGenotype(AllConcordanceType* input) {
  AllConcordanceType& data = *input;
  for (AllConcordanceType::iterator iter = data.begin();
       iter != data.end();
       iter ++) {
    ConcordanceType& v = iter->second;
    for (ConcordanceType::iterator i = v.begin(); 
         i != v.end();
         i++) {
      i->second.clearComparison();
    }
  }
};

void printComparision(AllConcordanceType& data) {
  for (AllConcordanceType::iterator iter = data.begin();
       iter != data.end();
       iter ++) {
    // count type (fill in 4 by 4 matrix) , we intentionally left row 0 and column 0 unused.
    int c[5][5] = {0}; // concordance matrix

    // calculate non-ref concordance
    // calculate discovery rate
    int discoveredVariant = 0;
    int undiscoveredVariant = 0;
    printf("Comparison for people %s\n", iter->first.c_str());
    for (ConcordanceType::iterator i = iter->second.begin();
         i != iter->second.end();
         i ++) {
      c[ (i->second)[Value::REFERENCE] ] [ (i->second)[Value::COMPARISON] ] ++;
      if ( (i->second)[Value::REFERENCE] != Value::HOMREF) {
        if ( (i->second)[Value::COMPARISON] == Value::HOMREF) {
          undiscoveredVariant ++;
        } else if ( (i->second)[Value::COMPARISON] == Value::HOMREF) {
          discoveredVariant ++;
        }
      }
    }

    // outputs
    int nonRefConcordNum = c[Value::HET][Value::HET] + 
                           c[Value::HOMALT][Value::HOMALT];
    int nonRefConcordDom = nonRefConcordNum +
                           c[Value::HOMREF][Value::HET] + 
                           c[Value::HOMREF][Value::HOMALT] + 
                           c[Value::HET][Value::HOMREF] + 
                           c[Value::HET][Value::HOMALT] +
                           c[Value::HOMALT][Value::HET] + 
                           c[Value::HOMALT][Value::HOMREF];
    printf("Nonref-Concordance= %10.2f \t DiscoveryRate = %10.2f \n", 
           1.0 * nonRefConcordNum / nonRefConcordDom,
           1.0 * discoveredVariant / (discoveredVariant + undiscoveredVariant));
  }
  return;
};

size_t set_union(const StringArray& input, 
               StringArray* output) {
  output->resize(input.size());
  StringArray::iterator i;
  i = set_union(input.begin(), input.end(), input.end(), input.end(), output->begin());
  output->resize( i - output->begin());
  return output->size();
};

size_t set_intersection(const StringArray& input1, 
                        const StringArray& input2, 
                        StringArray* output) {
  output->resize(std::max(input1.size(), input2.size()));
  StringArray::iterator i;
  i = set_union(input1.begin(), input1.end(), input2.begin(), input2.end(), output->begin());
  output->resize( i - output->begin());
  return output->size();
};


int main(int argc, char** argv){
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST(pl)
      ADD_PARAMETER_GROUP(pl, "Input/Output")
      ADD_STRING_PARAMETER(pl, inVcf, "--inVcf", "input VCF File")
      ADD_PARAMETER_GROUP(pl, "Site Filter")
      ADD_STRING_PARAMETER(pl, rangeList, "--rangeList", "Specify some ranges to use, please use chr:begin-end format.")
      ADD_STRING_PARAMETER(pl, rangeFile, "--rangeFile", "Specify the file containing ranges, please use chr:begin-end format.")
      END_PARAMETER_LIST(pl)
      ;

  pl.Read(argc, argv);
  pl.Status();

  REQUIRE_STRING_PARAMETER(FLAG_inVcf, "Please provide input file using: --inVcf");

  // load referene
  const char* fn = FLAG_inVcf.c_str();
  VCFInputFile vin(fn);
  StringArray refPeople;
  vin.getVCFHeader()->getPeopleName(&refPeople);
  
  VCFInputFile** compareVcfs = new VCFInputFile* [FLAG_REMAIN_ARG.size()];
  StringArray comparePeopleNames;
  for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
      compareVcfs[i] = new VCFInputFile(FLAG_REMAIN_ARG[i]);
      StringArray names;
      compareVcfs[i]->getVCFHeader()->getPeopleName(&names);
      for (unsigned int j = 0; j < names.size(); j++) {
          comparePeopleNames.push_back(names[j]);
      };
  };
  StringArray unionPeopleNames;
  StringArray commonNames;
  set_union(comparePeopleNames, &unionPeopleNames);
  set_intersection(refPeople, unionPeopleNames, &commonNames);
 
  vin.excludeAllPeople();
  vin.includePeople(commonNames);
  
  // set range filters here
  // e.g.
  // vin.setRangeList("1:69500-69600");
  vin.setRangeList(FLAG_rangeList.c_str());
  vin.setRangeFile(FLAG_rangeFile.c_str());
  AllConcordanceType data;
  loadGenotype(vin, &data, Value::REFERENCE);

  for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
    fprintf(stderr, "Process %s ... \n", FLAG_REMAIN_ARG[i].c_str());
    loadGenotype(*compareVcfs[i], &data, Value::COMPARISON);
    printComparision(data);
    clearGenotype(&data);
  }

  VCFInputFile** compareVcfs = new VCFInputFile* [FLAG_REMAIN_ARG.size()];
  for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
    delete compareVcfs[i];
  }
  delete[] compareVcfs;

  currentTime = time(0);
  fprintf(stderr, "Analysis end at: %s", ctime(&currentTime));
  return 0;
};
