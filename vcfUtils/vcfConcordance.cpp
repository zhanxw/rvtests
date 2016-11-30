#include "Argument.h"
#include "IO.h"
#include "tabix.h"

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Utils.h"
#include "VCFUtil.h"

#include "MathMatrix.h"
#include "MathVector.h"

#include "SiteSet.h"

class Value;

typedef std::unordered_map<std::string, Value> ConcordanceType;
typedef std::map<std::string, ConcordanceType> AllConcordanceType;
typedef std::vector<std::string> StringArray;

class Value {
 public:
  const static int HOMREF = 1;
  const static int HET = 2;
  const static int HOMALT = 3;
  const static int MISSING = 4;
  const static int UNDEF = 0;
  const static int REFERENCE = 0;
  const static int COMPARISON = 1;
  int& operator[](const int idx) {
    assert(0 <= idx && idx <= 1);
    return geno[idx];
  }
  void clearComparison() { geno[Value::COMPARISON] = Value::UNDEF; }
  Value() {
    geno[0] = UNDEF;
    geno[1] = UNDEF;
  };

 public:
  int geno[2];  // VALUE
};

void dump(const StringArray& a) {
  for (size_t i = 0; i < a.size(); i++) {
    fprintf(stderr, "[ %zu ] = \"%s\", ", i, a[i].c_str());
  }
  fprintf(stderr, "\n");
}

int loadGenotype(VCFInputFile& vin, AllConcordanceType* input, int idx) {
  AllConcordanceType& data = *input;
  std::string key;
  int lineNo = 0;
  while (vin.readRecord()) {
    lineNo++;
    key.clear();
    VCFRecord& r = vin.getVCFRecord();
    key += r.getChrom();
    key += ":";
    key += r.getPosStr();

    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    int GTidx = r.getFormatIndex("GT");
    if (GTidx < 0) continue;
    for (size_t i = 0; i < people.size(); ++i) {
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
  fprintf(stderr, "Total %d VCF records have read successfully\n", lineNo);
  return lineNo;
};

void clearGenotype(AllConcordanceType* input) {
  AllConcordanceType& data = *input;
  for (AllConcordanceType::iterator iter = data.begin(); iter != data.end();
       iter++) {
    ConcordanceType& v = iter->second;
    for (ConcordanceType::iterator i = v.begin(); i != v.end(); i++) {
      i->second.clearComparison();
    }
  }
};

void printComparision(const VCFInputFile& vin, AllConcordanceType& data,
                      const StringArray& names) {
  std::set<std::string> nameSet;
  for (unsigned int i = 0; i < names.size(); i++) {
    nameSet.insert(names[i]);
  }

  for (AllConcordanceType::iterator iter = data.begin(); iter != data.end();
       iter++) {
    if (nameSet.find(iter->first) == nameSet.end()) {
      // fprintf(stderr, "Skip people [ %s ] as it is not in reference VCF\n",
      // iter->first.c_str());
      continue;
    }

    // count type (fill in 4 by 4 matrix) , we intentionally left row 0 and
    // column 0 unused.
    int c[5][5] = {{0}, {0}, {0}, {0}, {0}};  // concordance matrix

    // calculate non-ref concordance
    // calculate discovery rate
    // calculate standard/input sites, and overlapping sites...
    int discoveredVariant = 0;
    int undiscoveredVariant = 0;
    for (ConcordanceType::iterator i = iter->second.begin();
         i != iter->second.end(); i++) {
      // concordance table
      c[(i->second)[Value::REFERENCE]][(i->second)[Value::COMPARISON]]++;
      // discovery rate
      if ((i->second)[Value::REFERENCE] == Value::HET ||
          (i->second)[Value::REFERENCE] ==
              Value::HOMALT) {  // ref is variant site
        if ((i->second)[Value::COMPARISON] == Value::HOMREF) {
          undiscoveredVariant++;
        } else if ((i->second)[Value::COMPARISON] == Value::HET ||
                   (i->second)[Value::COMPARISON] == Value::HOMALT) {
          discoveredVariant++;
        }
      }
    }

    // outputs
    // printf("Comparison for people %s\n", iter->first.c_str());
    int nonRefConcordNum =
        c[Value::HET][Value::HET] + c[Value::HOMALT][Value::HOMALT];
    int nonRefConcordDom =
        nonRefConcordNum + c[Value::HOMREF][Value::HET] +
        c[Value::HOMREF][Value::HOMALT] + c[Value::HET][Value::HOMREF] +
        c[Value::HET][Value::HOMALT] + c[Value::HOMALT][Value::HET] +
        c[Value::HOMALT][Value::HOMREF];
    // printf("Nonref-Concordance= %10f \t DiscoveryRate = %10f \n",
    //        1.0 * nonRefConcordNum / nonRefConcordDom,
    //        1.0 * discoveredVariant / (discoveredVariant +
    //        undiscoveredVariant));
    // printf("%d\t%d\n", discoveredVariant, undiscoveredVariant);

    // // print the 5 by 5 table
    // const char* title[] = {"Undef", "HomeRef", "Het", "Homalt", "Missing"};
    // for (int i = 0; i < 5; ++i) {
    //   printf("\t%s", title[i]);
    // }
    // printf("\n");
    // for (int i = 0; i <= 4; i++ ) {
    //   printf("%s", title[i]);
    //     for (int j = 0; j <= 4; j++ ){
    //         printf("\t%d" , c[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // following counts does not take 'missing' (bottom row and right most
    // column)
    // overlap: 3 (homref, het, homalt) by 3 matrix
    int overlap =
        c[Value::HOMREF][Value::HOMREF] + c[Value::HOMREF][Value::HET] +
        c[Value::HOMREF][Value::HOMALT] + c[Value::HET][Value::HOMREF] +
        c[Value::HET][Value::HET] + c[Value::HET][Value::HOMALT] +
        c[Value::HOMALT][Value::HOMREF] + c[Value::HOMALT][Value::HET] +
        c[Value::HOMALT][Value::HOMALT];
    // stdTotal: 3 (homref, het, homalt) by 4 matrix
    int stdTotal =
        c[Value::HOMREF][Value::UNDEF] + c[Value::HOMREF][Value::HOMREF] +
        c[Value::HOMREF][Value::HET] + c[Value::HOMREF][Value::HOMALT] +
        c[Value::HET][Value::UNDEF] + c[Value::HET][Value::HOMREF] +
        c[Value::HET][Value::HET] + c[Value::HET][Value::HOMALT] +
        c[Value::HOMALT][Value::UNDEF] + c[Value::HOMALT][Value::HOMREF] +
        c[Value::HOMALT][Value::HET] + c[Value::HOMALT][Value::HOMALT];
    // inputTotal: 4 (undef, hom, het, homalt) by 3 matrix
    int inputTotal =
        c[Value::UNDEF][Value::HOMREF] + c[Value::UNDEF][Value::HET] +
        c[Value::UNDEF][Value::HOMALT] + c[Value::HOMREF][Value::HOMREF] +
        c[Value::HOMREF][Value::HET] + c[Value::HOMREF][Value::HOMALT] +
        c[Value::HET][Value::HOMREF] + c[Value::HET][Value::HET] +
        c[Value::HET][Value::HOMALT] + c[Value::HOMALT][Value::HOMREF] +
        c[Value::HOMALT][Value::HET] + c[Value::HOMALT][Value::HOMALT];
    //     // 8 cells in overlap, removing bottom top one
    //     int nonRefOverlap = overlap - c[Value::HOMREF][Value::HOMREF];
    // stdVariant: 2 (het, homalt) by 4 matrix
    int stdVariant =
        c[Value::HET][Value::UNDEF] + c[Value::HET][Value::HOMREF] +
        c[Value::HET][Value::HET] + c[Value::HET][Value::HOMALT] +
        c[Value::HOMALT][Value::UNDEF] + c[Value::HOMALT][Value::HOMREF] +
        c[Value::HOMALT][Value::HET] + c[Value::HOMALT][Value::HOMALT];
    // inputTotal: 4 (undef, hom, het, homalt) by 2 matrix
    int inputVariant =
        c[Value::UNDEF][Value::HET] + c[Value::UNDEF][Value::HOMALT] +
        c[Value::HOMREF][Value::HET] + c[Value::HOMREF][Value::HOMALT] +
        c[Value::HET][Value::HET] + c[Value::HET][Value::HOMALT] +
        c[Value::HOMALT][Value::HET] + c[Value::HOMALT][Value::HOMALT];
    // stdVariantInInput: 2 (het, homalt) by 3 matrix (homref, het, homalt)
    int stdVariantInInput =
        c[Value::HET][Value::HOMREF] + c[Value::HET][Value::HET] +
        c[Value::HET][Value::HOMALT] + c[Value::HOMALT][Value::HOMREF] +
        c[Value::HOMALT][Value::HET] + c[Value::HOMALT][Value::HOMALT];
    // stdOnly: 3 by 1 matrix
    int stdOnly = c[Value::HOMREF][Value::UNDEF] + c[Value::HET][Value::UNDEF] +
                  c[Value::HOMALT][Value::UNDEF];
    // inputOnly: 1 by 3 matrix
    int inputOnly = c[Value::UNDEF][Value::HOMREF] +
                    c[Value::UNDEF][Value::HET] +
                    c[Value::UNDEF][Value::HOMALT];

    // printf("File\t"
    //        "PeopleId\t"
    //        "Overlap\t"
    //        "Std_total\t"
    //        "Input_total\t"
    //        "Std_only\t"
    //        "Input_only\t"
    //        "nonRefConcord_overlap\t"
    //        "Std_variants\t"
    //        "Input_variants\t"
    //        "Std_variants_in_Input\t"
    //        "Ptg_Std_variants_in_Input\t"
    //        "HomR/HomR\tHomR/Het\tHomR/HomA\t"
    //        "Het/HomR\tHet/Het\tHet/HomA\t"
    //        "HomA/HomR\tHomA/Het\tHomA/HomA\n");

    printf("%s\t", vin.getFileName());
    printf("%s\t", iter->first.c_str());
    printf("%d\t", overlap);
    printf("%d\t", stdTotal);
    printf("%d\t", inputTotal);
    printf("%d\t", stdOnly);
    printf("%d\t", inputOnly);

    printf("%.6f\t", 1.0 * nonRefConcordNum / nonRefConcordDom);
    printf("%d\t", stdVariant);
    printf("%d\t", inputVariant);
    printf("%d\t", stdVariantInInput);
    printf("%.6f\t", 1.0 * stdVariantInInput / stdVariant);
    printf("%d\t", c[Value::HOMREF][Value::HOMREF]);
    printf("%d\t", c[Value::HOMREF][Value::HET]);
    printf("%d\t", c[Value::HOMREF][Value::HOMALT]);
    printf("%d\t", c[Value::HET][Value::HOMREF]);
    printf("%d\t", c[Value::HET][Value::HET]);
    printf("%d\t", c[Value::HET][Value::HOMALT]);
    printf("%d\t", c[Value::HOMALT][Value::HOMREF]);
    printf("%d\t", c[Value::HOMALT][Value::HET]);
    printf("%d\n", c[Value::HOMALT][Value::HOMALT]);
  }
  return;
};

size_t set_union(const StringArray& input, StringArray* output) {
  // sort input
  StringArray tmp = input;
  sort(tmp.begin(), tmp.end());

  output->resize(tmp.size());
  StringArray::iterator i;
  i = set_union(tmp.begin(), tmp.end(), tmp.end(), tmp.end(), output->begin());
  output->resize(i - output->begin());
  return output->size();
};

size_t set_intersection(const StringArray& input1, const StringArray& input2,
                        StringArray* output) {
  StringArray t1 = input1;
  StringArray t2 = input2;
  sort(t1.begin(), t1.end());
  sort(t2.begin(), t2.end());

  output->resize(std::max(t1.size(), t2.size()));
  StringArray::iterator i;
  i = set_intersection(t1.begin(), t1.end(), t2.begin(), t2.end(),
                       output->begin());
  output->resize(i - output->begin());
  return output->size();
};

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(s, "-s", "standard VCF file (to which other files compare")
ADD_PARAMETER_GROUP("Site Filter")
ADD_STRING_PARAMETER(
    rangeList, "--rangeList",
    "Specify some ranges to use, please use chr:begin-end format.")
ADD_STRING_PARAMETER(
    rangeFile, "--rangeFile",
    "Specify the file containing ranges, please use chr:begin-end format.")
ADD_STRING_PARAMETER(
    siteFile, "--siteFile",
    "Specify the file containing chromosomal sites, please use chr:pos")
END_PARAMETER_LIST();

int main(int argc, char** argv) {
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  PARSE_PARAMETER(argc, argv);
  PARAMETER_STATUS();

  REQUIRE_STRING_PARAMETER(FLAG_s, "Please provide input file using: -s");

  // load referene
  const char* fn = FLAG_s.c_str();
  VCFInputFile vin(fn);
  StringArray refPeople;
  vin.getVCFHeader()->getPeopleName(&refPeople);

  VCFInputFile** compareVcfs = new VCFInputFile*[FLAG_REMAIN_ARG.size()];
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
  fprintf(stderr, "Total %d samples are included.\n", (int)commonNames.size());

  // set range filters here
  // e.g.
  // vin.setRangeList("1:69500-69600");
  vin.setRangeList(FLAG_rangeList.c_str());
  vin.setRangeFile(FLAG_rangeFile.c_str());
  vin.setSiteFile(FLAG_siteFile.c_str());
  AllConcordanceType data;
  loadGenotype(vin, &data, Value::REFERENCE);

  printf(
      "File\t"
      "PeopleId\t"
      "Overlap\t"
      "Std_total\t"
      "Input_total\t"
      "Std_only\t"
      "Input_only\t"
      "nonRefConcord_overlap\t"
      "Std_variants\t"
      "Input_variants\t"
      "Std_variants_in_Input\t"
      "Ptg_Std_variants_in_Input\t"
      "HomR/HomR\tHomR/Het\tHomR/HomA\t"
      "Het/HomR\tHet/Het\tHet/HomA\t"
      "HomA/HomR\tHomA/Het\tHomA/HomA\n");

  for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
    fprintf(stderr, "Process %s ... \n", FLAG_REMAIN_ARG[i].c_str());
    compareVcfs[i]->setRangeList(FLAG_rangeList.c_str());
    compareVcfs[i]->setRangeFile(FLAG_rangeList.c_str());
    compareVcfs[i]->setSiteFile(FLAG_siteFile.c_str());
    loadGenotype(*compareVcfs[i], &data, Value::COMPARISON);

    StringArray names;
    compareVcfs[i]->getVCFHeader()->getPeopleName(&names);
    printComparision(*compareVcfs[i], data, names);

    clearGenotype(&data);
  }

  for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
    delete compareVcfs[i];
  }
  delete[] compareVcfs;

  currentTime = time(0);
  fprintf(stderr, "Analysis end at: %s", ctime(&currentTime));
  return 0;
};
