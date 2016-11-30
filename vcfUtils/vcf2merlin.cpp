/**
   immediately TODO:
   8. force loading index when read by region.
   3. support tri-allelic (getAlt())
   4. speed up VCF parsing. (make a separate line buffer).
   5. loading phenotype and covariate (need tests now).
   6. do analysis. (test CMC for now)
   7. VT (combine Collapsor and ModelFitter)

   DONE:
   1. suppport PLINK output
   2. support access INFO tag
   5. give warnings for: Argument.h detect --inVcf --outVcf empty argument value
   after --inVcf
   8. Make code easy to use ( hide PeopleSet and RangeList)
   9. Inclusion/Exclusion set should be considered sequentially.
   futher TODO:
   1. handle different format GT:GD:DP ... // use getFormatIndex()

*/
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

int main(int argc, char** argv) {
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

  ////////////////////////////////////////////////
  BEGIN_PARAMETER_LIST()
  ADD_PARAMETER_GROUP("Input/Output")
  ADD_STRING_PARAMETER(inVcf, "--inVcf", "input VCF File")
  ADD_STRING_PARAMETER(outMerlin, "--outMerlin", "output prefix")
  ADD_PARAMETER_GROUP("People Filter")
  ADD_STRING_PARAMETER(peopleIncludeID, "--peopleIncludeID",
                       "give IDs of people that will be included in study")
  ADD_STRING_PARAMETER(
      peopleIncludeFile, "--peopleIncludeFile",
      "from given file, set IDs of people that will be included in study")
  ADD_STRING_PARAMETER(peopleExcludeID, "--peopleExcludeID",
                       "give IDs of people that will be included in study")
  ADD_STRING_PARAMETER(
      peopleExcludeFile, "--peopleExcludeFile",
      "from given file, set IDs of people that will be included in study")
  ADD_PARAMETER_GROUP("Site Filter")
  ADD_STRING_PARAMETER(
      rangeList, "--rangeList",
      "Specify some ranges to use, please use chr:begin-end format.")
  ADD_STRING_PARAMETER(
      rangeFile, "--rangeFile",
      "Specify the file containing ranges, please use chr:begin-end format.")
  END_PARAMETER_LIST();

  PARSE_PARAMETER(argc, argv);
  PARAMETER_STATUS();

  if (FLAG_REMAIN_ARG.size() > 0) {
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    fprintf(stderr, "\n");
    abort();
  }

  REQUIRE_STRING_PARAMETER(FLAG_inVcf,
                           "Please provide input file using: --inVcf");

  const char* fn = FLAG_inVcf.c_str();
  VCFInputFile vin(fn);

  // set range filters here
  // e.g.
  // vin.setRangeList("1:69500-69600");
  vin.setRangeList(FLAG_rangeList.c_str());
  vin.setRangeFile(FLAG_rangeFile.c_str());

  // set people filters here
  if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
    vin.excludeAllPeople();
    vin.includePeople(FLAG_peopleIncludeID.c_str());
    vin.includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
  }
  vin.excludePeople(FLAG_peopleExcludeID.c_str());
  vin.excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());

  // let's write it out.
  FILE* fMap;  // CHROMOSOME   MARKER          POSITION
  FILE* fDat;  // A some_disease\n
               // T some_trait
               // M some_marker
               // M another_marker
  FILE* fPed;  // first 5 column: FID, IID, PID, MID, SEX; then follow Dat file
  FILE* fPid;  // Person ID file, (extra for Merlin), including all people ID as
               // they are in PED file.

  fMap = fopen((FLAG_outMerlin + ".map").c_str(), "wt");
  fDat = fopen((FLAG_outMerlin + ".dat").c_str(), "wt");
  fPed = fopen((FLAG_outMerlin + ".ped").c_str(), "wt");
  fPid = fopen((FLAG_outMerlin + ".pid").c_str(), "wt");
  assert(fMap && fDat && fPed && fPid);

  std::string marker;  // marker x people
  std::vector<std::string> allMarker;
  Matrix geno;
  fputs("CHROMOSOME\tMARKER\tPOSITION\n", fMap);

  while (vin.readRecord()) {
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    // write map file
    marker = r.getID();
    if (marker == ".") {
      fprintf(fMap, "%s\t%s:%d\t%d\n", r.getChrom(), r.getChrom(), r.getPos(),
              r.getPos());
      fprintf(fDat, "M\t%s:%d\n", r.getChrom(), r.getPos());
    } else {
      fprintf(fMap, "%s\t%s\t%d\n", r.getChrom(), marker.c_str(), r.getPos());
      fprintf(fDat, "M\t%s\n", marker.c_str());
    }
    allMarker.push_back(marker);

    geno.Dimension(allMarker.size(), people.size());

    // e.g.: get TAG from INFO field
    // fprintf(stderr, "%s\n", r.getInfoTag("ANNO"));

    int m = allMarker.size() - 1;
    // e.g.: Loop each (selected) people in the same order as in the VCF
    for (int i = 0; i < people.size(); i++) {
      indv = people[i];
      // get GT index. if you are sure the index will not change, call this
      // function only once!
      int GTidx = r.getFormatIndex("GT");
      if (GTidx >= 0) {
        geno[m][i] = (*indv)[GTidx].getGenotype();
      } else {
        fprintf(stderr, "Cannot find GT field!\n");
        abort();
      }
    }
  };
  VCFHeader* h = vin.getVCFHeader();
  std::vector<std::string> peopleId;
  h->getPeopleName(&peopleId);

  // dump PED and PID file
  for (int p = 0; p < peopleId.size(); p++) {
    fprintf(fPed, "%s\t%s\t0\t0\t0", peopleId[p].c_str(), peopleId[p].c_str());
    for (int m = 0; m < allMarker.size(); m++) {
      int g = (int)geno[m][p];
      switch (g) {
        case 0:
          fputs("\t0/0", fPed);
          break;
        case 1:
          fputs("\t0/1", fPed);
          break;
        case 2:
          fputs("\t1/1", fPed);
          break;
        default:
          fputs("x/x", fPed);
          break;
      }
    }
    fputs("\n", fPed);

    fprintf(fPid, "%s\n", peopleId[p].c_str());
  }
  return 0;
};
