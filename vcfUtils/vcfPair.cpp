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

#include "IO.h"
#include "Regex.h"

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(inVcf, "--inVcf", "input VCF File")
ADD_STRING_PARAMETER(out, "--out", "output prefix")
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
ADD_PARAMETER_GROUP("Gene Extractor")
ADD_STRING_PARAMETER(
    geneFile, "--geneFile",
    "Specify the gene file (refFlat format), so we know gene start and end.")
ADD_STRING_PARAMETER(geneName, "--gene", "Specify the gene names to extract")
ADD_STRING_PARAMETER(annoType, "--annoType",
                     "Specify the type of annotation to extract")
ADD_PARAMETER_GROUP("Genotype Filter")
ADD_INT_PARAMETER(
    minGQ, "--minGQ",
    "Specify minimum GQ required (lower than this will be marked missing).")
ADD_PARAMETER_GROUP("Other Function")
ADD_BOOL_PARAMETER(variantOnly, "--variantOnly",
                   "Only variant sites from the VCF file will be processed.")
ADD_STRING_PARAMETER(updateId, "--update-id",
                     "Update VCF sample id using "
                     "given file (column 1 and 2 are "
                     "old and new id).")
END_PARAMETER_LIST();

int main(int argc, char** argv) {
  time_t currentTime = time(0);
  fprintf(stderr, "Analysis started at: %s", ctime(&currentTime));

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
  REQUIRE_STRING_PARAMETER(FLAG_out, "Please provide output file using: --out");

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

  // // let's write it out.
  // VCFOutputFile* vout = NULL;
  // PlinkOutputFile* pout = NULL;
  // if (FLAG_outVcf.size() > 0) {
  //     vout = new VCFOutputFile(FLAG_outVcf.c_str());
  // };
  // if (FLAG_outPlink.size() > 0) {
  //     pout = new PlinkOutputFile(FLAG_outPlink.c_str());
  // };
  // if (!vout && !pout) {
  //     vout = new VCFOutputFile("temp.vcf");
  // }

  if (FLAG_updateId != "") {
    int ret = vin.updateId(FLAG_updateId.c_str());
    fprintf(stdout, "%d samples have updated id.\n", ret);
  }

  // load gene ranges
  std::map<std::string, std::string> geneRange;
  if (FLAG_geneName.size()) {
    if (FLAG_geneFile.size() == 0) {
      fprintf(stderr, "Have to provide --geneFile to extract by gene.\n");
      abort();
    }
    LineReader lr(FLAG_geneFile);
    std::vector<std::string> fd;
    while (lr.readLineBySep(&fd, "\t ")) {
      if (FLAG_geneName != fd[0]) continue;
      fd[2] = chopChr(fd[2]);  // chop "chr1" to "1"
      if (geneRange.find(fd[0]) == geneRange.end()) {
        geneRange[fd[0]] = fd[2] + ":" + fd[4] + "-" + fd[5];
      } else {
        geneRange[fd[0]] += "," + fd[2] + ":" + fd[4] + "-" + fd[5];
      }
    };
  }
  std::string range;
  for (std::map<std::string, std::string>::iterator it = geneRange.begin();
       it != geneRange.end(); it++) {
    if (range.size() > 0) {
      range += ",";
    }
    range += it->second;
  };
  if (!range.empty()) {
    // fprintf(stdout, "range = %s\n", range.c_str());
    vin.setRangeList(range.c_str());
  }
  Regex regex;
  if (FLAG_annoType.size()) {
    regex.readPattern(FLAG_annoType);
  }

  // real working park
  // if (vout) vout->writeHeader(vin.getVCFHeader());
  // if (pout) pout->writeHeader(vin.getVCFHeader());
  std::vector<std::string> name;
  vin.getVCFHeader()->getPeopleName(&name);
  int num = name.size();
  Matrix p0(num, num);  // store counts when same: 00 - 11
  Matrix p1(num, num);  // store counts when share 1 allele: 00 - 01, 01 - 11
  Matrix p2(
      num,
      num);  // store counts when share 2 allels: 00 - 00, 01 - 01 , 11 - 11
  Matrix p9(num, num);  // other
  p0.Zero();
  p1.Zero();
  p2.Zero();
  p9.Zero();
  std::vector<int> g;  // store genotype matrix
  g.reserve(num);
  std::string foutName = FLAG_out + ".pairDiff.log";
  FILE* fLog = fopen(foutName.c_str(), "wt");
  if (!fLog) {
    fprintf(stderr, "Cannot open %s for write!\n", foutName.c_str());
    return -1;
  };

  int lineNo = 0;
  int nonVariantSite = 0;
  while (vin.readRecord()) {
    lineNo++;
    VCFRecord& r = vin.getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    if (FLAG_variantOnly) {
      bool hasVariant = false;
      int geno;
      int GTidx = r.getFormatIndex("GT");
      for (size_t i = 0; i < people.size(); i++) {
        indv = people[i];
        geno = indv->justGet(GTidx).getGenotype();
        if (geno != 0 && geno != MISSING_GENOTYPE) hasVariant = true;
      }
      if (!hasVariant) {
        nonVariantSite++;
        continue;
      }
    }

    if (FLAG_annoType.size()) {
      bool isMissing = false;
      const char* tag = r.getInfoTag("ANNO", &isMissing).toStr();
      if (isMissing) continue;
      // fprintf(stdout, "ANNO = %s", tag);
      bool match = regex.match(tag);
      // fprintf(stdout, " %s \t", match ? "match": "noMatch");
      // fprintf(stdout, " %s \n", exists ? "exists": "missing");
      if (!match) {
        continue;
      }
    }

    // calculate pair wise difference
    int GTidx = r.getFormatIndex("GT");
    int GQidx = r.getFormatIndex("GQ");
    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];
      int gq = indv->justGet(GQidx).toInt();
      if (FLAG_minGQ > 0) {
        if (gq >= FLAG_minGQ) {
          g[i] = indv->justGet(GTidx).getGenotype();
        } else {
          g[i] = MISSING_GENOTYPE;
        }
      } else {
        g[i] = indv->justGet(GTidx).getGenotype();
      }
    }
    int diff;
    for (int i = 0; i < num; ++i) {
      for (int j = 0; j < i; ++j) {  // only use half of the matrix
        if (g[i] < 0 || g[j] < 0 || g[i] > 2 || g[j] > 2) {
          p9[i][j]++;
          continue;
        } else {
          diff = (int)abs(g[i] - g[j]);
          switch (diff) {
            case 0:
              p2[i][j]++;
              break;
            case 1:
              p1[i][j]++;
              break;
            case 2:
              p0[i][j]++;
              break;
            default:
              fprintf(stderr, "Something wrong!\n");
              exit(1);
          }
        }
      }  // end for j
    }    // end for i
    fprintf(fLog, "Processed %s:%d\n", r.getChrom(), r.getPos());
    // if (vout) vout->writeRecord(& r);
    // if (pout) pout ->writeRecord(& r);
  };

  // open output file
  foutName = FLAG_out + ".pairDiff";
  FILE* fp = fopen(foutName.c_str(), "wt");
  if (!fp) {
    fprintf(stderr, "Cannot open %s for write!\n", foutName.c_str());
    return -1;
  };

  for (int i = 0; i < num; ++i) {
    for (int j = 0; j < i; ++j) {  // only use half of the matrix
      fprintf(fp, "%s\t%s\t%g\t%g\t%g\t%g\n", name[i].c_str(), name[j].c_str(),
              p0[i][j], p1[i][j], p2[i][j], p9[i][j]);
    }
  }
  fclose(fp);

  fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);
  if (FLAG_variantOnly) {
    fprintf(stdout, "Skipped %d non-variant VCF records\n", nonVariantSite);
  }

  // if (vout) delete vout;
  // if (pout) delete pout;

  return 0;
};
