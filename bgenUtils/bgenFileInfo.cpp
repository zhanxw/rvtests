#include <stdio.h>
#include <stdlib.h>

#include "base/Argument.h"
#include "libBgen/BGenFile.h"

#define DEBUG
#undef DEBUG

void printVariant(const BGenVariant& var, bool hideVarId, bool hideGT,
                  bool showDosage, FILE* fout) {
  fputs(var.chrom.c_str(), fout);
  fputs("\t", fout);
  fprintf(fout, "%d", var.pos);
  fputs("\t", fout);
  fputs(var.rsid.c_str(), fout);
  if (!hideVarId && !var.varid.empty()) {
    fprintf(fout, ",%s", var.varid.c_str());
  }
  fputs("\t", fout);
  fputs(var.alleles[0].c_str(), fout);
  fputs("\t", fout);
  for (size_t i = 1; i != var.alleles.size(); ++i) {
    if (i != 1) {
      fputs(",", fout);
    }
    fputs(var.alleles[i].c_str(), fout);
  }
  fputs("\t.\t.\t.", fout);

  if (var.isPhased) {  // phased info
    if (hideGT) {
      fputs("\tHP", fout);
    } else {
      fputs("\tGT:HP", fout);
    }
    if (showDosage) {
      fputs(":DS", fout);
    }
    for (size_t j = 0; j < var.N; ++j) {
      fputs("\t", fout);
      if (!hideGT) {
        var.printGT(j, fout);
        fputs(":", fout);
      }
      var.printHP(j, fout);
      if (showDosage) {
        fputs(":", fout);
        var.printDosage(j, fout);
      }
    }
  } else {  // genotypes
    if (hideGT) {
      fputs("\tGP", fout);
    } else {
      fputs("\tGT:GP", fout);
    }
    if (showDosage) {
      fputs(":DS", fout);
    }
    for (size_t j = 0; j < var.N; ++j) {
      fputs("\t", fout);
      if (!hideGT) {
        var.printGT(j, fout);
        fputs(":", fout);
      }
      var.printGP(j, fout);
      if (showDosage) {
        fputs(":", fout);
        var.printDosage(j, fout);
      }
    }
  }
  fputs("\n", fout);
}

//////////////////////////////////////////////////
// Parameter list
//////////////////////////////////////////////////
BEGIN_PARAMETER_LIST();
ADD_PARAMETER_GROUP("Basic Input/Output");
ADD_STRING_PARAMETER(inBgen, "--inBgen", "Input BGEN File");
ADD_BOOL_PARAMETER(help, "--help", "Print detailed help message");
END_PARAMETER_LIST();

int main(int argc, char* argv[]) {
  PARSE_PARAMETER(argc, argv);

  if (FLAG_help) {
    PARAMETER_HELP();
    return 0;
  }

  PARAMETER_STATUS();
  if (FLAG_REMAIN_ARG.size() > 0) {
    fprintf(stderr, "Unparsed arguments: ");
    for (unsigned int i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
      fprintf(stderr, " %s", FLAG_REMAIN_ARG[i].c_str());
    }
    exit(1);
  }

  REQUIRE_STRING_PARAMETER(FLAG_inBgen,
                           "Please provide input file using: --inBgen");

  BGenFile read(FLAG_inBgen);
  read.printInfo();
  return 0;
}
