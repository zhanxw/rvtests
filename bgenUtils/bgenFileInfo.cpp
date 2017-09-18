#include <stdio.h>
#include <stdlib.h>

#include "base/Argument.h"
#include "libBgen/BGenFile.h"

#define DEBUG
#undef DEBUG

void printVariant(const BGenVariant& var, bool hideVarId, bool hideGT,
                  bool showDosage, FileWriter* fout) {
  fout->write(var.chrom.c_str());
  fout->write("\t");
  fout->printf("%d", var.pos);
  fout->write("\t");
  fout->write(var.rsid.c_str());
  if (!hideVarId && !var.varid.empty()) {
    fout->printf(",%s", var.varid.c_str());
  }
  fout->write("\t");
  fout->write(var.alleles[0].c_str());
  fout->write("\t");
  for (size_t i = 1; i != var.alleles.size(); ++i) {
    if (i != 1) {
      fout->write(",");
    }
    fout->write(var.alleles[i].c_str());
  }
  fout->write("\t.\t.\t.");

  if (var.isPhased) {  // phased info
    if (hideGT) {
      fout->write("\tHP");
    } else {
      fout->write("\tGT:HP");
    }
    if (showDosage) {
      fout->write(":DS");
    }
    for (size_t j = 0; j < var.N; ++j) {
      fout->write("\t");
      if (!hideGT) {
        var.printGT(j, fout);
        fout->write(":");
      }
      var.printHP(j, fout);
      if (showDosage) {
        fout->write(":");
        var.printDosage(j, fout);
      }
    }
  } else {  // genotypes
    if (hideGT) {
      fout->write("\tGP");
    } else {
      fout->write("\tGT:GP");
    }
    if (showDosage) {
      fout->write(":DS");
    }
    for (size_t j = 0; j < var.N; ++j) {
      fout->write("\t");
      if (!hideGT) {
        var.printGT(j, fout);
        fout->write(":");
      }
      var.printGP(j, fout);
      if (showDosage) {
        fout->write(":");
        var.printDosage(j, fout);
      }
    }
  }
  fout->write("\n");
}

//////////////////////////////////////////////////
// Parameter list
//////////////////////////////////////////////////
BEGIN_PARAMETER_LIST();
ADD_PARAMETER_GROUP("Basic Input/Output");
ADD_STRING_PARAMETER(inBgen, "--inBgen", "Input BGEN File");
ADD_STRING_PARAMETER(inSample, "--inSample", "Input Sample File");
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
  if (!FLAG_inSample.empty()) {
    if (read.loadSampleFile(FLAG_inSample)) {
      fprintf(stderr, "ERROR: failed to sample file [ %s ]!\n",
              FLAG_inSample.c_str());
      exit(1);
    }
  }
  read.printInfo();
  return 0;
}
