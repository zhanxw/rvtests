#include <stdio.h>
#include <stdlib.h>

#include "base/Argument.h"
#include "libBgen/BGenFile.h"

// #include <stdint.h>  // uint32_t
// #include <cassert>
// #include <string>
// #include <vector>

// #include <sqlite3.h>
// #include <zlib.h>
// #include "zstd/lib/zstd.h"

#define DEBUG
#undef DEBUG

void printVCFMeta(FILE* fout) {
  // GP is between 0 and 1 in VCF v4.3, but phred-scaled value in VCF v4.2
  fprintf(fout, "##fileformat=VCFv4.3\n");
  if (false) {
    // it is not straightforward to know which tag to output without examining
    // all variant blocks
    fprintf(fout,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype "
            "calls\">\n");
    fprintf(fout,
            "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype call "
            "probabilities\">\n");
    fprintf(fout,
            "##FORMAT=<ID=HP,Number=.,Type=Float,Description=\"Haplotype "
            "probabilities\">\n");
  }
}

void printVCFHeader(const std::vector<std::string>& sm, FILE* fout) {
  fprintf(fout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  for (size_t i = 0; i != sm.size(); ++i) {
    fprintf(fout, "\t%s", sm[i].c_str());
  }
  fprintf(fout, "\n");
}

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
ADD_STRING_PARAMETER(outPrefix, "--out", "Output prefix");

ADD_PARAMETER_GROUP("Specify Options");
ADD_BOOL_PARAMETER(hideVarId, "--hideVarId", "Do not output Variant ID");
ADD_BOOL_PARAMETER(hideGT, "--hideGT", "Do not call genotypes");
ADD_BOOL_PARAMETER(showDS, "--showDS", "Print bi-allelic dosage");
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

  if (!FLAG_outPrefix.size()) FLAG_outPrefix = "rvtest";

  REQUIRE_STRING_PARAMETER(FLAG_inBgen,
                           "Please provide input file using: --inBgen");

  FILE* fout = fopen((FLAG_outPrefix + ".vcf").c_str(), "wt");
  BGenFile read(FLAG_inBgen);
  const int N = read.getNumSample();
  const int M = read.getNumMarker();
  const std::vector<std::string> sm = read.getSampleIdentifier();
  printVCFMeta(fout);
  printVCFHeader(sm, fout);

  const BGenVariant& var = read.getVariant();

  while (read.readRecord()) {
    printVariant(var, FLAG_hideVarId, FLAG_hideGT, FLAG_showDS, fout);
  }  // loop marker

  fclose(fout);
  printf("Sample = %d, #SampleIdentifier = %d, marker = %d\n", N,
         (int)sm.size(), M);
  printf("Conversion succeed!\n");
  return 0;
}
