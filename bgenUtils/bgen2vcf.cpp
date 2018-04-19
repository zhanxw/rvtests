#include <stdio.h>
#include <stdlib.h>

#include "base/Argument.h"
#include "base/IO.h"
#include "libBgen/BGenFile.h"

#define DEBUG
#undef DEBUG

void printVCFMeta(FileWriter* fout) {
  // GP is between 0 and 1 in VCF v4.3, but phred-scaled value in VCF v4.2
  fout->write("##fileformat=VCFv4.3\n");
  if (false) {
    // it is not straightforward to know which tag to output without examining
    // all variant blocks
    fout->write(
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype "
        "calls\">\n");
    fout->write(
        "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype call "
        "probabilities\">\n");
    fout->write(
        "##FORMAT=<ID=HP,Number=.,Type=Float,Description=\"Haplotype "
        "probabilities\">\n");
  }
}

void printVCFHeader(const BGenFile& read, const std::vector<std::string>& sm,
                    FileWriter* fout) {
  const size_t sampleSize = read.getNumEffectiveSample();
  fout->write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

  for (size_t i = 0; i != sampleSize; ++i) {
    fout->printf("\t%s", sm[read.getEffectiveIndex(i)].c_str());
  }
  fout->write("\n");
}

void printVariant(const BGenFile& read, bool hideVarId, bool hideGT,
                  bool showDosage, FileWriter* fout) {
  const size_t sampleSize = read.getNumEffectiveSample();
  const BGenVariant& var = read.getVariant();

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
    for (size_t jIdx = 0; jIdx < sampleSize; ++jIdx) {
      const int j = read.getEffectiveIndex(jIdx);
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
    for (size_t jIdx = 0; jIdx < sampleSize; ++jIdx) {
      const int j = read.getEffectiveIndex(jIdx);
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
ADD_STRING_PARAMETER(inBgenSample, "--inBgenSample",
                     "Input SAMPLE file for the BGEN File");
ADD_STRING_PARAMETER(outPrefix, "--out", "Output prefix");

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
ADD_STRING_PARAMETER(
    siteFile, "--siteFile",
    "Specify the file containing site to extract, please use chr:pos format.")

ADD_PARAMETER_GROUP("Specify Options");
ADD_BOOL_PARAMETER(hideVarId, "--hideVarId",
                   "Do not output Variant ID (only output rsid)");
ADD_BOOL_PARAMETER(hideGT, "--hideGT",
                   "Do not call genotypes by skipping the GT tag");
ADD_BOOL_PARAMETER(showDS, "--showDS",
                   "Calculate bi-allelic dosage using the DS tag");
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

  std::string outFileName = FLAG_outPrefix + ".vcf.gz";
  FileWriter fout(outFileName, BGZIP);
  BGenFile read(FLAG_inBgen);
  if (!FLAG_inBgenSample.empty()) {
    read.loadSampleFile(FLAG_inBgenSample);
  }

  const int N = read.getNumSample();
  const int M = read.getNumMarker();
  std::vector<std::string> sm = read.getSampleIdentifier();
  if (sm.empty()) {
    fprintf(stderr,
            "Cannot find sample names from BGEN input file. Sample names wth "
            "be set at sample_0, sample_1, ...\n");
    char buf[1024];
    for (int i = 0; i < N; ++i) {
      sprintf(buf, "sample_%d", i);
      sm.push_back(buf);
    }
  }

  read.setRangeList(FLAG_rangeList.c_str());
  read.setRangeFile(FLAG_rangeFile.c_str());
  read.setSiteFile(FLAG_siteFile.c_str());
  // set people filters here
  if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
    read.excludeAllPeople();
    read.includePeople(FLAG_peopleIncludeID.c_str());
    read.includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
  }
  read.excludePeople(FLAG_peopleExcludeID.c_str());
  read.excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());

  printVCFMeta(&fout);
  printVCFHeader(read, sm, &fout);

  fprintf(stderr, "BGEN File has [ %d ] samples, [ %d ] markers\n", N, M);
  fprintf(stderr, "Effective sample size is [ %d ]\n",
          read.getNumEffectiveSample());
  while (read.readRecord()) {
    printVariant(read, FLAG_hideVarId, FLAG_hideGT, FLAG_showDS, &fout);
  }  // loop marker

  fprintf(stderr, "Total %d sample and %d variatns processed\n", N, M);
  fprintf(stderr, "Conversion succeed!\n");
  return 0;
}
