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

extern double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);

/**
 * PLINK allow these chromosomes (1-22, X, Y, XY, MT, 0)
 * according to http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml
 */
bool isPlinkCompatibleChrom(const std::string &chrom) {
  std::string c = toupper(chopChr(chrom));
  if (c == "X" || c == "Y" || c == "XY" || c == "MT") {
    return true;
  }
  int i;
  if (!str2int(c, &i) || !isdigit(c)) {  // not numeric
    return false;
  }
  if (0 <= i && i <= 25) {
    return true;
  }
  return false;
}
bool isACGT(char c) {
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T') return true;
  if (c == 'a' || c == 'c' || c == 'g' || c == 't') return true;
  return false;
}
/**
 * Check if the VCF site is a bi-allelic SNP site
 */
bool isBiallelicSite(const char *ref, const char *alt) {
  if (strlen(ref) != 1 || strlen(alt) != 1) return false;
  if (isACGT(ref[0]) && isACGT(alt[0])) return true;
  return false;
}

////////////////////////////////////////////////
BEGIN_PARAMETER_LIST()
ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(inVcf, "--inVcf", "input VCF File")
ADD_STRING_PARAMETER(outVcf, "--outVcf", "output prefix")
ADD_STRING_PARAMETER(outPlink, "--make-bed", "output prefix")
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
ADD_PARAMETER_GROUP("Gene Extractor")
ADD_STRING_PARAMETER(
    geneFile, "--geneFile",
    "Specify the gene file (refFlat format), so we know gene start and end.")
ADD_STRING_PARAMETER(geneName, "--gene", "Specify the gene names to extract")
ADD_STRING_PARAMETER(annoType, "--annoType",
                     "Specify the type of annotation to extract")
ADD_PARAMETER_GROUP("Quality Filter")
ADD_DOUBLE_PARAMETER(minSiteQual, "--minSiteQual", "Specify minimum site qual")
ADD_DOUBLE_PARAMETER(minGQ, "--minGQ",
                     "Specify the minimum genotype "
                     "quality, otherwise marked as "
                     "missing genotype")
ADD_DOUBLE_PARAMETER(
    minGD, "--minGD",
    "Specify the minimum genotype depth, otherwise marked as missing genotype")
ADD_PARAMETER_GROUP("Filter Option")
ADD_BOOL_PARAMETER(passFilter, "--pass",
                   "Only output variants that pass filters")
ADD_PARAMETER_GROUP("Quality Controls")
ADD_BOOL_PARAMETER(biallelic, "--biallelic",
                   "Specify only include bi-allelic SNP sites")
ADD_BOOL_PARAMETER(
    plinkChrom, "--plinkChrom",
    "Specify only include chromosome 1-22, X, Y, MT as PLINK required")
ADD_DEFAULT_DOUBLE_PARAMETER(minMAF, -1, "--minMAF",
                             "Specify minimum accepted minor allele frequency")
ADD_DEFAULT_DOUBLE_PARAMETER(maxMAF, -1, "--maxMAF",
                             "Specify maximum accepted minor allele frequency")
ADD_DEFAULT_DOUBLE_PARAMETER(minHWE, -1, "--minHWE",
                             "Specify minimum accepted HWE p-value")
ADD_DEFAULT_DOUBLE_PARAMETER(minCallRate, -1, "--minCallRate",
                             "Specify minimum call rate threshold")
ADD_PARAMETER_GROUP("Other Function")
ADD_BOOL_PARAMETER(variantOnly, "--variantOnly",
                   "Only variant sites from the VCF file will be processed.")
ADD_STRING_PARAMETER(updateId, "--update-id",
                     "Update VCF sample id using "
                     "given file (column 1 and 2 are "
                     "old and new id).")
END_PARAMETER_LIST();

int main(int argc, char **argv) {
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

  const char *fn = FLAG_inVcf.c_str();
  VCFInputFile vin(fn);

  // set range filters here
  // e.g.
  // vin.setRangeList("1:69500-69600");
  vin.setRangeList(FLAG_rangeList.c_str());
  vin.setRangeFile(FLAG_rangeFile.c_str());
  vin.setSiteFile(FLAG_siteFile.c_str());
  // set people filters here
  if (FLAG_peopleIncludeID.size() || FLAG_peopleIncludeFile.size()) {
    vin.excludeAllPeople();
    vin.includePeople(FLAG_peopleIncludeID.c_str());
    vin.includePeopleFromFile(FLAG_peopleIncludeFile.c_str());
  }
  vin.excludePeople(FLAG_peopleExcludeID.c_str());
  vin.excludePeopleFromFile(FLAG_peopleExcludeFile.c_str());

  // let's write it out.
  VCFOutputFile *vout = NULL;
  PlinkOutputFile *pout = NULL;
  if (FLAG_outVcf.size() > 0) {
    vout = new VCFOutputFile(FLAG_outVcf.c_str());
  };
  if (FLAG_outPlink.size() > 0) {
    pout = new PlinkOutputFile(FLAG_outPlink.c_str());
  };
  if (!vout && !pout) {
    vout = new VCFOutputFile("temp.vcf");
  }

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

  int lowSiteFreq = 0;  // counter of low site qualities
  int lowMAF = 0;
  int highMAF = 0;
  int lowCallRate = 0;
  int lowHWE = 0;
  int nonSnp = 0;  // counter of non-SNP site
  int incompatibleChrom =
      0;  // counter of chrom names that are incompatible to PLINK
  // int lowGDFreq = 0;
  // int lowGQFreq = 0;
  const std::string PassFilter = "PASS";
  // real working park
  if (vout) vout->writeHeader(vin.getVCFHeader());
  if (pout) pout->writeHeader(vin.getVCFHeader());
  int lineNo = 0;
  int nonVariantSite = 0;
  while (vin.readRecord()) {
    lineNo++;
    VCFRecord &r = vin.getVCFRecord();
    VCFPeople &people = r.getPeople();
    VCFIndividual *indv;
    if (FLAG_passFilter && PassFilter != r.getFilt()) continue;

    // check if only keep PLINK compatible chromosome names
    if (FLAG_plinkChrom && !isPlinkCompatibleChrom(r.getChrom())) {
      ++incompatibleChrom;
    }
    // check if the site is snp
    if (FLAG_biallelic && !isBiallelicSite(r.getRef(), r.getAlt())) {
      ++nonSnp;
      continue;
    }

    // examine individual genotype
    if (FLAG_variantOnly || FLAG_minMAF >= 0. || FLAG_maxMAF >= 0. ||
        FLAG_minHWE >= 0. || FLAG_minCallRate >= 0.0) {
      bool hasVariant = false;
      int geno;
      int homRef = 0, het = 0, homAlt = 0, missing = 0;
      int GTidx = r.getFormatIndex("GT");
      for (size_t i = 0; i < people.size(); i++) {
        indv = people[i];
        geno = indv->justGet(GTidx).getGenotype();
        if (geno != 0 && geno != MISSING_GENOTYPE) hasVariant = true;
        switch (geno) {
          case 0:
            ++homRef;
            break;
          case 1:
            ++het;
            break;
          case 2:
            ++homAlt;
            break;
          default:
            break;
        }
      }
      int total = homRef + het + homAlt + missing;
      if (!hasVariant) {  // check variant site
        nonVariantSite++;
        continue;
      }
      if (FLAG_minMAF >= 0. || FLAG_maxMAF >= 0.) {
        double maf = 0.0;
        if (homRef + het + homAlt > 0) {
          maf = 0.5 * (het + homAlt + homAlt) / (homRef + het + homAlt);
        }
        if (FLAG_minMAF >= 0. && maf < FLAG_minMAF) {
          ++lowMAF;
          continue;
        }
        if (FLAG_maxMAF >= 0. && maf > FLAG_maxMAF) {
          ++highMAF;
          continue;
        }
      }
      if (FLAG_minHWE >= 0.) {
        double hwe = 0.0;
        if (homRef > 0 || het > 0 || homAlt > 0) {
          hwe = SNPHWE(het, homAlt, homRef);
        }
        if (hwe < FLAG_minHWE) {
          ++lowHWE;
          continue;
        }
      }
      if (FLAG_minCallRate >= 0.) {
        double cr = 0.0;
        if (total - missing > 0) {
          cr = 1.0 - 1.0 * missing / total;
        }
        if (cr < FLAG_minCallRate) {
          ++lowCallRate;
          continue;
        }
      }
    }
    if (FLAG_minSiteQual > 0 && r.getQualDouble() < FLAG_minSiteQual) {
      ++lowSiteFreq;
      continue;
    }
    if (FLAG_annoType.size()) {
      bool isMissing = false;
      const char *tag = r.getInfoTag("ANNO", &isMissing).toStr();
      if (isMissing) continue;
      // fprintf(stdout, "ANNO = %s", tag);
      bool match = regex.match(tag);
      // fprintf(stdout, " %s \t", match ? "match": "noMatch");
      // fprintf(stdout, " %s \n", exists ? "exists": "missing");
      if (!match) {
        continue;
      }
    }
    if (FLAG_minHWE > 0) {
    }
    if (FLAG_minGD > 0 || FLAG_minGQ > 0) {
      if (vout) {
        vout->writeRecordWithFilter(&r, FLAG_minGD, FLAG_minGQ);
      }
      if (pout) {
        pout->writeRecordWithFilter(&r, FLAG_minGD, FLAG_minGQ);
      }
    } else {
      if (vout) vout->writeRecord(&r);
      if (pout) pout->writeRecord(&r);
    }
  };
  fprintf(stdout, "Total %d VCF records have converted successfully\n", lineNo);

  if (incompatibleChrom) {
    fprintf(stdout,
            "Skipped %d variants that are not located on PLINK "
            "compatible chromosomes (1-22, X, Y, XY, MT, 0)\n",
            incompatibleChrom);
  }

  if (nonSnp) {
    fprintf(stdout, "Skipped %d non-SNP VCF records\n", nonSnp);
  }

  if (nonVariantSite) {
    fprintf(stdout, "Skipped %d non-variant VCF records\n", nonVariantSite);
  }
  if (lowSiteFreq) {
    fprintf(stdout, "Skipped %d sites due to site quality lower than %f\n",
            lowSiteFreq, FLAG_minSiteQual);
  }
  if (lowMAF) {
    fprintf(stdout,
            "Skipped %d sites due to Minor Allele Frequency "
            "lower than %f\n",
            lowMAF, FLAG_minMAF);
  }
  if (highMAF) {
    fprintf(stdout,
            "Skipped %d sites due to Minor Allele Frequency "
            "higher than %f\n",
            highMAF, FLAG_maxMAF);
  }
  if (lowHWE) {
    fprintf(stdout, "Skipped %d sites due to HWE P-values lower than %f\n",
            lowHWE, FLAG_minHWE);
  }
  if (lowCallRate) {
    fprintf(stdout, "Skipped %d sites due to call rate lower than %f\n",
            lowCallRate, FLAG_minCallRate);
  }
  if (vout) delete vout;
  if (pout) delete pout;

  return 0;
};

// obtained from http://www.sph.umich.edu/csg/abecasis/Exact/snp_hwe.c
// 2012-11-29 Xiaowei
/*
  // This code implements an exact SNP test of Hardy-Weinberg Equilibrium as
  described in
  // Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
  // Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 -
  000
  //
  // Written by Jan Wigginton
  */

/**
 * NOTE (by zhanxw)
 * !! Makesure not all parameters equal to 0, or the program will crash.
 */
double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2) {
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) {
    printf(
        "FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) "
        "includes a"
        " negative count",
        obs_hets, obs_hom1, obs_hom2);
    exit(EXIT_FAILURE);
  }

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes = obs_hets + obs_homc + obs_homr;

  double *het_probs =
      (double *)malloc((size_t)(rare_copies + 1) * sizeof(double));
  if (het_probs == NULL) {
    printf(
        "FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote "
        "probabilities");
    exit(EXIT_FAILURE);
  }

  int i;
  for (i = 0; i <= rare_copies; i++) het_probs[i] = 0.0;

  /* start at midpoint */
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

  /* check to ensure that midpoint and rare alleles have same parity */
  if ((rare_copies & 1) ^ (mid & 1)) mid++;

  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
    het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets *
                               (curr_hets - 1.0) /
                               (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
    sum += het_probs[curr_hets - 2];

    /* 2 fewer heterozygotes for next iteration -> add one rare, one common
     * homozygote */
    curr_homr++;
    curr_homc++;
  }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
    het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr *
                               curr_homc /
                               ((curr_hets + 2.0) * (curr_hets + 1.0));
    sum += het_probs[curr_hets + 2];

    /* add 2 heterozygotes for next iteration -> subtract one rare, one common
     * homozygote */
    curr_homr--;
    curr_homc--;
  }

  for (i = 0; i <= rare_copies; i++) het_probs[i] /= sum;

  /* alternate p-value calculation for p_hi/p_lo
        double p_hi = het_probs[obs_hets];
           for (i = obs_hets + 1; i <= rare_copies; i++)
                p_hi += het_probs[i];

                   double p_lo = het_probs[obs_hets];
                      for (i = obs_hets - 1; i >= 0; i--)
                            p_lo += het_probs[i];


                               double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 *
     p_lo;
  */

  double p_hwe = 0.0;
  /*  p-value calculation for p_hwe  */
  for (i = 0; i <= rare_copies; i++) {
    if (het_probs[i] > het_probs[obs_hets]) continue;
    p_hwe += het_probs[i];
  }

  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}
