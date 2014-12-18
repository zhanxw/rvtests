//: src:kbac.cpp
// kbac test
// Copyright 2011 Gao Wang
// Modified 2014 Xiaowei Zhan
#include "kbac.h"
//#include <stdio.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_gamma.h"

namespace {
  const double AFFECTED = 1.0, UNAFFECTED = 0.0, HOMO_ALLELE = 2.0, 
        MINOR_ALLELE = 1.0, MAJOR_ALLELE = 0.0, MISSING_ALLELE = -9.0;
}

KbacTest::KbacTest (int* nn, int* qq, double* aa, double* mafUpper, double* xdatIn, double* ydatIn, double* mafIn, int* xcol, int* ylen) 
{

  this->__nPermutations = *nn;
  this->__quiet = *qq;
  this->__alpha = *aa;
  this->__mafUpper = *mafUpper;
  this->__mafLower = 0.0;
  
  if (__alpha >= 1.0) {
    __adaptive = 0;
  }
  else {
    __adaptive = 5000;
  }

  for (int i = 0; i != *ylen; ++i) { 
    __ydat.push_back(ydatIn[i]);
  }
  for (int i = 0; i != *xcol; ++i) { 
    __observedMafs.push_back(mafIn[i]);
  }

  int start = 0;
  int nInvalidGenotypes = 0;
  for (int i = 0; i != *ylen; ++i) {
    std::vector<double> tmpvct(0);
    for (int j = 0; j != *xcol; ++j) {
      if (xdatIn[start+j] != MAJOR_ALLELE && xdatIn[start+j] != MINOR_ALLELE && xdatIn[start+j] != HOMO_ALLELE) {
        ++nInvalidGenotypes;
        tmpvct.push_back(MAJOR_ALLELE);
      }
      else {
        tmpvct.push_back(xdatIn[start+j]);
        // if (xdatIn[start+j])
        //   std::cout << "xdatIn[start+j]" << start+j << " " <<xdatIn[start+j] << "\n";
      }
    }
    __xdat.push_back(tmpvct);
    start += *xcol;
  }
     
  if (nInvalidGenotypes > 0) {
    std::cout << "**Warning: Invalid genotype codings detected on " << nInvalidGenotypes << " variant sites (codings must be 0 = wild-type, 1 = heterozygous, or 2 = homozygous. KBAC will treat them as 0 = wild-type" << std::endl;
  }

  if (__quiet == false) {
    std::cout << "Read data Y length " << __ydat.size() << " and X dimension " << __xdat.size() << " by " << __xdat[0].size() << std::endl;
  }

  bool isInputOk = (start == *xcol * *ylen && __xdat.size() != 0 && __ydat.size() == __xdat.size() 
      && __observedMafs.size() == __xdat[0].size() && __mafLower >= 0.0 
      && __mafUpper <= 1.0 && __mafUpper > __mafLower && __alpha > 0.0);

  if (isInputOk);
  else {
    std::cerr << "**Error** Input data problem in KbacTest::KbacTest(). Now Quit." << std::endl;
    exit(-1);
  }
}


KbacTest::~KbacTest() {};


void KbacTest::calcKbacP(double* pvalue, int* sided)
{
  /*! * the KBAC Statistic: sum of genotype pattern frequencies differences, weighted by hypergeometric kernel. <br>
   * * It is a permutation based two-sided test.  <br>
   * * See <em> Liu DJ 2010 PLoS Genet. </em> <br><br>
   *Implementation:
   */
  if (*sided == 0 || *sided > 2)
    *sided = 1;

  //!- trim data by mafs upper-lower bounds
  m_trimXdat();


  unsigned int sampleSize = __ydat.size();
  // sample size
  unsigned int regionLen = __xdat[0].size();
  // candidate region length
  unsigned int nCases = 0; 
  // case size

  for(unsigned int i = 0; i !=  __ydat.size(); ++i) {
    if(__ydat[i] == AFFECTED) 
      ++nCases;
  }
  unsigned int nCtrls = sampleSize - nCases;

  if (__quiet == false) {
    std::cout << "Sample size: " << sampleSize << std::endl;
    std::cout << "Number of cases: " << nCases << std::endl;
    std::cout << "Candidate region length: " << regionLen << std::endl;
  }

  std::vector<double> genotypeId(sampleSize);   
  //!-Compute unique genotype patterns (string) as ID scores (double) 

  for (unsigned int i = 0; i != sampleSize; ++i) {

    double vntIdL = 0.0; 
    double vntIdR = 0.0;
    const double ixiix= pow(9.0, 10.0);
    unsigned int lastCnt = 0;
    unsigned int tmpCnt = 0;

    for (unsigned int j = 0; j != regionLen; ++j) { 

      if (__xdat[i][j] != MISSING_ALLELE && __xdat[i][j] != MAJOR_ALLELE) 
        vntIdR += pow(3.0, 1.0 * (j - lastCnt)) * __xdat[i][j];
      else 
        continue;
      if (vntIdR >= ixiix) {
        vntIdL = vntIdL + 1.0;
        vntIdR = vntIdR - ixiix;
        lastCnt = lastCnt + tmpCnt + 1;
        tmpCnt = 0;
        continue;
      }
      else { 
        ++tmpCnt; 
        continue; 
      }
    }

    genotypeId[i] = vntIdL + vntIdR * 1e-10;
    // one-to-one "ID number" for a genotype pattern
  }

  // unique genotype patterns
  // use "list" and some STL algorithms

  std::list<double> uniqueId(genotypeId.begin(), genotypeId.end());
  uniqueId.remove(MAJOR_ALLELE);   

  if (uniqueId.size() == 0) {
    std::cout << "**Warning** non-wildtype genotype data is empty. KBAC has nothing to work on. Return p-value 1.0" << std::endl;
    *pvalue = 1.0;
    return;
  }

  uniqueId.sort(); 
  uniqueId.unique();
  // remove wildtype and get unique genotype patterns 
  //!- Unique genotype patterns that occur in the sample


  std::vector<double> uniquePattern(uniqueId.size());
  copy(uniqueId.begin(), uniqueId.end(), uniquePattern.begin());
  uniqueId.clear();

  if (__quiet == false) {
    std::cout << "All individual genotype patterns: " << "\n" << std::endl;
    std::cout << genotypeId << std::endl;
    std::cout << "Unique individual genotype patterns: " << "\n" << std::endl;
    std::cout << uniquePattern << std::endl;
  }

  // count number of sample individuals for each genotype pattern
  unsigned int uniquePatternCounts[uniquePattern.size()];
  for (unsigned int u = 0; u != uniquePattern.size(); ++u) 
    uniquePatternCounts[u] = 0;

  for (unsigned int i = 0; i != sampleSize; ++i) {
    // for each sample, identify/count its genotype pattern

    for (unsigned int u = 0; u != uniquePattern.size(); ++u) {

      if (genotypeId[i] == uniquePattern[u]) {
        // genotype pattern identified
        ++uniquePatternCounts[u];
        // count this genotype pattern
        break;
      }
      else;
      // genotype pattern not found -- move on to next pattern
    }
  }

  if (__quiet == false) {
    std::cout << "Number of each unique individual genotype patterns (totaling " << uniquePattern.size() << " patterns excluding wildtype): " << "\n" << std::endl;
    for (unsigned int u = 0; u != uniquePattern.size(); ++u) std::cout << uniquePatternCounts[u] << ", "; 
    std::cout << "\n" << std::endl;
  }

  unsigned int iPermutation = 0;
  unsigned int permcount1 = 0, permcount2 = 0;
  double observedStatistic = 0.0;
  while (iPermutation <= __nPermutations) {

    // the KBAC statistic. Will be of length 1 or 2
    std::vector<double> kbacStatistics(0);
    // two models
    for (int s = 0; s != *sided; ++s) {

      //!- count number of sample cases (for the 1st model, or ctrls for the 2nd model) for each genotype pattern
      unsigned int uniquePatternCountsSub[uniquePattern.size()];
      for (unsigned int u = 0; u != uniquePattern.size(); ++u) 
        uniquePatternCountsSub[u] = 0;
      // genotype pattern counts in cases (for the 1st model, or ctrls for the 2nd model) 

      for (unsigned int i = 0; i != sampleSize; ++i) {
        if ( __ydat[i] == (AFFECTED - 1.0 * s) ) {
          // for each "case (for the 1st model, or ctrls for 2nd model)", identify/count its genotype pattern
          for (unsigned int u = 0; u != uniquePattern.size(); ++u) {
            if (genotypeId[i] == uniquePattern[u]) {
              // genotype pattern identified in cases (for the 1st model, or ctrls for 2nd model)
              ++uniquePatternCountsSub[u];
              // count this genotype pattern
              break;
            }
            else;
            // genotype pattern not found -- move on to next pattern
          }
        }
        else;
      }

      //!- KBAC weights
      double uniquePatternWeights[uniquePattern.size()];
      // genotype pattern weights, the hypergeometric distribution cmf
      for (unsigned int u = 0; u != uniquePattern.size(); ++u) 
        uniquePatternWeights[u] = 0.0;

      for (unsigned int u = 0; u != uniquePattern.size(); ++u) {
        if (s == 0) 
          uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);
          //uniquePatternWeights[u] = gw_hypergeometric_cmf(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCases);
        else
          uniquePatternWeights[u] = gsl_cdf_hypergeometric_P(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCtrls);
          //uniquePatternWeights[u] = gw_hypergeometric_cmf(uniquePatternCountsSub[u], uniquePatternCounts[u], sampleSize - uniquePatternCounts[u], nCtrls);
      }

      //!- KBAC statistic: sum of genotype pattern frequencies differences in cases vs. controls, weighted by the hypergeometric distribution kernel
      double kbac = 0.0;
      for (unsigned int u = 0; u != uniquePattern.size(); ++u) { 
        if (s == 0)
          kbac = kbac + ( (1.0 * uniquePatternCountsSub[u]) / (1.0 * nCases) - (1.0 * (uniquePatternCounts[u] - uniquePatternCountsSub[u])) / (1.0 * nCtrls) ) *  uniquePatternWeights[u];
        else
          kbac = kbac + ( (1.0 * uniquePatternCountsSub[u]) / (1.0 * nCtrls) - (1.0 * (uniquePatternCounts[u] - uniquePatternCountsSub[u])) / (1.0 * nCases) ) *  uniquePatternWeights[u];
      }

      //std::cout << kbac << std::endl;

      //FIXME
      //gw_round(kbac, 0.0001);
      kbacStatistics.push_back(kbac);
    }

    double statistic = 0.0;
    //!- one model statistic
    if (kbacStatistics.size() == 1) {
      statistic = kbacStatistics[0];
    }
    //!- two model statistic
    else if (kbacStatistics.size() == 2) {
      statistic = fmax(kbacStatistics[0], kbacStatistics[1]);
    }
    else {
      std::cerr << "**Error KBAC statistic (Error code -5)" << std::endl;
      exit(-1);
    }

    if (iPermutation == 0) 
      observedStatistic = statistic;
    else {
      if (statistic >= observedStatistic) 
        ++permcount1;
      if (statistic <= observedStatistic)
        ++permcount2;
      if (__adaptive != 0)
        *pvalue = m_checkAdaptivePvalue(permcount1, permcount2, iPermutation, __adaptive, 0);
    }
    if (*pvalue <= 1.0)
      break;
    //!- Permutation
    random_shuffle(__ydat.begin(), __ydat.end());
    ++iPermutation;
  }
  if (*pvalue <= 1.0);
  else {
    *pvalue = (1.0 * permcount1 + 1.0) / (1.0 * __nPermutations + 1.0);
  }
  return;
}


double KbacTest::m_checkAdaptivePvalue(unsigned int permcount1, unsigned int permcount2, unsigned int currentIdx, unsigned int checkPoint, unsigned int alternative) const
{
  if (currentIdx % checkPoint == 0 && checkPoint > 5) {
    //!- adaptive p-value calculation, at an interval of #checkPoint permutations 
    // apply the "six-sigma" rule

    double adaptivePvalue = 1.0;
    if (alternative == 1 || alternative == 0) { 
      adaptivePvalue = (1.0 * permcount1 + 1.0) / (1.0 * currentIdx + 1.0);
    }
    else {
      double permcount = (permcount1 < permcount2) ? permcount1 : permcount2;
      adaptivePvalue = (2.0 * permcount + 1.0) / (1.0 * currentIdx + 1.0);
    }

    double sd = sqrt(adaptivePvalue * (1.0 - adaptivePvalue) / (1.0 * currentIdx));
    double sixsigma = adaptivePvalue - 6.0 * sd;

    if (__quiet == false) {
      std::cout << "nPerm" << "\t" << "Current.P" << "\t" << "St.Error" << "\t" << "6-sigma.lower.bound" << std::endl;
      std::cout << currentIdx << "\t" << adaptivePvalue << "\t" << sd << "\t" << sixsigma << std::endl;
    }

    if (sixsigma > __alpha) 
      return adaptivePvalue;
    else
      return 9.0;
  }
  else 
    return 9.0;
}


void KbacTest::m_trimXdat()
{
  std::vector< std::vector<double> > xdat = __xdat;
  __xdat.clear();
  __xdat.resize(xdat.size());

  for (unsigned int j = 0; j != __observedMafs.size(); ++j) {
    if (__observedMafs[j] <= __mafLower || __observedMafs[j] > __mafUpper) 
      continue;

    else {
      for (unsigned int i = 0; i != xdat.size(); ++i)
        __xdat[i].push_back(xdat[i][j]);
    }
  }  
  return;
}
