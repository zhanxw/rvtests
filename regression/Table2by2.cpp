/*
 *  Table2by2.cpp
 *  LogitReg
 *
 *  Created by Youna Hu on 11/21/10.
 *  Modified by Xiaowei Zhan 2014
 *  Copyright 2010 U of Michigan. All rights reserved.
 *
 */

#include "Table2by2.h"
#include "Error.h"
#include <cmath>
#include <cstring>
#include <gsl/gsl_sf_gamma.h>


Table2by2::Table2by2()
{
  // initialize the cell count to be 0
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      m_n[i][j] = 0;
    }
  }
  // update the marginal columns
  for (int i = 0; i < 2; i ++) {
    m_marginCol[i] = 0;
    for (int j = 0; j < 2; j ++) {
      m_marginCol[i] += m_n[i][j];
    }
  }
  // Here, can I use the UpdateMarginSum() function instead???
  // update the marginal rows
  for (int j = 0; j < 2; j ++) {
    m_marginRow[j] = 0;
    for (int i = 0; i < 2; i ++) {
      m_marginCol[j] += m_n[i][j];
    }
  }
  m_sum = 0;
  InitOtherVariables();
}

void Table2by2::InitOtherVariables(){
  chisqStat = -1.0;
  pChisq = -1.0;
  pExact = -1.0;
  pGeneral = -1.0;
  pExactTwoSided = -1;
  pExactOneSidedLess = -1;
  pExactOneSidedGreater = -1;
  lowerBound = -1;
  upperBound = -1;
}

void Table2by2::InitializeFromNums(int n00, int n01, int n10, int n11){
  m_n[0][0] = n00;
  m_n[0][1] = n01;
  m_n[1][0] = n10;
  m_n[1][1] = n11;
}

Table2by2::Table2by2(int n00, int n01, int n10, int n11)
{
  InitializeFromNums(n00,n01,n10,n11);
  UpdateMarginSum();
  InitOtherVariables();
}

Table2by2::Table2by2(Vector &yIn, Vector &xIn)
{
  InitializeFromVectors(yIn, xIn);
  UpdateMarginSum();
  InitOtherVariables();
}

void Table2by2::reset() {
  m_n[0][0] = 0;
  m_n[0][1] = 0;
  m_n[1][0] = 0;
  m_n[1][1] = 0;
};

void Table2by2::Increment(int rowIndex, int colIndex)
{
  //      printf("m_n[%d,%d] = %d\n",rowIndex, colIndex, m_n[rowIndex][colIndex]);
  m_n[rowIndex][colIndex]  ++;
}

void Table2by2::UpdateMarginSum()
{
  // update the marginal columns
  // Row sum: Row[i] = sum(m_n[i][j])
  for (int i = 0; i < 2; i ++) {
    m_marginRow[i] = 0;
    for (int j = 0; j < 2; j ++) {
      m_marginRow[i] += m_n[i][j];
    }
  }

  // update the marginal rows and total sum
  // Col sum: Col[j] = sum(m_n[i][j])
  m_sum = 0;
  for (int j = 0; j < 2; j ++) {
    m_marginCol[j] = 0;
    for (int i = 0; i < 2; i ++) {
      m_marginCol[j] += m_n[i][j];
      m_sum += m_n[i][j];
    }
  }
}

void Table2by2::InitializeFromVectors(Vector &yIn, Vector &xIn)
{
  // output the 2 by 2 table with (case, ctrl) * (0,1)
  if (yIn.Length() != xIn.Length()) {
    error("Error: input vectors are of different length!\n");
  }
  int xyLen = yIn.Length();
  // here do I need to initialize the table again?
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      m_n[i][j] = 0;
    }
  }
  for (int i = 0; i < xyLen; i ++ ) {
    if (xIn[i] != 0.0) {
      Increment((int) (1-yIn[i]),1);
    } else {
      Increment((int) (1-yIn[i]) , (int) xIn[i]);
    }
  }
  UpdateMarginSum();
}

void Table2by2::InitializeFromVecMat(Vector &yIn, Matrix &xIn)
{
  // output the 2 by 2 table with (case, ctrl) * (0,1)
  if (yIn.Length() != xIn.rows) {
    error("Error: input length is not the same as matrix row number. \n");
  }
  int xyLen = yIn.Length();
  // here do I need to initialize the table again?
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      m_n[i][j] = 0;
    }
  }
  for (int i = 0; i < xyLen; i ++ ) {
    if (xIn[i][0] != 0.0) { // code nonzero xIn as 1
      Increment((int) (1-yIn[i]),1);
    } else {
      Increment((int) (1-yIn[i]) , (int) xIn[i][0]);
    }
  }
  UpdateMarginSum();
  chisqStat = 0.0;
  pChisq = -1.0;
  pExact = -1.0;
  pGeneral = -1.0;
}


bool Table2by2::MustExactTest()
{
  double expected2by2[2][2];
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      expected2by2[i][j] = (double)m_marginRow[i]* (double) m_marginCol[j]/(double) m_sum;
      if (expected2by2[i][j] < 5.0 || m_n[0][0] == 0 || m_n[0][1] == 0
          || m_n[1][0] == 0 || m_n[1][1] == 0) {
        return (true);
      }
    }
  }
  return (false);
}

void Table2by2::ChisqTest(bool & valid)
{
  valid = true;
  double expected2by2[2][2];
  UpdateMarginSum();
  //      PrintMargin();
  chisqStat = 0.0;
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      expected2by2[i][j] = (double)m_marginRow[i]* (double) m_marginCol[j]/(double) m_sum;
      if (expected2by2[i][j] < 5.0) {
        valid = false;
        return;
      }
      chisqStat       += ((double) m_n[i][j] - expected2by2[i][j])*( (double) m_n[i][j] - expected2by2[i][j])/expected2by2[i][j];
      //                      printf("expected2by2 = %f, %f\n",expected2by2[i][j], ((double) m_n[i][j] - expected2by2[i][j])*( (double) m_n[i][j] - expected2by2[i][j])/expected2by2[i][j]);
    }
  }
  // p value (df for IbyJ table the df is (I-1)*(J-1))
  //      printf("pChisq = %f,chisqStat = %f\n",pChisq,chisqStat);
  pChisq = chidist(chisqStat,1);
}

void Table2by2::PearsonChisq()
{
  double expected2by2[2][2];
  UpdateMarginSum();
  for     (int i = 0; i < 2; i ++){
    if (m_marginCol[i] == 0 || m_marginRow[i] == 0){
      chisqStat = 0;
      return;
    }
  }
  //  PrintMargin();
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      expected2by2[i][j] = (double)m_marginRow[i]* (double) m_marginCol[j]/ (double) m_sum;
      //                      printf("expected2by2 = %f, %f\n",expected2by2[i][j], ((double) m_n[i][j] - expected2by2[i][j])*( (double) m_n[i][j] - expected2by2[i][j])/expected2by2[i][j]);
      chisqStat       += ((double) m_n[i][j] - expected2by2[i][j])*( (double) m_n[i][j] - expected2by2[i][j])/expected2by2[i][j];
    }
  }
  //      printf("chisqStat = %f\n",chisqStat);
}


void Table2by2::FisherExactTest()
{
  //Update the margin and sum
  UpdateMarginSum();
  CalculateBoundsIn00ForFisher();
  // p value for observed = choose(n1+,n11)*choose(n2+,n21)/choose(n++,n+1)
  double pObserved = gsl_sf_choose(m_marginRow[0],m_n[0][0])*gsl_sf_choose(m_marginRow[1], m_n[1][0])/gsl_sf_choose(m_sum, m_marginCol[0]);
  pExact = 0;
  //      printf("n1+ = %d, n11 = %d, sum = %d, pObserved = %f\n",m_marginCol[0], m_n[0][0], m_sum, pObserved);
  //      printf("lowerbound = %d, upperbound = %d\n",lowerbound, upperBound);
  if (upperBound != lowerBound) {
    for (int i = lowerBound; i <= upperBound; i ++) {
      if (i != m_n[0][0]) {
        double pPossible = gsl_sf_choose(m_marginRow[0],i)*gsl_sf_choose(m_marginRow[1], m_marginCol[0] - i)/gsl_sf_choose(m_sum, m_marginCol[0]);
        //                      printf("i = %d, pPossible = %f\n",i, pPossible);
        if (pPossible <= pObserved) {
          pExact += pPossible;
        }
      }
    }
    pExact += pObserved;
  }
  else {
    pExact = pObserved;
  }
  //      printf("pExact = %f\n",pExact);
}

void Table2by2::GeneralTest()
{
  bool useChisqValid;
  ChisqTest(useChisqValid);
  if (!useChisqValid) {
    FisherExactTest();
    pGeneral = pExact;
  } else {
    pGeneral = pChisq;
  }

  /*      printf("valid to use pearson chisq test? %s\n",useChisqValid?"true":"false");
          printf("p value for general test = %f\n",pGeneral); */
}

void Table2by2::CalculateBoundsIn00ForFisher(){

  // uppderBound = min(n1+,n+1)
  upperBound = m_marginRow[0];
  if (upperBound > m_marginCol[0]) upperBound = m_marginCol[0];

  // lowerbound = max(0,n1+ + n+1 - sum)
  lowerBound = m_marginRow[0] + m_marginCol[0] - m_sum;
  if (lowerBound < 0) lowerBound = 0;
}

// this logFacs is used to store log(n!),
void Table2by2::initLogFacs(double *logFacs, int n){
  logFacs[0] = 0;
  for(int i=1; i < n+1; ++i) {
    logFacs[i] = logFacs[i-1] + log((double)i); // only n times of log() calls
  }
}

void Table2by2::initLogFacs(double *logFacs){
  logFacs[0] = 0;
  for (int i = 1; i < m_sum + 1; i ++) {
    logFacs[i] = logFacs[i-1] + log((double)i);
  }
}

double Table2by2::logHypergeometricProb(double *logFacs, int a, int b, int c, int d){
  return (logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d]
          - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d]);
}

double Table2by2::logHypergeometricProb(double *logFacs){
  int a = m_n[0][0], b = m_n[0][1], c = m_n[1][0], d = m_n[1][1];
  return (logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d]
          - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d]);
}

void Table2by2::FullFastFisherExactTest(){
  double* storeLogFacs = new double[m_sum + 1];

  // calculate and store all the factorials
  initLogFacs(storeLogFacs);

  // calcualte the observed p value
  double logpCutoff = logHypergeometricProb(storeLogFacs);
  double pFraction = 0;
  double pFractionLess = 0;
  double pFractionGreater = 0;

  // calculate the boundaries for m_n[0][0] given the four marginals
  CalculateBoundsIn00ForFisher();

  for (int i = lowerBound; i <= upperBound; i ++) {
    double logpPossible = logHypergeometricProb(storeLogFacs,
                                                i,m_marginRow[0] - i, m_marginCol[0] - i,
                                                m_marginRow[1] + i - m_marginCol[0]);
    if (logpPossible <= logpCutoff) {
      pFraction += exp(logpPossible - logpCutoff);
    }
    if (i <= m_n[0][0]) {
      pFractionLess += exp(logpPossible - logpCutoff);
    }
    if (i >= m_n[0][0]) {
      pFractionGreater += exp(logpPossible - logpCutoff);
    }
  }

  double logpValue, logpValueLess, logpValueGreater;

  logpValue = logpCutoff + log(pFraction);
  logpValueLess = logpCutoff + log(pFractionLess);
  logpValueGreater = logpCutoff + log(pFractionGreater);

  // update the p values
  pExactTwoSided = exp(logpValue);
  pExactOneSidedLess = exp(logpValueLess);
  pExactOneSidedGreater = exp(logpValueGreater);
  delete[] storeLogFacs;
}


// get values from the 2 by 2 table
int Table2by2::GetValue(int rowIndex, int colIndex)
{
  return(m_n[rowIndex][colIndex]);
}

double Table2by2::GetpChisq()
{
  return(pChisq);
}

double Table2by2::GetpExact()
{
  return(pExact);
}

double Table2by2::GetstatChisq()
{
  return(chisqStat);
}

double Table2by2::GetpGeneral()
{
  return(pGeneral);
}

void Table2by2::Print()
{
  printf("+++ 2 by 2 table +++\n");
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      printf("%d\t",m_n[i][j]);
    }
  }
  printf("\n");
}

void Table2by2::PrintMargin()
{
  for (int i = 0; i < 2; i ++) {
    printf("col margin: %d\t",m_marginCol[i]);
  }
  printf("\t");
  for (int i = 0; i < 2; i ++) {
    printf("row margin: %d\t",m_marginRow[i]);
  }
  printf("\n");
}

void Table2by2::GetTable( int b[2][2]){
  for (int i = 0; i < 2; i ++) {
    for (int j = 0; j < 2; j ++) {
      b[i][j] = m_n[i][j];
    }
  }
}

double Table2by2::GetpFasterFisherExact(const char* testIndicator){
  if (strcmp(testIndicator, "TwoSided") == 0) {
    return (pExactTwoSided);
  } else if (strcmp(testIndicator, "Less") == 0) {
    return (pExactOneSidedLess);
  } else if (strcmp(testIndicator, "Greater") == 0) {
    return (pExactOneSidedGreater);
  } else {
    error("Inappropriate test side specification! 0: two sided; 1: one sided greater; -1: one sided less\n");
  }

  return 0.0;
}
