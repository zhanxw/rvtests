/*
 *  Table2by2.h
 *  LogitReg
 *
 *  Created by Youna Hu on 11/21/10.
 *  Copyright 2010 U of Michigan. All rights reserved.
 *
 */
#ifndef __TABLE2by2_h__
#define __TABLE2by2_h__

#include "MathStats.h"

class Table2by2{
public:
	// constructors and initializers
	Table2by2();
	Table2by2(int n00, int n01, int n10, int n11);
    void reset();
	void InitOtherVariables();
	void InitializeFromNums(int n00, int n01, int n10, int n11);
	void InitializeFromVectors(Vector &yIn, Vector &xIn);
	void InitializeFromVecMat(Vector &yIn, Matrix &xIn); // here xIn is a nSample by 1 matrix
	Table2by2(Vector & yIn, Vector & xIn);

	// a tool to increment the count in the 2 by 2 table
	void Increment(int rowIndex, int colIndex);
	void UpdateMarginSum();
	int GetValue( int rowIndex, int colIndex);
	void ChisqTest(bool &valid);
	void PearsonChisq();

	// tests for 2 by 2 tables
	bool MustExactTest();
	void FisherExactTest();
	void GeneralTest();

	// Biostat 815: Implementation of fast fisher's exact test
	void CalculateBoundsIn00ForFisher();
	void initLogFacs(double *logFacs);
	void initLogFacs(double *logFacs, int n);
	double logHypergeometricProb(double *logFacs);
	double logHypergeometricProb(double *logFacs, int a, int b, int c, int d);
	void FullFastFisherExactTest();

	// get p values
	double GetpChisq();
	double GetpExact();
	double GetstatChisq();
	double GetpGeneral();
	double GetpFasterFisherExact(const char* testIndicator); // testIndicator ("TwoSided","Less","Greater")

	// other tools
	void	Print();
	void PrintMargin();
	void GetTable(int b[2][2]);
    int & Get00() {return(m_n[0][0]);};
    int & Get01() {return(m_n[0][1]);};
    int & Get10() {return(m_n[1][0]);};
    int & Get11() {return(m_n[1][1]);};

	double getPExactTwoSided() const {return this->pExactTwoSided; };
	double getPExactOneSidedLess() const {return this->pExactOneSidedLess; };
	double getPExactOneSidedGreater() const {return this->pExactOneSidedGreater; };

private:
	int m_n[2][2];
	int m_marginCol[2];
	int m_marginRow[2];
	int m_sum;
	int lowerBound,upperBound;
	double chisqStat;
	double pChisq;
	double pExact;
	double pGeneral;
	double pExactTwoSided;
	double pExactOneSidedLess;
	double pExactOneSidedGreater;
};

#endif
