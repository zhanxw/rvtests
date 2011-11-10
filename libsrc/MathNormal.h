////////////////////////////////////////////////////////////////////// 
// libsrc/MathNormal.h 
// (c) 2000-2010 Goncalo Abecasis
// 
// This file is distributed as part of the Goncalo source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile Goncalo.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Sunday May 02, 2010
// 
 
#ifndef __NORMALEQUATIONS_H__
#define __NORMALEQUATIONS_H__

#include "IntArray.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "MathCholesky.h"

#define NORMAL_AMOEBA_MIN    0
#define NORMAL_POWELL_MIN    1
#define NORMAL_FLETCHER_MIN  2

class NormalEquations
   {
   public:
      Vector   means, variances;

      Matrix * varComponents;
      Matrix   linearModel;
      Vector   scores;

      double   likelihood;

      NormalEquations();
      virtual ~NormalEquations();

      void   Dimension(int vcCount);

      virtual void   Prepare();
      virtual void   SetParameters(Vector & means, Vector & variances);
      virtual double Evaluate();

      Cholesky  cholesky;
      Matrix    varMatrix;
      Vector    residuals;
      double    constant;
      bool      includeLikelihoodConstant;
      int       multiple;

      bool operator == (const NormalEquations & rhs);

      void      EnableConstant();
      void      DisableConstant();

      // Diagnostic statistics
      // see JL Hopper and JD Matthews Ann Hum Genet (1992) 46:373 - 383
      double rawQ;      // This is a chi-square with n degrees of freedom
      double Q;         // This is Q1 and has a standard normal distribution
      Vector Qi;        // Each Qi is approximately chi-square with 1 df
      void   Diagnostics();

   protected:
      void Free();

      void CalculateResiduals();
      void CalculateCovariances();

      bool      meanChange, varChange, init;
      IntArray  meanFlags;
   };

class NormalSet
   {
   public:
      NormalEquations ** sets;
      Vector    weights;
      IntArray  operators;

      double precision;
      int    numericMinimizer;
      int    size;
      int    count;
      int    maxThreads;
      double likelihood;
      Vector variances, means;

      // Number of function evaluations
      int    evaluations;

      NormalSet(int threads = 0);

      virtual ~NormalSet() { Free(); }

      void        Dimension(int setCount, int vcCount, int vcDerived = 0);
      double      Evaluate();
      void        SelectPoint(Vector & v);
      void        Solve();
      int         CountObservations();
      virtual int CountParameters();

      NormalEquations & operator [] (int n)
         { return *(sets[n]); }

      void EnableConstant();
      void DisableConstant();

      // Vector for storing intermediate likelihoods
      Vector recordedLikelihoods;

      // This function should be over-ridden to calculate constrained variance
      // components appropriately
      virtual void CalculateConstrainedVariances();

   protected:
      // for multi-threading
      static void * EvaluateOneSet(void * which);

      // house-keeping
      void  Free();
      virtual void AllocateSets();

      // Helpers for solver
      void         EditLinearDegenerates();
      virtual void GetStartingPoint(Vector & startPoint);
      void         RemoveRedundancy();

      // Intermediate results when calculating likelihoods
      Vector logLikelihoods;

      // How many variance components should be estimated?
      int  vcEstimated;

      // And how many variance components are constrained by the other parameters?
      int  vcConstrained;

      // relative size of variance components
      double varScale;
   };

class NonLinearNormalSet : public NormalSet
   {
   public:
      virtual void CalculateConstrainedVariances();

      IntArray  nonLinearVariances;
      IntArray  component1, component2;
   };

class NormalSolver : public VectorFunc
   {
   public:
      NormalSet * normal;

      NormalSolver(NormalSet * n) : VectorFunc()
         { normal = n; }

      virtual ~NormalSolver() { }

      virtual double Evaluate(Vector & point);
   };

// Constants for setting elements of operations array
// Which tell normal set class how to combine partial likelihoods
//

#define NORMAL_OP_MASK    7
#define NORMAL_NOP        0
#define NORMAL_MUL_LK     1
#define NORMAL_SCALE_LLK  2
#define NORMAL_SUM_LK     3
#define NORMAL_DIV_LK     4
#define NORMAL_POP        5
#define NORMAL_RECORD_LLK 6

#define NORMAL_OP(a,b,c,d)  ((a) | ((b) << 3) | ((c) << 6) | ((d) << 9))
#define NORMAL_LAST_OP(a)   ((a) << 12)

#endif


 
