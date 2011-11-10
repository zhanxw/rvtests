////////////////////////////////////////////////////////////////////// 
// libsrc/MathNormal.cpp 
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
 
#include "MathNormal.h"
#include "MathSVD.h"
#include "MathGenMin.h"
#include "MathStats.h"

#include <math.h>

NormalEquations::NormalEquations() :
   means("linear"), variances("variances"),
   linearModel("linear"), scores("scores"),
   cholesky(), varMatrix("variances"), residuals("residuals"),
   Qi("Qi")
   {
   variances.Dimension(0);
   varComponents = NULL;
   includeLikelihoodConstant = false;
   multiple = 1;
   }

NormalEquations::~NormalEquations()
   {
   Free();
   }

void NormalEquations::Free()
   {
   if (variances.Length())
      delete [] varComponents;
   varComponents = NULL;
   variances.Dimension(0);
   }

void NormalEquations::CalculateResiduals()
   {
   // deviation from predicted scores in linear model
   residuals.Product(linearModel, means);
   residuals.Subtract(scores);
   }

void NormalEquations::CalculateCovariances()
   {
   // Construct the covariance matrix by summing its
   // component parts and then multiply the residuals
   // by inverse of covariance matrix
   //    - Inv(A)*b = x where A*x = b and x = BackSubst(b)
   //    - cholesky decomposition only requires the upper
   //      triangle of the positive definite A
   varMatrix.Zero();
   for (int j = 0; j < variances.dim; j++)
      // Ignore empty 0x0 matrices
      if (varComponents[j].rows && varComponents[j].cols)
         for (int r = 0; r < varMatrix.rows; r++)
            for (int c = r; c < varMatrix.rows; c++)
               varMatrix[r][c] += variances[j] * varComponents[j][r][c];

   /* The following two lines achieve the same in a
      prettier but more inefficient manner:

      for (int j = 0; j < variances.dim; j++)
         varMatrix.AddMultiple(variances[j], varComponents[j]);
   */
   }

double NormalEquations::Evaluate()
   {
   // If nothing changed ...
   if (!meanChange && !varChange)
      return likelihood;

   // (2pi)^(-k/2)
   likelihood = includeLikelihoodConstant ? constant : 0.0;

   // deviation from predicted scores in linear model
   if (meanChange) CalculateResiduals();

   // The variance matrix A (varMatrix)
   if (varChange)  CalculateCovariances();

   // and its cholesky decomposition are updated
   if (varChange)  cholesky.FastDecompose(varMatrix);

   // This is cholesky.x = Inv(A) * residuals
   cholesky.BackSubst(residuals);

   // Add e ^ [- 0.5 * (yi - ui)T * A^-1 * (yi - ui)]
   likelihood -= 0.5 * residuals.InnerProduct(cholesky.x);

   // Det(A[i]^-1/2) = -Det(sqrt(A)) = -Det(L)
   likelihood -= cholesky.lnDeterminantL();

   // Flag that allows us to check for redundant calculations
   init = true;

   return likelihood;
   }

void NormalEquations::Diagnostics()
   {
   CalculateResiduals();
   CalculateCovariances();

   cholesky.BackSubst(residuals);

   rawQ = residuals.InnerProduct(cholesky.x);
   Q = sqrt(2.0 * rawQ) - sqrt(2.0 * residuals.dim - 1.0);

   Qi.Dimension(residuals.dim);

   Matrix M;
   Vector v, r;
   Cholesky chol;

   for (int i = 0; i < residuals.dim; i++)
      {
      M = varMatrix;

      v = M[i];
      M.DeleteColumn(i);
      M.DeleteRow(i);
      v.DeleteDimension(i);
      chol.Decompose(M);

      chol.BackSubst(v);

      double var = varMatrix[i][i] - v.InnerProduct(chol.x);

      r = residuals;
      r.DeleteDimension(i);

      chol.BackSubst(r);
      Qi[i] = residuals[i] - v.InnerProduct(chol.x);
      Qi[i] = Qi[i] * Qi[i] / var;
      }
   }

void NormalEquations::Prepare()
   {
   varMatrix.Dimension(scores.dim, scores.dim);
   residuals.Dimension(scores.dim);
   means.Dimension(linearModel.cols);

   meanFlags.Dimension(linearModel.cols);
   meanFlags.Zero();
   for (int i = 0; i < linearModel.cols; i++)
      for (int j = 0; j < linearModel.rows; j++)
         if (linearModel[j][i] != 0.0)
            {
            meanFlags[i] = 1;
            break;
            }

   constant = log(2 * M_PI) * -0.5 * scores.dim;
   init = false;
   }

void NormalEquations::SetParameters(Vector & mu, Vector & sigma)
   {
   if (!init)
      {
      means = mu;
      variances = sigma;
      varChange = meanChange = true;
      return;
      }

   varChange = meanChange = false;

   for (int i = 0; i < mu.Length(); i++)
      if (mu[i] != means[i] && meanFlags[i])
         {
         means[i] = mu[i];
         meanChange = true;
         }

   for (int i = 0; i < sigma.Length(); i++)
      if (sigma[i] != variances[i] && varComponents[i].rows)
         {
         variances[i] = sigma[i];
         varChange = true;
         }
   }

void NormalEquations::Dimension(int vcCount)
   {
   if (variances.Length() != vcCount)
      {
      int swapCount = vcCount > variances.Length() ? variances.Length() : vcCount;

      Matrix * newVarComponents = new Matrix[vcCount];
      for (int i = 0; i < swapCount; i++)
         newVarComponents[i].Swap(varComponents[i]);
      if (varComponents != NULL)
         delete [] varComponents;

      variances.Dimension(vcCount);
      varComponents = newVarComponents;
      }
   }

bool NormalEquations::operator == (const NormalEquations & rhs)
   {
   if (scores != rhs.scores) return false;
   if (linearModel != rhs.linearModel) return false;
   for (int i = 0; i < variances.dim; i++)
      if (varComponents[i] != rhs.varComponents[i])
         return false;
   return true;
   }

void NormalEquations::EnableConstant()
   {
   if (!includeLikelihoodConstant)
      {
      includeLikelihoodConstant = true;

      if (init) likelihood += constant;
      }
   }

void NormalEquations::DisableConstant()
   {
   if (includeLikelihoodConstant)
      {
      includeLikelihoodConstant = false;

      if (init) likelihood -= constant;
      }
   }

// Normal Set class
//

NormalSet::NormalSet(int threads)
   : variances("variances"), means("means")
   {
   maxThreads = threads;
   count = size = 0;
   numericMinimizer = 2;
   precision = 1e-8;
   varScale = 1.0;
   }

void NormalSet::Free()
   {
   if (size)
      {
      for (int i = 0; i < size; i++)
         delete sets[i];
      delete [] sets;
      }
   size = 0;
   }

void NormalSet::Dimension(int setCount, int vcCount, int vcDerived)
   {
   weights.Dimension(setCount);
   operators.Dimension(setCount);

   for (int i = size; i < setCount; i++)
      {
      weights[i] = 1.0;
      operators[i] = NORMAL_MUL_LK;
      }

   if (setCount > size)
      {
      Free();

      size = setCount;
      sets = new NormalEquations * [size];

      AllocateSets();
      }

   count = setCount;
   for (int i = 0; i < count; i++)
      sets[i]->Dimension(vcCount);

   vcEstimated = vcCount - vcDerived;
   vcConstrained = vcDerived;
   }

void NormalSet::AllocateSets()
   {
   for (int i = 0; i < size; i++)
      sets[i] = new NormalEquations;
   }

double NormalSet::Evaluate()
   {
   evaluations++;

   for (int i = 0; i < count; i++)
        sets[i]->Evaluate();

   logLikelihoods.Clear();
   logLikelihoods.Push(0.0);
   recordedLikelihoods.Clear();

   for (int i = 0; i < count; i++)
      {
      logLikelihoods.Push(sets[i]->likelihood);

      int op = operators[i];

      do
         switch (op & NORMAL_OP_MASK)
            {
            case NORMAL_SUM_LK:
               {
               double lk1 = logLikelihoods.Pop();
               double lk2 = logLikelihoods.Pop();

               if (lk1 > lk2)
                 logLikelihoods.Push(lk1 + log(1.0 + exp(lk2 - lk1)));
               else
                 logLikelihoods.Push(lk2 + log(1.0 + exp(lk1 - lk2)));
               }
               break;
            case NORMAL_MUL_LK:
               {
               double lk = logLikelihoods.Pop();
               logLikelihoods.Last() += lk;
               }
               break;
            case NORMAL_DIV_LK:
               {
               double lk = logLikelihoods.Pop();
               logLikelihoods.Last() -= lk;
               }
               break;
            case NORMAL_SCALE_LLK:
               logLikelihoods.Last() += weights[i];
               break;
            case NORMAL_POP:
               logLikelihoods.Pop();
               break;
            case NORMAL_RECORD_LLK:
               recordedLikelihoods.Push(logLikelihoods.Last());
               break;
            }
      while ((op >>= 3) != 0);
      }

   if (logLikelihoods.Length() != 1)
      error("Internal error in likelihood formulation");

   likelihood = -logLikelihoods[0];

   return likelihood;
   }

void * NormalSet::EvaluateOneSet(void * which)
   {
   NormalEquations * w = (NormalEquations *) which;

   w->Evaluate();

   return NULL;
   }

void NormalSet::SelectPoint(Vector & v)
   {
   for (int i = 0; i < means.dim; i++)
      means[i] = v[i];

   for (int i = 0, j = means.dim; i < vcEstimated; i++, j++)
      if (v[j] > 16 || v[j] < -16)
         variances[i] = (v[j] > 0 ? 1e7 : 1e-7) * varScale;
      else
         variances[i] = exp(v[j]) * varScale;

   if (vcConstrained) CalculateConstrainedVariances();

   for (int i = 0; i < count; i++)
      sets[i]->SetParameters(means, variances);
   }

void NormalSet::Solve()
   {
   EditLinearDegenerates();

   GeneralMinimizer * solver = NULL;   // Initialization avoids compiler warnings

   switch (numericMinimizer)
      {
      case NORMAL_AMOEBA_MIN : solver = new AmoebaMinimizer(); break;
      case NORMAL_POWELL_MIN : solver = new PowellMinimizer(); break;
      case NORMAL_FLETCHER_MIN : solver = new FletcherMinimizer(); break;
      }

   solver->func = new NormalSolver(this);

   // set the number of parameters to minimize
   int parameters = CountParameters();
   solver->Reset(parameters);

   // Reset the number of likelihood evaluations
   evaluations = 0;

   // If we are not using the Nelder-Mead minimizer, use it
   // to conduct a rough pre-optimization.
   if (solver != NORMAL_AMOEBA_MIN)
      {
      AmoebaMinimizer presolver;

      presolver.func = solver->func;
      presolver.Reset(parameters);
      GetStartingPoint(presolver.point);
      presolver.Minimize(precision * 1000);

      solver->point = presolver.point;
      }
   else
      GetStartingPoint(solver->point);

   // Two rounds of minimizing to be safe
   solver->Minimize(precision);

   double lastmin;
   double scale = 2.0;

   do {
      // Find the largest variance ...
      double varMax = solver->point[means.Length()];
      for (int i = means.Length() + vcEstimated - 1; i > means.Length(); i--)
         if (solver->point[i] > varMax)
            varMax = solver->point[i];

      // Check that none of the variances is effectively zero...
      for (int i = means.Length() + vcEstimated - 1; i >= means.Length(); i--)
         if (solver->point[i] < (varMax - 5.0))
            solver->point[i] = varMax - 5.0;

      lastmin = solver->fmin;
      scale *= -0.5;
      solver->Reset(parameters, scale);
      solver->Minimize(precision);
      }
   while (solver->fmin > precision &&
         (lastmin - solver->fmin)/solver->fmin > precision);

   SelectPoint(solver->point);

   delete solver->func;
   delete solver;
   }

int NormalSet::CountObservations()
   {
   int rows = 0;
   for (int i = 0; i < count; i++)
      rows += sets[i]->linearModel.rows;

   return rows;
   }

void NormalSet::EditLinearDegenerates()
   {
   // create a matrix M with the linear side of the model
   // for each of the normal equations
   int    rows = CountObservations();
   Matrix m(rows, sets[0]->linearModel.cols);;
   Vector b(rows);

   for (int i = 0, total = 0; i < count; i++)
      for (int j = 0; j < sets[i]->linearModel.rows; j++)
         {
         b[total] = sets[i]->scores[j];
         m[total++] = sets[i]->linearModel[j];
         }

   // now we decompose m successively and remove redundant
   // columns from m and each individual set...
   SVD engine;

   for (int col = 1; col <= m.cols; )
      {
      engine.Decompose(m, m.rows, col);
      engine.Edit();

      bool redundant = false;
      for (int j = 0; j < col; j++)
         if (engine.w[j] == 0.0)
            redundant = true;

      if (redundant)
         {
         for (int i = 0; i < count; i++)
            sets[i]->linearModel.DeleteColumn(col - 1);
         m.DeleteColumn(col - 1);

         // Check if last column was deleted ...
         if (col > m.cols) engine.Decompose(m, m.rows, m.cols);
         }
      else
         col++;
      }

  // at this point we might as well use the SVD information
  // to choose good starting points
  for (int i = 0; i < count; i++)
     sets[i]->Prepare();

  variances.Dimension(sets[0]->variances.dim);

  engine.BackSubst(b);

  // Assume all variances * are * the total variance
  // and that SVD estimates are good for the linear side

  means = engine.x;
  means.Dimension(m.cols);

  // make sure the variances are positive despite roundoff
  double _TINY_ = 1e-20;
  double residualVar = fabs(engine.RSS(m, b)) + _TINY_;
  residualVar /= rows;

  // Calculate weights for each variance, so that starting estimates
  // assume all observations have variance approximately equal to estimate
  // for pooled sample.
  Vector   var_weights(variances.Length()); var_weights.Zero();
  IntArray var_counts(variances.Length());  var_counts.Zero();
  Vector   scaled_variances(variances.Length());

  for (int i = 0; i < count; i++)
    for (int j = 0; j < sets[i]->scores.Length(); j++)
      {
      for (int v = 0; v < variances.Length(); v++)
        {
        if (sets[i]->varComponents[v].rows)
          scaled_variances[v] = sets[i]->varComponents[v][j][j];
        else
          scaled_variances[v] = 0.0;
        }
      scaled_variances *= 1.0 / scaled_variances.Sum();

      for (int v = 0; v < variances.Length(); v++)
         if (scaled_variances[v] > 0.0)
            {
            var_weights[v] += scaled_variances[v];
            var_counts[v] ++;
            }
      }

  for (int v = 0; v < vcEstimated; v++)
    variances[v] = residualVar * var_weights[v] / var_counts[v];

  varScale = residualVar;
  }

int NormalSet::CountParameters()
   {
   return means.dim + vcEstimated;
   };

void NormalSet::GetStartingPoint(Vector & startPoint)
   {
   for (int i = 0; i < means.dim; i++)
      startPoint[i] = means[i];
   for (int i = 0; i < variances.dim; i++)
      startPoint[i + means.dim] = log(variances[i] / varScale);
   }

void NormalSet::RemoveRedundancy()
   {
   for (int i = count - 1; i >= 0; i--)
      for (int j = 0; j < i; j++)
         if (*(sets[i]) == *(sets[j]))
            {
            sets[j]->multiple++;
            count--;
            NormalEquations * temp = sets[i];
            for ( j = i; j < count; j++)
               sets[j] = sets[j + 1];
            sets[j] = temp;
            }
   }

void NormalSet::EnableConstant()
   {
   for (int i = 0; i < count; i++)
      sets[i]->EnableConstant();
   }

void NormalSet::DisableConstant()
   {
   for (int i = 0; i < count; i++)
      sets[i]->DisableConstant();
   }

void NormalSet::CalculateConstrainedVariances()
   {
   error("This function must be over-ridden");
   }

// NonLinearNormalSet
//

void NonLinearNormalSet::CalculateConstrainedVariances()
   {
   int vc = variances.Length();

   variances[vc - 1] = sqrt(variances[vc - 2] * variances[vc - 3]);

/*
   for (int i = 0; i < nonLinearVariances.Length(); i++)
       variances[i] = sqrt(variances[component1[i]] * variances[component2[i]]);
*/
   }

// Normal Solver classes
//

double NormalSolver::Evaluate(Vector & point)
   {
   normal->SelectPoint(point);
   return normal->Evaluate();
   };

 
