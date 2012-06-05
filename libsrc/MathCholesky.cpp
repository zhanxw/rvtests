////////////////////////////////////////////////////////////////////// 
// libsrc/MathCholesky.cpp 
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
 
#include "MathCholesky.h"
#include "Error.h"

#include <math.h>

void Cholesky::Decompose(Matrix & A)
   {
   L.Dimension(A.rows, A.rows);
   L.Zero();
   FastDecompose(A);
   }

void Cholesky::FastDecompose(Matrix & A)
   {
   if (A.rows != A.cols)
      error("Cholesky.Decompose: Matrix %s is not square",
            (const char *) A.label);

   L.Dimension(A.rows, A.rows);

   for (int i=0; i<L.rows; i++)
      for (int j=i; j<L.rows; j++)
         {
         double sum = A[i][j];
         for (int k = i - 1; k >= 0; k--)
            sum -= L[i][k] * L[j][k];
         if (i == j)
            if (sum <= 0.0)
               error("Cholesky - matrix %s is not positive definite",
                     (const char *) A.label);
            else
               L[i][i] = sqrt(sum);
         else
            L[j][i] = sum / L[i][i];
         }
   }

bool Cholesky::TryDecompose(Matrix & A)
   {
   L.Dimension(A.rows, A.rows);
   L.Zero();

   if (A.rows != A.cols)
     {
       printf("The matrix to be decomposed is not square!");
       return false;
     }


   L.Dimension(A.rows, A.rows);

   for (int i=0; i<L.rows; i++)
      for (int j=i; j<L.rows; j++)
         {
         double sum = A[i][j];
         for (int k = i - 1; k >= 0; k--)
            sum -= L[i][k] * L[j][k];
         if (i == j)
            if (sum <= 0.0)
	      {
#ifndef NDEBUG            
            printf("%s:%d sum < = 0.0\n", __FILE__, __LINE__);
#endif
		return false;
	      }
            else
               L[i][i] = sqrt(sum);
         else
            L[j][i] = sum / L[i][i];
         }

   return true;
   }

void Cholesky::BackSubst(Vector & b)
   {
   x.Dimension(L.rows);

   // Solve L*v = b (store v in x)
   for (int i = 0; i < L.rows; i++)
      {
      double sum = b[i];
      for (int k = i-1; k>=0; k--)
         sum -= L[i][k] * x[k];
      x[i] = sum / L[i][i];
      }

   // Solve transpose(L)*x = v
   // End result is ... A*x = L*t(L)*x = L*v = b
   for (int i=L.rows-1; i>=0; i--)
      {
      double sum = x[i];
      for (int k = i+1; k < L.rows; k++)
         sum -= L[k][i] * x[k];
      x[i] = sum / L[i][i];
      }

   // Done!
   }

void Cholesky::Invert()
   {
   inv.Dimension(L.rows, L.rows);

   inv.Identity();

   for(int i = 0; i < L.rows; i++)
      {
      BackSubst(inv[i]);
      inv[i] = x;
      }
   }

double Cholesky::lnDeterminantL()
   {
   double sum = 0;
   for (int i = 0; i < L.rows; i++)
      sum += log(L[i][i]);
   return sum;
   }

double Cholesky::DeterminantL()
   {
   double product = 1;
   for (int i=0; i<L.rows; i++)
      product *= L[i][i];
   return product;
   }

 
