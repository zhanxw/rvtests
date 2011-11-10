////////////////////////////////////////////////////////////////////// 
// libsrc/MathSVD.cpp 
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
 
#include "MathSVD.h"

#include <math.h>

SVD::SVD() : u("svd.u"), w("svd.w"), v("svd.v"), x("svd.x")
   {
   m = n = 0;
   }

SVD::~SVD()
   {
   Empty();
   }

void SVD::Empty()
   {
   m = n = 0;
   }

void SVD::Decompose(Matrix & a, int mp, int np)
   {
   int      flag, i, its, j, jj, k, l, nm;
   double   c, f, g, h, s, scale, x, y, z;
   double   anorm;
   Vector   rv1;

   Empty();

#ifndef __BORLANDC__
   l = nm = 0; // Initialization avoids compiler warnings (except for BORLAND C!)
#endif

   m = (mp == -1) ? a.rows : mp;
   n = (np == -1) ? a.cols : np;

//   a.Print(stdout,m,n);

   u.Copy(a);
   u.Dimension(m, n);
   v.Dimension(n, n);
   w.Dimension(n);

   rv1.Dimension(n);

   g = scale = anorm = 0.0;      // Householder reduction to
                                 // bidiagonal form
   for (i=0; i<n; i++)
      {
      l = i + 1;
      rv1[i] = scale * g;
      g = s = scale = 0.0;
      if (i < m)
         {
         for (k=i; k<m; k++) scale += fabs(u[k][i]);
         if (scale > FPMIN)
            {
            for (k=i; k<m; k++)
               {
               u[k][i] /= scale;
               s += u[k][i]*u[k][i];
               }
            f = u[i][i];
            // g = -SIGN(sqrt(s), f)
            g = f < 0.0 ? sqrt(s) : -sqrt(s);
            h = f * g - s;
            u[i][i] = f - g;
            for (j=l; j<n; j++)
               {
               for ( s=0.0, k=i; k<m; k++)
                  s += u[k][i]*u[k][j];
               f = s/h;
               for ( k=i; k<m; k++)
                  u[k][j] += f*u[k][i];
               }
            for (k=i; k<m; k++) u[k][i] *= scale;
            }
         }
      w[i] = scale * g;
      g = s = scale = 0.0;
      if ((i < m) && (i != n - 1))
         {
         for ( k=l; k<n; k++)
            scale += fabs(u[i][k]);
         if (scale > FPMIN)
            {
            for (k=l; k<n; k++)
               {
               u[i][k] /= scale;
               s += u[i][k]*u[i][k];
               }
            f = u[i][l];
            // g = -SIGN(sqrt(s), f)
            g = f < 0.0 ? sqrt(s) : -sqrt(s);
            h = f * g - s;
            u[i][l] = f - g;
            for (k=l; k<n; k++)
               rv1[k]=u[i][k]/h;
            for (j=l; j<m; j++)
               {
               for (s=0.0, k=l; k<n; k++)
                  s+= u[j][k]*u[i][k];
               for (k=l; k<n; k++)
                  u[j][k] += s * rv1[k];
               }
            for (k=l; k<n; k++)
               u[i][k] *= scale;
            }
         }
      x = fabs(w[i]) + fabs(rv1[i]);
      if (anorm < x)
          anorm = x;
      }

   // accumulation of right-hand transformations
   for (i=n-1; i>=0; i--)
      {
      if (i < n-1)
         {
         if (fabs(g) > FPMIN)
            {
            // double division to avoid possible underflow
            for (j=l; j<n; j++)
               v[j][i]=(u[i][j]/u[i][l])/g;
            for (j=l; j<n; j++)
               {
               for (s=0.0,k=l; k<n; k++)
                  s += u[i][k]*v[k][j];
               for (k=l; k<n; k++)
                  v[k][j] += s * v[k][i];
               }
            }
         for (j=l; j<n; j++) v[i][j]=v[j][i]=0.0;
         }
         v[i][i]=1.0;
         g=rv1[i];
         l=i;
      }

   // accumulation of left-hand transformations
   for (i=min(m,n)-1; i>=0; i--)
      {
      l = i + 1;
      g = w[i];
      for (j=l; j<n; j++) u[i][j] = 0.0;
      if (fabs(g) > FPMIN)
         {
         g = 1.0 / g;
         for (j=l; j<n; j++)
            {
            for (s = 0.0, k=l; k<m; k++)
               s += u[k][i] * u[k][j];
            f = (s/u[i][i])*g;
            for (k=i;k<m;k++) u[k][j] += f*u[k][i];
            }
         for (j=i;j<m;j++) u[j][i] *= g;
         }
      else
         for (j=i;j<m;j++) u[j][i] = 0.0;
      ++u[i][i];
      }

   // Diagonalization of the bi-diagonal form:
   // Loop over singular values and over allowed iterations
   for (k=n-1;k>=0;k--)
      {
      for (its=1; its<=30; its++)
         {
         flag=1;
         for (l=k; l>=0; l--)
            {
            // Test for splitting
            // note that rv1[1] is always zero.
            nm = l - 1;
            if (check_equality(fabs(rv1[l])+anorm, anorm))
               {
               flag = 0;
               break;
               }
            if (check_equality(fabs(w[nm])+anorm, anorm))
               break;
            }

//         for (int index = l; index <= k; index++)
//            printf("rv1[%d] = %g\n", index, rv1[index]);
//         printf("rv1[l] = %f, anorm = %f, l = %d, flag = %d\n", rv1[l], anorm, l, flag);

         if (flag)
            {
            c = 0.0;       // cancellation of rv1[l], if l > 1
            s = 1.0;
            for (i=l; i<k; i++)
               {
               f = s*rv1[i];
               rv1[i] = c*rv1[i];
               if (check_equality(fabs(f)+anorm,anorm))
                  break;
               g = w[i];
               h = pythag(f,g);
               w[i] = h;
               h = 1.0/h;
               c = g*h;
               s = -f*h;
               for (j=0; j<m; j++)
                  {
                  y = u[j][nm];
                  z = u[j][i];
                  u[j][nm] = y*c + z*s;
                  u[j][i] = z*c - y*s;
                  }
               }
            }
         z = w[k];
         if (l==k)         // Convergence
            {
            if (z < 0.0)   // Singular value is made nonnegative
               {
               w[k] = -z;
               for (j=0; j<n; j++)
                  v[j][k] = -v[j][k];
               }
            break;
            }
         if (its == 30)
            {
            FILE * f = fopen("SVD-failure.txt", "wt");

            if (f != NULL)
               {
               a.Print(f, m, n);
               fclose(f);

               error("No convergence in 30 SVD Decomposition iterations\n\n"
                     "Problematic matrix could help debug the problem,\n"
                     "saved to file [SVD-failure.txt]\n");
               }

            error("No convergence in 30 SVD Decomposition iterations");
            }

         // shift from bottom 2-by-2 minor
         x = w[l];
         nm = k - 1;
         y = w[nm];
         g = rv1[nm];
         h = rv1[k];
         f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
         g = pythag(f, 1.0);
         f = (
             (x-z)*(x+z) +
             // h * ((y/(f+SIGN(g,f)))-h)
             h*((y/(f+(f<0?-g:g)))-h)
             )/x;
         c = s = 1.0;
         // next QR transformation
         for (j=l;j<=nm;j++)
            {
            i = j+1;
            g = rv1[i];
            y = w[i];
            h = s*g;
            g = c*g;
            z = pythag(f,h);
            rv1[j] = z;
            c = f/z;
            s = h/z;
            f = x*c + g*s;
            g = g*c - x*s;
            h = y*s;
            y *= c;
            for (jj=0; jj<n; jj++)
               {
               x = v[jj][j];
               z = v[jj][i];
               v[jj][j] = x*c + z*s;
               v[jj][i] = z*c - x*s;
               }
            z = pythag(f,h);
            w[j] = z;         // Rotation can be arbitrary if z=0
            if (fabs(z) > FPMIN)
               {
               z = 1.0 / z;
               c = f * z;
               s = h * z;
               }
            f = c*g + s*y;
            x = c*y - s*g;
            for (jj=0; jj<m; jj++)
               {
               y = u[jj][j];
               z = u[jj][i];
               u[jj][j] = y*c + z*s;
               u[jj][i] = z*c - y*s;
               }
            }
         rv1[l]=0.0;
         rv1[k]=f;
         w[k]=x;
         }
      }
   }

void SVD::Edit(double tol)
   {
   int j;
   double wmax, wmin;

//   w.Print(stdout);

   wmax = 0.0;

   for (j=0; j<n; j++)
      if (w[j] > wmax)
         wmax = w[j];

   wmin = wmax * tol;

   for (j=0; j<n; j++)
      if (w[j] < wmin)
         w[j] = 0.0;
   }

void SVD::BackSubst(Vector & b)
   {
   int jj, j, i;
   double s;
   Vector tmp;

   x.Dimension(n);
   tmp.Dimension(n);

   // calculate U^T * B
   for (j=0; j<n; j++)
      {
      s = 0.0;
      if (w[j])      // Nonzero result only if wj is nonzero
         {
         for (i=0; i<m; i++)
            s += u[i][j] * b[i];
         s /= w[j];
         }
      tmp[j] = s;
      }

   // Matrix multiply by V to get answer
   for (j=0; j<n; j++)
      {
      s = 0.0;
      for (jj=0; jj<n; jj++)
         s += v[j][jj] * tmp[jj];
      x[j] = s;
      }
   }

double SVD::RSS(Matrix & M, Vector & b)
   {
   double rss = 0.0;
   for (int i=0; i<m; i++)
      {
      double partial = 0.0;
      for (int j=0; j<n; j++)
         partial += x[j] * M[i][j];
      partial -= b[i];
      rss += partial * partial;
      }
   return rss;
   }

void SVD::Residuals(Matrix & M, Vector & b, Vector & delta)
   {
   delta.Dimension(m);

   for (int i=0; i<m; i++)
      {
      double partial = 0.0;
      for (int j=0; j<n; j++)
         partial += x[j] * M[i][j];
      delta[i] = b[i] - partial;
      if (fabs(delta[i]) < ZEPS) delta[i] = 0.0;
      }
   }

void SVD::Covariances()
   {
   cov.Dimension(m, m);

   Vector wti(m);
   wti.Zero();

   for (int i=0; i < m; i++)
      if (w[i]) wti[i] = 1.0 / ( w[i] * w[i]);

   // Sum contributions to covariance matrix
   for (int i=0; i < m; i++)
      for (int j=0; j<= i; j++) {
         double sum = 0.0;
         for ( int k = 0; k < m; k++)
            sum += v[i][k] * v[j][k] * wti[k];
         cov[j][i] = cov[i][j] = sum; }
   }

// calculates SQRT(a*a + b*b) safely
//

double SVD::pythag(double a, double b)
   {
   double absa, absb;
   absa = fabs(a);
   absb = fabs(b);
   if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
      else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa/absb)*(absa/absb)));
   }

// We use a wrapper function to prevent the compiler from using
// higher precision

bool SVD::check_equality(double a, double b)
   {
   return a == b;
   }

// The following function replaces a matrix a with its inverse
//

void SVD::InvertInPlace(Matrix & a)
   {
   Matrix ainv;
   ainv.Dimension(a.cols,a.rows);
   ainv.Zero();
   Decompose(a);
   //Edit();

   Matrix vwinv; // V*W^(-1)
   vwinv.Dimension(n,m);
   vwinv.Zero();

   for (int i = 0; i < n; i++)
     for (int j = 0; j <  m; j++)
      if (w[j] != 0)
         vwinv[i][j] = v[i][j] / w[j];

   // multiply U^t (notice implicit transpose of u);
   for (int i = 0; i < n; i++)
     for (int j = 0; j < m; j++)
        for (int k = 0; k < m; k++)
           ainv[i][j] += vwinv[i][k] * u[j][k];

   a = ainv;
   }

 
