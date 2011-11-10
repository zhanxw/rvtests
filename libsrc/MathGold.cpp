////////////////////////////////////////////////////////////////////// 
// libsrc/MathGold.cpp 
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
 
#include "MathGold.h"
#include "Error.h"

#include <math.h>

#define GLIMIT             100      // Maximum magnification, per round
#define SHIFT(a, b, c, d)  (a)=(b); (b)=(c); (c)=(d)
#define SHFT3(a, b, c)     (a)=(b); (b)=(c)

void ScalarMinimizer::Bracket(double lo, double hi)
// Try and bracket the minimum of function f(x), using
// a and b as initial guesstimates.
   {
   a = lo;
   b = hi;
   fa = f(a);
   fb = f(b);

   if (fb > fa)
      {
      double temp;
      SHIFT(temp, a, b, temp);
      SHIFT(temp, fa, fb, temp);
      }

   c = b + (GOLD + 1.0) * (b - a);
   fc = f(c);

   while (fb > fc)      // when we stop going down, the minimum is bracketed
      {                 // parabolic extrapolation
      double r = (b - a) * (fb - fc);
      double q = (b - c) * (fb - fa);
      double u =  b - ((b - c) * q - (b - a) * r) /
                      (2.0 * sign(max(fabs(q - r), TINY), q - r));

      double ulim = b + GLIMIT*(c - b);
      double fu;

      if ((b - u) * (u - c) > 0.0)  // Parabolic between u and b?
         {
         fu = f(u);
         if (fu < fc)               // Minimum between b and c
            {
            SHFT3(a, b, u);
            SHFT3(fa, fb, fu);
            return;
            }
         else if (fu > fb)          // Minimum between a and u
            {
            c = u;
            fc = fu;
            return;
            }
         u = c + (GOLD + 1.0) * (c - b);    // Magnify!
         fu = f(u);
         }
      else if ((c - u) * (u - ulim) > 0.0)   // Parabolic between c and ulim
         {
         fu = f(u);
         if (fu < fc)
            {
            SHIFT(b, c, u, c + (GOLD + 1.0) * (c - b));
            SHIFT(fb, fc, fu, f(u));
            }
         }
      else if ((u - ulim)*(ulim - c) >= 0.0) // Constrained by ulim
         {
         u = ulim;
         fu = f(u);
         }
      else                                    // Magnify!
         {
         u = c + (GOLD + 1.0)*(c - b);
         fu = f(u);
         }
      SHIFT(a, b, c, u);
      SHIFT(fa, fb, fc, fu);
      }
   }

double ScalarMinimizer::Brent(double tol)
   {
   double temp;

   if (a > c)
      {
      SHIFT(temp, a, c, temp);
      SHIFT(temp, fa, fc, temp);
      }

   min = b; fmin = fb;
   double w = b, v = b;
   double fw = fb, fv = fb;

   double delta = 0.0;         // distance moved in step before last

   double u, fu, d = 0.0;      // Initializing d is not necesary, but avoids warnings
   for (int iter = 1; iter <= ITMAX; iter++)
      {
      double middle = 0.5 * (a + c);
      double tol1 = tol * fabs(min) + ZEPS;
      double tol2 = 2.0 * tol1;

      if (fabs(min - middle) <= (tol2 - 0.5 * (c - a)))
         return fmin;

      if (fabs(delta) > tol1)
         // Try a parabolic fit
         {
         double r = (min - w) * (fmin - fv);
         double q = (min - v) * (fmin - fw);
         double p = (min - v) * q - (min - w) * r;

         q = 2.0 * (q - r);
         if (q > 0.0) p = -p;
         q = fabs(q);

         temp = delta;
         delta = d;

         if (fabs(p) >= fabs(0.5 * q * temp) ||
             p <= q * (a - min) || p >= q * (c - min))
            // parabolic doesn't look like a good idea
            {
            delta = min >= middle ? a - min : c - min;
            d = CGOLD * delta;
            }
         else
            // parabolic fit is the way to go
            {
            d = p / q;
            u = min + d;
            if (u - a < tol2 || c - u < tol2)
               d = sign(tol1, middle - min);
            }
         }
      else
         {
         // Golden ratio for first step
         delta = min >= middle ? a - min : c - min;
         d = CGOLD * delta;
         }

      // Don't take steps smaller than tol1
      u = fabs(d) >= tol1 ? min + d : min + sign(tol1, d);
      fu = f(u);

      if (fu <= fmin)
         {
         if (u >= min)
            a = min;
         else
            c = min;
         SHIFT(v, w, min, u);
         SHIFT(fv, fw, fmin, fu);
         }
      else
         {
         if (u < min)
            a = u;
         else
            c = u;
         if (fu <= fw || w == min)
            {
            SHFT3(v, w, u);
            SHFT3(fv, fw, fu);
            }
         else if (fu <= fv || v == min || v == w)
            {
            v = u;
            fv = fu;
            }
         }
      }
      numerror("ScalarMinimizer::Brent got stuck");
      return fmin;
   }

double ScalarMinimizer::f(double x)
   {
   return func(x);
   };

// Scalar minimizer class
//

LineMinimizer::LineMinimizer()
   : ScalarMinimizer(), line(), point(), temp()
   {
   garbage = false;
   func = NULL;
   }

LineMinimizer::LineMinimizer(double (*myfunc)(Vector & v))
   : ScalarMinimizer(), line(), point(), temp()
   {
   garbage = true;
   func = new VectorFunc(myfunc);
   }

LineMinimizer::LineMinimizer(VectorFunc & vfunc)
   : ScalarMinimizer(), line(), point(), temp()
   {
   garbage = false;
   func = &vfunc;
   }

double LineMinimizer::f(double x)
   {
   temp = point;
   temp.AddMultiple(x, line);
   return func->Evaluate(temp);
   };

double LineMinimizer::Brent(double tol)
   {
   ScalarMinimizer::Brent(tol);

   line.Multiply(min);
   point.Add(line);

   return fmin;
   };


 
