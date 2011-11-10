////////////////////////////////////////////////////////////////////// 
// libsrc/MathGenMin.cpp 
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
 
#include "MathGenMin.h"

#include <math.h>
#include <time.h>

// GeneralMinimizer
//

GeneralMinimizer::GeneralMinimizer() :
   directions(), point()
   {
   func = NULL;
   fmin = FPMAX;      // Very large - not a likely minimum
   };

void GeneralMinimizer::Reset(int n, double scale)
   {
   directions.Dimension(n, n);
   point.Dimension(n);

   directions.Identity();
   directions.Multiply(scale);
   fmin = FPMAX;
   }

// PowellMinimizer
//

double PowellMinimizer::Minimize(double ftol)
   {
   LineMinimizer linmin(*func);
   Vector pstart;

   fmin = f(point);

   iter = 1;
   while (true)
      {
      pstart = point;
      double fstart = fmin;

      double ybig = 0;
      int    ibig = 0;

      // Try all the directions in the current set
      for (int i = 0; i < directions.rows; i++)
         {
         linmin.point = point;
         linmin.line = directions[i];
         linmin.Bracket(0, 1);
         linmin.Brent(ftol);

         if ((fmin - linmin.fmin) > ybig)
            {
            ibig = i;
            ybig = fmin - linmin.fmin;
            }

         fmin = linmin.fmin;
         point = linmin.point;
         }

      // are we done?
      if (2 * (fstart - fmin) <= ftol * (fabs(fstart) + fabs(fmin)))
         return fmin;

      if (iter == ITMAX)
         error("Powell.Minimize - Exceeding maximum iterations");

      // Average direction moved
      linmin.line = point;
      linmin.line.Subtract(pstart);

      // Extrapolated point
      linmin.point.Add(linmin.line);

      double fextra = f(linmin.point);

      if (fextra < fstart)
         {
         double avgdir = fstart - fmin - ybig;
         double newdir = fstart - fextra;
         double promise = 2.0 * (fstart - 2.0 * fmin + fextra) *
                        avgdir * avgdir - ybig * newdir * newdir;

         if (promise < 0.0)
            {
            linmin.Bracket(0,1);
            linmin.Brent(ftol);

            directions.SwapRows(ibig, directions.rows);
            directions[directions.rows - 1] = linmin.line;
            }
         }

      iter++;
      }
   }

// SAMinimizer
//

SAMinimizer::SAMinimizer() :
   GeneralMinimizer(), y(), simplex(), psum(), ptry()
   {
   time_t t;
   rand = new Random((long)time(&t));
   freeRand = true;

   Constructor();
   }

SAMinimizer::SAMinimizer(Random & r) :
   GeneralMinimizer(), y(), simplex(), psum(), ptry()
   {
   freeRand = false;
   rand = &r;

   Constructor();
   }

SAMinimizer::~SAMinimizer()
   {
   if (freeRand) delete rand;
   }

void SAMinimizer::Constructor()
   {
   T = minT = maxT = 0.0;
   iter = Titer = Tcycles = 0;
   }

void SAMinimizer::Reset(int ndim, double scale)
   {
   GeneralMinimizer::Reset(ndim, scale);

   simplex.Dimension(ndim + 1, ndim);
   y.Dimension(ndim + 1);
   psum.Dimension(ndim);
   ptry.Dimension(ndim);

   maxT = 100;
   Tcycles = 20;
   Titer = 1000;
   }

double SAMinimizer::Minimize(double ftol)
   {
   // Setup the simplex, and evaluate f at each vortex
   for (int i=0; i < point.dim; i++)
      {
      simplex[i] = point;
      simplex[i].Add(directions[i]);
      y[i] = f(simplex[i]);
      }

   simplex[simplex.rows - 1] = point;
   y[simplex.rows - 1] = f(simplex[simplex.rows - 1]);

   T = maxT;

   do
      {
      iter = Titer;
      MinimizeLoop(ftol);
      T -= (maxT - minT) / Tcycles;
      }
   while (T > minT);

   return fmin;
   }

double SAMinimizer::MinimizeLoop(double )
   {
   int i, ihi, ilo, m, nvertex = point.dim + 1;
   double ylo, ynhi, ysave, yt, ytry;

   // Calculate psum
   for (psum = simplex[0], m = 1; m < nvertex; m++)
      psum.Add(simplex[m]);

   while (true)
      {
      // Determine which point is highest (worst),
      // next highest, and lowest (best)
      //
      ilo = 0;
      ihi = 1;

      // whenever we look at a vertex it gets a
      // random 'thermal' fluctuation
      ynhi = ylo = y[0] - T * log(rand->Next());
      yhi = y[1] - T * log(rand->Next());
      if (ylo > yhi)
         {
         ihi = 1;
         ilo = 0;
         ynhi = yhi;
         yhi = ylo;
         ylo = ynhi;
         }

      for (i = 2; i < nvertex; i++)
         {
         yt = y[i] - T * log(rand->Next());
         if (yt <= ylo)
            {
            ilo = i;
            ylo = yt;
            }
         else if (yt > yhi)
            {
            ynhi = yhi;
            ihi = i;
            yhi = yt;
            }
         else if (yt > ynhi)
            ynhi = yt;
         }

      // Compute the fractional range from highest to lowest
      // and return if satisfactory
      // rtol = 2.0 * (yhi - ylo) / (fabs(yhi) + fabs(ylo));
      // if ((rtol < ftol) || (iter < 0))
      if (iter < 0)
         {
         // put best point and value in slot 1 before return
         y.Swap(0, ilo);
         simplex.SwapRows(0, ilo);
         return fmin;
         }

      // Begin a new iteration:
      // First. Extrapolate by a factor -1 through the face of the simplex
      //        across from the high point (ie. reflect from the high point)

      iter -= 2;

      ytry = Amoeba(ihi, -1.0);
      if (ytry < ylo)
         // gives a better result than the best point, so go even further
         Amoeba(ihi, 2.0);
      else if (ytry >= ynhi)
         {
         // the reflected point is worse than the second highest so
         // look for an intermediate (do a one-dimensional contration)
         ysave = yhi;
         ytry = Amoeba(ihi, 0.5);
         if (ytry >= ysave)
            {
            // Can't get rid of worst point, so contract around
            // lowest point
            for (i = 0; i < nvertex; i++)
               if ( i != ilo )
                  {
                  simplex[i].Add(simplex[ilo]);
                  simplex[i].Multiply(0.5);
                  y[i] = f(simplex[i]);

                  if (y[i] <= fmin)
                     {
                     point = simplex[i];
                     fmin = y[i];
                     }
                  }
            // Calculate psum
            for (psum = simplex[0], m = 1; m < nvertex; m++)
               psum.Add(simplex[m]);

            iter -= point.dim;
            }
         }
      else
         iter++;        // correct evaluation count
      }
   }

double SAMinimizer::Amoeba(int ihi, double factor)
   {
   double fac, yflu, ytry;

   fac = (1.0 - factor) / point.dim;
   ptry.SetMultiple(fac, psum);
   ptry.AddMultiple(factor - fac, simplex[ihi]);
   ytry = f(ptry);

   if (ytry <= fmin)
      {
      // save the best ever
      point = ptry;
      fmin = ytry;
      }

   // we added a thermal fluctuation to all
   // current vertices, but we subtract here
   // the simplex likes to accept any change
   yflu = ytry + T * log(rand->Next());
   if (yflu < yhi)
      {
      y[ihi] = ytry;
      yhi = yflu;
      psum.Subtract(simplex[ihi]);
      psum.Add(simplex[ihi] = ptry);
      }

   return yflu;
   }

// AmoebaMinimizer
//

AmoebaMinimizer::AmoebaMinimizer() : GeneralMinimizer(), psum(), ptry(), y()
   { cycleMax = 50000; }

void AmoebaMinimizer::Reset(int ndim, double scale)
   {
   GeneralMinimizer::Reset(ndim, scale);

   simplex.Dimension(ndim + 1, ndim);
   y.Dimension(ndim + 1);
   psum.Dimension(ndim);
   ptry.Dimension(ndim);
   }

double AmoebaMinimizer::Minimize(double ftol)
   {
   int i, ilo, ihi, inhi, m, nvertex = point.dim + 1;
   double rtol, ysave, ytry;

   if (point.dim == 0)
      return fmin = f(point);

   // Setup the simplex, and evaluate f at each vortex
   for (i=0; i < point.dim; i++)
      {
      simplex[i] = point;
      simplex[i].Add (directions[i]);
      y[i] = f(simplex[i]);
      if (y[i] < fmin) fmin = y[i];
      }

   simplex[nvertex - 1] = point;
   y[nvertex - 1] = f(simplex[nvertex - 1]);
   if (y[nvertex - 1] < fmin) fmin = y[nvertex - 1];

   cycleCount = nvertex;

   // Calculate psum
   for (psum = simplex[0], m = 1; m < nvertex; m++)
      psum.Add(simplex[m]);

   while (true)
      {
      // First we must determine which is the highest, next
      // highest and lowest point in the simplex
      ihi = y[0] > y[1] ? (ilo = inhi = 1, 0) : (ilo = inhi = 0, 1);

      for (i = 2; i < nvertex; i++)
         {
         if (y[i] <= y[ilo])
            ilo = i;
         else if (y[i] > y[ihi])
            {
            inhi = ihi;
            ihi = i;
            }
         else if (y[i] > y[inhi])
            inhi = i;
         }

      // Compute the range from highest to lowest, return if satisfactory
      rtol = 2 * fabs(y[ihi] - y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+ZEPS);
      if (rtol < ftol)
         {
         point = simplex[ilo];
         return fmin = y[ilo];
         }

      if (cycleCount > cycleMax)
         error("Amoeba.Minimize - Couldn't converge in %d cycles", cycleMax);

      // begin a new iteration...
      // first extrapolate by a factor of -1 through the face of the simplex
      // (reflection)

      cycleCount += 2;
      ytry = Amoeba(ihi, -1.0);

      if (ytry <= y[ilo])
         // Better than previous best, keep going
         Amoeba(ihi, 2.0);
      else if (ytry >= y[inhi])
         {
         // still worse than next highest point, so do a one-
         // dimensional contraction to find an intermediate
         // lower point
         ysave = y[ihi];
         ytry = Amoeba(ihi, 0.5);
         if (ytry >= ysave)
            {
            // Can't get rid of that high point
            // Better to contract around lowest point
            for (i = 0; i < nvertex; i++)
               if (i != ilo)
                  {
                  simplex[i].Add(simplex[ilo]);
                  simplex[i].Multiply(0.5);
                  y[i] = f(simplex[i]);
                  }

            cycleCount += point.dim;     // keep track of function evaluations

            for (psum = simplex[0], m = 1; m < nvertex; m++)
               psum.Add(simplex[m]);
            }
         }
      else
         cycleCount --;
      }
   }

double AmoebaMinimizer::Amoeba(int ihi, double factor)
   {
   double fac, ytry;

   fac = (1.0 - factor)/point.dim;
   ptry.SetMultiple(fac, psum);
   ptry.AddMultiple(factor - fac, simplex[ihi]);
   ytry = f(ptry);

   if (ytry < y[ihi])
      {
      // if this is better than the highest, replace it
      y[ihi] = ytry;
      psum.Subtract(simplex[ihi]);
      psum.Add(simplex[ihi] = ptry);
      }

   return ytry;
   }

// FletcherMinimizer
//

void FletcherMinimizer::Reset(int n, double )
   {
   point.Dimension(n);
   g.Dimension(n);
   h.Dimension(n);

   fmin = FPMAX;
   }

double FletcherMinimizer::Minimize(double ftol)
   {
   LineMinimizer linmin(*func);

   if (point.Length() == 0)
      return fmin = f(point);

   // Initialize direction vector
   linmin.line.Dimension(point.Length());

   // Evaluate function and derivatives
   fmin = f(point);

   while (true)
      {
      df(point, linmin.line, 0.1);

      // Minimizing so go against gradient
      linmin.line.Negate();
      h = linmin.line;
      g = linmin.line;

      iter = 1;
      while (true)
         {
         linmin.point = point;
         linmin.Bracket(0, 1);
         linmin.Brent(ftol);

         // If values tested while evaluating derivative
         // are closer to minimum than those found by
         // descending gradient, use those as starting points!

         if (func->dfmin < linmin.fmin)
            {
            linmin.point = point;
            linmin.line = func->dpmin;
            linmin.line.Subtract(point);
            linmin.Bracket(0,1);
            linmin.Brent(ftol);
            point = linmin.point;
            fmin = linmin.fmin;
            break;
            };

         point = linmin.point;

         if (2.0 * fabs(linmin.fmin - fmin) <= ftol * (fabs(linmin.fmin) + fabs(fmin) + ZEPS))
            return fmin = linmin.fmin;

         fmin = linmin.fmin;
         df(point, linmin.line, 0.1);

         double dgg = 0.0, gg = 0.0;

         for (int j = 0; j < linmin.line.Length(); j++)
            {
            gg  += g[j] * g[j];
            dgg += (linmin.line[j] + g[j]) * linmin.line[j];
            }

         // Gradient is exactly zero
         if (gg == 0.0)
            return fmin;

         double gamma = dgg / gg;

         for (int j = 0; j < linmin.line.Length(); j++)
            {
            g[j] = -linmin.line[j];
            h[j] = g[j] + gamma * h[j];
            linmin.line[j] = h[j];
            }

         // Normalize direction vector to avoid degenerate cases
         double norm = linmin.line.SumSquares();

         if (norm > 1.0)
            {
            norm = sqrt (1.0 / norm);
            g.Multiply(norm);
            h.Multiply(norm);
            linmin.line.Multiply(norm);
            }

          if (iter == ITMAX)
            error("Fletcher.Minimizer -- Failed to converge after %d iterations", ITMAX);

         iter++;
         }
      }
   }

// Evolutionary Minimizer
//

EvolutionaryMinimizer::EvolutionaryMinimizer()
   {
   Init(globalRandom);
   }

EvolutionaryMinimizer::EvolutionaryMinimizer(Random & randomSeries)
   {
   Init(randomSeries);
   }

void EvolutionaryMinimizer::Init(Random & randomSeries)
   {
   generate_random_points = true;

   max_generations = 400;

   crossover = 0.0;
   step_size = 0.5;
   multiples = 8;

   rand = &randomSeries;
   }

void EvolutionaryMinimizer::Reset(int dimensions, double scale)
   {
   GeneralMinimizer::Reset(dimensions, scale);

   int NP2 = 8 + dimensions * multiples * 2;

   points.Dimension(NP2, dimensions);
   y.Dimension(NP2);

   generate_random_points = true;
   }

double EvolutionaryMinimizer::Minimize(double /* ftol */)
   {
   int NP = point.dim * multiples + 4;
   int D  = point.dim;

   // Trivial case where D == 0
   if (D == 0)
      return fmin = f(point);

   // Generate a random population of points
   if (generate_random_points)
      {
      double scale = directions[0].SumSquares();

      for (int i = 0; i < NP; i++)
         for (int j = 0; j < D; j++)
            points[i][j] = point[j] + rand->Uniform(-scale, scale);
      }

   // Evaluate the function at each point
   for (int i = 0; i < NP; i++)
      y[i] = f(points[i]);

   generations = 0;

   // These variables are used to decide when to stop early
   // int    strikes = 0;
   double optimum = fmin;

   // Loop until done
   while (generations < max_generations)
      {
      // Process each vector in the population
      for (int i = 0, a, b, c; i < NP; i++)
         {
         // Pick ancestors for the new point
         do { a = rand->NextInt() % NP; } while (a == i);
         do { b = rand->NextInt() % NP; } while (b == a || b == i);
         do { c = rand->NextInt() % NP; } while (c == b || c == a || c == i);

         int mutation = rand->NextInt() % D;

         // generate a new point
         for (int j = 0; j < D; j++)
            if (j == mutation || rand->Next() < crossover)
               points[NP + i][j] = points[c][j] + step_size * (points[a][j] - points[b][j]);
            else
               points[NP + i][j] = points[i][j];

         // Evaluate function at the new point
         y[NP + i] = f(points[NP + i]);
         }

      // Select improved solutions into the population
      for (int i = 0; i < NP; i++)
         if (y[i] > y[NP + i])
            points.SwapRows(i, NP + i),
            y[i] = y[NP + i];

      // For debugging purposes, print the minimum after each round
      // printf("\t\tMinimum after %d rounds is %.3f\n", generations, y.Min());

      generations++;

/*    double best = y.Min();
      double worst = y.Max();

      if ( (worst - best) < (fabs(worst + best) * ftol + ftol))
         { if (++strikes == 3) { printf("\nHIT!!\n\n"); break;} }
      else strikes = 0; */

      fmin = optimum;
      }

   // Retrieve the best point from the population
   int min = 0;

   for (int i = 1; i < NP; i++)
      if (y[i] < y[min])
         min = i;

   point = points[min];
   return fmin = y[min];
   }



 
