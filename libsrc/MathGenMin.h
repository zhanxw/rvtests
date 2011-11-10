////////////////////////////////////////////////////////////////////// 
// libsrc/MathGenMin.h 
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
 
#ifndef __MATHPOWELL_H__
#define __MATHPOWELL_H__

#include "MathGold.h"
#include "MathVector.h"
#include "MathMatrix.h"
#include "Random.h"

// Multidimensional minimization of a continuous function
// starting with a user supplied starting point and
// direction vector
//
class GeneralMinimizer
   {
   public:
      VectorFunc * func;               // Function to be minimized
      Matrix   directions;
      Vector   point;
      double   fmin;

      // Setup matrices assuming ndim point
      virtual void   Reset(int ndim, double scale = 1.0);

      // Find a minimum using direction set and starting point
      virtual double Minimize(double ftol = TOL) = 0;

      GeneralMinimizer();
      virtual ~GeneralMinimizer() { }

      double f(Vector & v)
         { return func->Evaluate(v); }

      void df(Vector & v, Vector & d, double scale = 1.0)
         { func->Derivative(v, d, scale); }
   };

// Powell's conjugate direction method
// After each round, the direction of largest decrease is replaces
// its biggest component among the original directions
//

class PowellMinimizer : public GeneralMinimizer
   {
   public:
      int iter;

      virtual ~PowellMinimizer() { }
      virtual double Minimize(double ftol = TOL);
   };


// Simulated annealing using simplex method of Nelder and Mead
//
class SAMinimizer : public GeneralMinimizer
   {
   public:
      int      iter;
      bool     freeRand;
      Random * rand;

      Vector   y;                   // evaluation of entropy at y
      Matrix   simplex;             // volume in n dimensions (n+1) points

      SAMinimizer();
      SAMinimizer(Random & rand);

      virtual ~SAMinimizer();

      virtual void   Reset(int ndim, double scale = 1.0);

      // Lowers temperature T from maxT to minT in Tcycles linear decay cycles
      // Titer iterations at each temperature
      virtual double Minimize(double ftol = TOL);
      double MinimizeLoop(double ftol = TOL);

      double   T, maxT, minT;      // Temperature
      int      Tcycles, Titer;     // Cycling parameters

   private:
      Vector psum;
      Vector ptry;
      double yhi;

      void   Constructor();
      double Amoeba(int ihi, double factor);
   };

// Multidimensional minimization of a continuous function by
// the down-hill simplex method of Nelder and Mead
// (Computer Journal 1965)
//
class AmoebaMinimizer : public GeneralMinimizer
   {
   public:
      Matrix      simplex;
      long        cycleCount, cycleMax;       // number of function evaluations

      AmoebaMinimizer();
      virtual ~AmoebaMinimizer() { }

      virtual void Reset(int dimensions, double scale = 1.0);
      virtual double Minimize(double ftol = TOL);

   private:
      Vector psum, ptry, y;

      double Amoeba(int ihi, double factor);
   };

// Differential Evolution minimizer
// A stochastic minimizer based on the algorithm of Storn and Price, 1996

class EvolutionaryMinimizer : public GeneralMinimizer
   {
   public:
      Matrix   points;
      Vector   y;

      double crossover;       // This is the CR parameter of Storn and Price
      double step_size;       // This is the L parameter of Storn and Price
      int    multiples;       // The NP paraemter of Storn and Price will be dimensions * multiple

      Random * rand;

      bool   generate_random_points;

      int    generations;
      int    max_generations;

      EvolutionaryMinimizer();
      EvolutionaryMinimizer(Random & randomSeries);

      ~EvolutionaryMinimizer() { }

      virtual void Reset(int dimensions, double scale = 1.0);
      virtual double Minimize(double ftol = TOL);

   private:
      void Init(Random & randomSeries);
   };

// Conjugate gradient minimizer
// Polak-Ribiere improvement on Fletcher-Reeves algorithm for
// multidimensional minimization.
//

class FletcherMinimizer : public GeneralMinimizer
   {
   public:
      int iter;

     FletcherMinimizer() { }

      virtual void    Reset(int dimensions, double scale = 1.0);
      virtual double  Minimize(double ftol = TOL);

   private:
      Vector g, h;
   };

#endif


 
