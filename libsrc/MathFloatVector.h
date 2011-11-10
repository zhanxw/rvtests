////////////////////////////////////////////////////////////////////// 
// libsrc/MathFloatVector.h 
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
 
#ifndef __MATHFLOATVECTOR_H__
#define __MATHFLOATVECTOR_H__

#include "StringBasics.h"

#include <stdio.h>
#include <assert.h>

class Matrix;

class FloatVector
   {
   public:
      int         dim, size;
      float *    data;
      String      label;

   FloatVector()
      { Init(); }
   FloatVector(FloatVector & v)
      { Init(); Copy(v); }
   FloatVector(int d)
      { Init(); Dimension(d); }
   FloatVector(const char * text)
      { Init(); label = text; }
   FloatVector(const char * text, int d)
      { Init(); label = text; Dimension(d); }
   FloatVector(const char * text, FloatVector & v)
      { Init(); label = text; Copy(v); }

   ~FloatVector();

   void   Dimension(int d);
   int    Length() const { return dim; }

   void   SetLabel(const char * text) { label = text; }

   void   Zero();
   void   Set(double k);
   void   Set(FloatVector & v) { Copy(v); };
   void   SetMultiple(double k, FloatVector & v);

   void   Negate();
   void   Add(double n);
   void   Multiply(double k);

   double InnerProduct(FloatVector & v);
   void   Copy(const FloatVector & v);
   void   Add(FloatVector & v);
   void   AddMultiple(double k, FloatVector & v);
   void   Subtract(FloatVector & v);

   void Product(Matrix & m, FloatVector & v);

   float & operator [] (int n)
      { assert(n < dim);  return data[n]; }
   float operator [] (int n) const
      { assert(n < dim);  return data[n]; }

   float operator [] (double fraction)
      { return data[(int) (dim * fraction)]; }
   float & operator [] (double fraction) const
      { return data[(int) (dim * fraction)]; }

   FloatVector & operator = (const FloatVector & v);
   bool operator == (const FloatVector & v) const;
   bool operator != (const FloatVector & v) const { return !(*this == v); }

   void Swap(int i, int j)
      { double swap = data[i]; data[i] = data[j]; data[j] = swap; }
   void Swap(FloatVector & rhs);

   FloatVector & operator *= (double rhs) { Multiply(rhs); return *this; }
   FloatVector & operator += (double rhs) { Add(rhs); return *this; }
   FloatVector & operator -= (double rhs) { return *this += -rhs; }
   FloatVector & operator /= (double rhs) { return *this *= 1/rhs; }

   void DeleteDimension (int n);
   void Delete(int n) { DeleteDimension(n); }
   void Insert(int n, double value);

   // Calculates average and variance
   void   AveVar(double & ave, double & var) const;
   double Average() const;
   double Var() const;

   // Common descriptive functions
   double Sum() const;
   double SumSquares() const;
   double Product() const;

   // Find extreme values
   double Min() const;
   double Max() const;

   // Return the number of elements in a subset
   int  CountIfGreater(double treshold) const;
   int  CountIfGreaterOrEqual(double treshold) const;

   // Append another vector to the end
   void Stack(const FloatVector & v);

   void Print(int maxDim = -1) { Print(stdout, maxDim); }
   void Print(FILE * output, int maxDim = -1);

   // Routines for creating and searching through sorted vectors
   void Sort();
   void Reverse();
   void Sort(FloatVector & freeRider);
   int  BinarySearch(double element);

   // Remove consecutive duplicate elements from FloatVector
   void RemoveDuplicates();

   // Query first and last elements
   //

   float & First() { return data[0]; }
   float & Last()  { return data[dim - 1]; }

   // Routines for using a vector as a stack of doubles
   //

   void   Clear()      { dim = 0; }
   void   Push(double value);
   double Pop()        { return data[--dim]; }
   double Peek() const { return data[dim-1]; }

   // This routine adds items to a sorted list
   //

   void   InsertInSortedList(int item);

   bool   isAscending();
   bool   isDescending();

   // Routines for dealing with vectors that include missing data
   //

   int SafeCount() const;
   double SafeMin() const;
   double SafeMax() const;

   private:
      static int CompareFloat(const float * a, const float * b);
      void Init();
   };

#endif



 
