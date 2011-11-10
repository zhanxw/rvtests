////////////////////////////////////////////////////////////////////// 
// libsrc/MathFloatVector.cpp 
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
 
#include "MathFloatVector.h"
#include "MathVector.h"
#include "MathMatrix.h"
#include "MathConstant.h"
#include "Sort.h"
#include "Error.h"

#ifdef  _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <string.h>
#include <math.h>

void FloatVector::Init()
   {
   dim = size = 0;
   label = "Unknown";
   data = NULL;
   }

FloatVector::~FloatVector()
   {
   // printf(" Deleting vector %s ...\n", (const char *) label);
   if (data != NULL) delete [] data;
   }

void FloatVector::Dimension(int d)
   {
   if (d > size)
      if (size < 1024)
         {
         size = (d + Vector::alloc) / Vector::alloc * Vector::alloc;
         float * newData = new float [size];
         if (data != NULL)
            {
            for (int i = 0; i < dim; i++)
               newData[i] = data[i];
            delete [] data;
            }
         data = newData;
         }
      else
         {
         while (size <= d)
            size *= 2;

         float * newData = new float [size];
         if (data != NULL)
            {
            for (int i = 0; i < dim; i++)
               newData[i] = data[i];
            delete [] data;
            }
         data = newData;
         }
   dim = d;
   }

void FloatVector::Negate()
   {
   for (int i = 0; i < dim; i++)
      data[i] = -data[i];
   }

void FloatVector::Add(double n)
   {
   for (int i = 0; i< dim; i++)
      data[i] += n;
   }

void FloatVector::Multiply(double k)
   {
   for (int i = 0; i < dim; i++)
      data[i] *= k;
   }

void FloatVector::Copy(const FloatVector & v)
   {
   Dimension(v.dim);

   if (v.data != NULL)
      for (int i=0; i < dim; i++)
         data[i] = v.data[i];
   }

FloatVector & FloatVector::operator = (const FloatVector & rhs)
   {
   Copy(rhs);
   return *this;
   }

void FloatVector::Add(FloatVector & v)
   {
   if (dim != v.dim)
      error("FloatVector::Add - vectors have different dimensions\n"
            "Vectors     - %s [%d] + %s [%d] ",
            (const char *) label, dim, (const char  *) v.label, v.dim);

   for (int i = 0; i < dim; i++)
      data[i] += v.data[i];
   }

void FloatVector::AddMultiple(double k, FloatVector & v)
   {
   if (dim != v.dim)
      error("FloatVector::AddMultiple - vectors are incompatible\n"
            "Vectors             - %s [%d] + %s [%d] ",
            (const char  *) label, dim, (const char  *) v.label, v.dim);

   for (int i = 0; i < dim; i++)
      data[i] += k * v.data[i];
   }


void FloatVector::Subtract(FloatVector & v)
   {
   if (dim != v.dim)
      error("FloatVector::Subtract - vectors have different dimensions\n"
            "Vectors          - %s [%d] + %s [%d] ",
            (const char  *) label, dim, (const char  *) v.label, v.dim);

   for (int i = 0; i < dim; i++)
      data[i] -= v.data[i];
   }


void FloatVector::Zero()
   {
   for (int i = 0; i < dim; i++)
      data[i] = 0.0;
   }

void FloatVector::Set(double k)
   {
   for (int i = 0; i < dim; i++)
      data[i] = k;
   }

void FloatVector::SetMultiple(double k, FloatVector & v)
   {
   Dimension(v.dim);

   for (int i = 0; i < dim; i++)
      data[i] = k * v[i];
   }

double FloatVector::InnerProduct(FloatVector & v)
   {
   if (dim != v.dim)
      error("FloatVector::InnerProduct - vectors have different dimensions\n"
            "Vectors              - %s[%d] * %s[%d] ",
            (const char  *) label, dim, (const char  *) v.label, v.dim);

   double sum = 0.0;
   for (int i = 0; i < dim; i++)
      sum += data[i] * v.data[i];

   return sum;
   }

void FloatVector::Insert(int n, double value)
   {
   Dimension(dim + 1);

   for (int i = dim - 1; i > n; i--)
      data[i] = data[i - 1];
   data[n] = value;
   }

void FloatVector::DeleteDimension(int n)
   {
   for (int i = n; i < dim - 1; i++)
      data[i] = data[i + 1];
   dim --;
   }

void FloatVector::Product(Matrix & m, FloatVector & v)
   {
   if (m.cols != v.dim)
      error ("FloatVector::Product - Cannot Multiply Matrix by FloatVector\n"
             "Vectors         - %s [%d, %d] * %s [%d]\n",
             (const char  *) m.label, m.rows, m.cols,
             (const char  *) v.label, v.dim);

   Dimension(m.rows);
   Zero();

   for(int i = 0; i < m.rows; i++)
      for (int j = 0; j < m.cols; j++)
         data[i] += m[i][j] * v[j];
   }

double FloatVector::Average() const
   {
   if (dim == 0)
      error("Average undefined for null FloatVector %s",
            (const char  *) label);

   return Sum() / dim;
   }

double FloatVector::Product() const
   {
   double product = 1.0;

   for (int j = 0; j < dim; j++)
      product *= data[j];

   return product;
   }

double FloatVector::Sum() const
   {
   double sum = 0.0;

   for (int j=0; j<dim; j++)
      sum += data[j];

   return sum;
   }

double FloatVector::SumSquares() const
   {
   double sum = 0.0;

   for (int j=0; j<dim; j++)
      sum += data[j] * data[j];

   return sum;
   }

void FloatVector::AveVar(double & ave, double & var) const
   {
   // uses a two pass method to correct for
   // round-off errors

   if (dim == 0)
      error("Average and Variance undefined for null vector %s",
            (const char  *) label);

   double s, ep;

   ave = var = ep = 0.0;

   for (int j=0; j<dim; j++)
      ave += data[j];

   ave /= dim;

   for (int j=0; j<dim; j++)
      {
      s = data[j] - ave;
      ep += s;
      var += s*s;
      }

   if (dim > 1)
      var = (var - ep*ep/dim)/(dim-1);
   }

double FloatVector::Var() const
   {
   double mean, var;
   AveVar(mean, var);
   return var;
   }

void FloatVector::Print(FILE * f, int d)
   {
   if (d == -1 || d > dim) d = dim;

   fprintf(f, "%.15s : ", (const char  *) label);
   for (int i = 0; i < d; i++)
      fprintf(f, "%7.3f ", data[i]);
   fprintf(f, "\n");
   }

int FloatVector::CompareFloat(const float * a, const float * b)
   {
   if (*a < *b) return -1;
   if (*a > *b) return 1;
   return 0;
   }

void FloatVector::Sort()
   {
   QuickSort(data, dim, sizeof(float), COMPAREFUNC CompareFloat);
   }

void FloatVector::Sort(FloatVector & freeRider)
   {
   QuickSort2(data, freeRider.data, dim, sizeof(float),
              COMPAREFUNC CompareFloat);
   }

int FloatVector::BinarySearch(double element)
   {
   void * pointer = ::BinarySearch
         (&element, data, dim, sizeof(float), COMPAREFUNC CompareFloat);

   if (pointer == NULL)
      return -1;

   return ((float *) pointer) - data;
   }

void FloatVector::RemoveDuplicates()
   {
   int out = 0;

   for (int in = 1; in < Length(); in++)
      if (data[in] != data[out])
         data[++out] = data[in];

   Dimension(out + 1);
   }

bool FloatVector::operator == (const FloatVector & rhs) const
   {
   if (rhs.dim != dim) return false;

   for (int i = 0; i < dim; i++)
      if (data[i] != rhs[i])
         return false;
   return true;
   }

// These functions are useful for simulation
//

int FloatVector::CountIfGreater(double threshold) const
   {
   int count = 0;

   for (int i = 0; i < dim; i++)
      if (data[i] > threshold)
         count++;

   return count;
   }

int FloatVector::CountIfGreaterOrEqual(double treshold) const
   {
   int count = 0;

   for (int i = 0; i < dim; i++)
      if (data[i] >= treshold)
         count++;

   return count;
   }

// Min and max functions
//

double FloatVector::Min() const
   {
   if (dim == 0)
      return 0.0;

   double min = data[0];

   for (int i = 1; i < dim; i++)
      if (data[i] < min)
         min = data[i];

   return min;
   }

double FloatVector::Max() const
   {
   if (dim == 0)
      return 0.0;

   double max = data[0];

   for (int i = 1; i < dim; i++)
      if (data[i] > max)
         max = data[i];

   return max;
   }

// Push and Pop functions for using vector as a stack
//

void FloatVector::Push(double value)
   {
   Dimension(dim + 1);
   data[dim - 1] = value;
   }

void FloatVector::Stack(const FloatVector & v)
   {
   int end = dim;

   Dimension(dim + v.dim);

   for (int i = 0; i < v.dim; i++)
      data[i + end] = v[i];
   }

// Check if all values are in ascending or descending order
//

bool FloatVector::isAscending()
   {
   for (int i = 1; i < dim; i++)
      if (data[i] < data[i - 1])
         return false;
   return true;
   }

bool FloatVector::isDescending()
   {
   for (int i = 1; i < dim; i++)
      if (data[i] > data[i - 1])
         return false;
   return true;
   }

int FloatVector::SafeCount() const
   {
   int nonMissing = dim;

   for (int i = 0; i < dim; i++)
      if (data[i] == _NAN_)
         nonMissing--;

   return nonMissing;
   }

double FloatVector::SafeMin() const
   {
   double min = _NAN_;
   int i;

   for (i = 0; i < dim; i++)
      if (data[i] != _NAN_)
       {
       min = data[i];
       break;
       }

   for (; i < dim; i++)
      if (data[i] != _NAN_ && data[i] < min)
         min = data[i];

   return min;
   }

double FloatVector::SafeMax() const
   {
   double max = _NAN_;
   int i;

   for (i = 0; i < dim; i++)
      if (data[i] != _NAN_)
       {
       max = data[i];
         break;
       }

   for (; i < dim; i++)
      if (data[i] != _NAN_ && data[i] > max)
         max = data[i];

   return max;
   }

void FloatVector::Reverse()
   {
   for (int i = 0, j = dim - 1; i < j; i++, j--)
      Swap(i, j);
   }

void FloatVector::InsertInSortedList(int value)
   {
   // Skip through large elements
   int pos = dim - 1;

   while (pos >= 0 && data[pos] > value)
      pos--;

   // If the value is already in the list, we are done
   if (pos >= 0 && data[pos] == value)
      return;

   // Otherwise we need to grow array
   Dimension(dim + 1);

   // And then shift larger elements to the right
   pos++;
   for (int i = dim - 1; i > pos; i--)
      data[i] = data[i - 1];

   data[pos] = value;
   }

void FloatVector::Swap(FloatVector & rhs)
   {
   float * temp = rhs.data;
   rhs.data = data;
   data = temp;

   int swap = rhs.dim;
   rhs.dim = dim;
   dim = swap;

   swap = rhs.size;
   rhs.size = size;
   size = swap;
   }


 
