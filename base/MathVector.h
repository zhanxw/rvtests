//////////////////////////////////////////////////////////////////////
// libsrc/MathVector.h
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

#ifndef __MATHVECTOR_H__
#define __MATHVECTOR_H__
#include <vector>

#define DECLARE_EIGEN_VECTOR(matRef, varName)               \
  Eigen::Map<Eigen::MatrixXd> varName((matRef).data.data(), \
                                      (matRef).data.size(), 1);

#define DECLARE_EIGEN_CONST_VECTOR(matRef, varName)               \
  Eigen::Map<const Eigen::MatrixXd> varName((matRef).data.data(), \
                                            (matRef).data.size(), 1);

class Vector {
 public:
  Vector() {}
  Vector(int n) { data.resize(n); }
  double operator[](int idx) const { return data[idx]; }
  double& operator[](int idx) { return data[idx]; }
  int Length() const { return data.size(); }
  void Dimension(int n);
  void Dimension(int n, double val);
  void Fill(double val);
  void Zero() { Fill(0); }
  double Sum() const;
  double Average() const;
  double Min() const;
  double Max() const;
#if 0
  operator Eigen::Map<Eigen::MatrixXd>() {
    Eigen::Map<Eigen::MatrixXd> ret(data.data(), data.size(), 1);
    return ret;
  }
  operator Eigen::Map<const Eigen::MatrixXd>() const {
    Eigen::Map<const Eigen::MatrixXd> ret(data.data(), data.size(), 1);
    return ret;
  }
#endif
 public:
  std::vector<double> data;
};

#if 0
#include "StringBasics.h"

#include <assert.h>
#include <stdio.h>

class Matrix;

class Vector {
 public:
  int dim, size;
  double *data;
  String label;

  Vector() { Init(); }
  Vector(Vector &v) {
    Init();
    Copy(v);
  }
  Vector(int d) {
    Init();
    Dimension(d);
  }
  Vector(const char *text) {
    Init();
    label = text;
  }
  Vector(const char *text, int d) {
    Init();
    label = text;
    Dimension(d);
  }
  Vector(const char *text, Vector &v) {
    Init();
    label = text;
    Copy(v);
  }

  ~Vector();

  void Dimension(int d);
  void Dimension(int d, double value);

  void GrowTo(int d) { Dimension(d > dim ? d : dim); }
  void GrowTo(int d, double value) { Dimension(d > dim ? d : dim, value); }

  int Length() const { return dim; }

  void SetLabel(const char *text) { label = text; }

  void Zero();
  void Set(double k);
  void Set(Vector &v) { Copy(v); };
  void SetMultiple(double k, Vector &v);

  void Negate();
  void Add(double n);
  void Multiply(double k);

  double InnerProduct(Vector &v);
  void Copy(const Vector &v);
  void Add(const Vector &v);
  void AddMultiple(double k, const Vector &v);
  void Subtract(Vector &v);

  void Product(Matrix &m, Vector &v);

  double &operator[](int n) {
    assert(n < dim);
    return data[n];
  }
  double operator[](int n) const {
    assert(n < dim);
    return data[n];
  }

  double operator[](double fraction) { return data[(int)(dim * fraction)]; }
  double &operator[](double fraction) const {
    return data[(int)(dim * fraction)];
  }

  Vector &operator=(const Vector &v);
  bool operator==(const Vector &v) const;
  bool operator!=(const Vector &v) const { return !(*this == v); }

  void Swap(int i, int j) {
    double swap = data[i];
    data[i] = data[j];
    data[j] = swap;
  }
  void Swap(Vector &rhs);

  Vector &operator*=(double rhs) {
    Multiply(rhs);
    return *this;
  }
  Vector &operator+=(double rhs) {
    Add(rhs);
    return *this;
  }
  Vector &operator-=(double rhs) { return *this += -rhs; }
  Vector &operator/=(double rhs) { return *this *= 1 / rhs; }

  Vector &operator+=(Vector &rhs) {
    Add(rhs);
    return *this;
  }
  Vector &operator-=(Vector &rhs) {
    Subtract(rhs);
    return *this;
  }

  void DeleteDimension(int n);
  void Delete(int n) { DeleteDimension(n); }
  void Insert(int n, double value);

  // Calculates average and variance
  void AveVar(double &ave, double &var) const;
  double Average() const;
  double Var() const;

  double Average(double returnIfNull);
  double Var(double returnIfNull);

  // Common descriptive functions
  double Sum() const;
  double SumSquares() const;
  double Product() const;

  // Find extreme values
  double Min() const;
  double Max() const;

  // Return the number of elements in a subset
  int CountIfGreater(double treshold) const;
  int CountIfGreaterOrEqual(double treshold) const;

  // Append another vector to the end
  void Stack(const Vector &v);

  void Print(int maxDim = -1) { Print(stdout, maxDim); }
  void Print(FILE *output, int maxDim = -1);

  // Routines for creating and searching through sorted vectors
  void Sort();
  void Reverse();
  void Sort(Vector &freeRider);
  int BinarySearch(double value);
  int FastFind(double value) { return BinarySearch(value); }

  // Remove consecutive duplicate elements from vector
  void RemoveDuplicates();

  // Query first and last elements
  //

  double &First() { return data[0]; }
  double &Last() { return data[dim - 1]; }

  // Routines for using a vector as a stack of doubles
  //

  void Clear() { dim = 0; }
  void Push(double value);
  double Pop() { return data[--dim]; }
  double Peek() const { return data[dim - 1]; }

  // This routine adds items to a sorted list
  //

  void InsertInSortedList(int item);

  static int alloc;

  bool isAscending();
  bool isDescending();

  // Routines for dealing with vectors that include missing data
  //

  int SafeCount() const;
  double SafeMin() const;
  double SafeMax() const;

 private:
  static int CompareDouble(const double *a, const double *b);
  void Init();
};

class VectorFunc
// Wrapper for multi-dimensional functions
// so that they can be used as parameters
// and keep private data
{
 private:
  double (*f)(Vector &);

 public:
  // Constructors
  VectorFunc();
  VectorFunc(double (*func)(Vector &));

  // Virtual destructor ensures that dynamic objects are
  // handled correctly
  virtual ~VectorFunc() {}

  virtual double Evaluate(Vector &v);

  // Calculate derivatives along each direction. Delta is a guess value
  // for the initial stepsize in numerical derivation
  virtual void Derivative(Vector &point, Vector &d, double delta = 1.0);

  // Minimum function value found while evaluating derivative
  // and its location...
  double dfmin;
  Vector dpmin;
};
#endif
#endif
