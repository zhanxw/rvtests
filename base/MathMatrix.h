//////////////////////////////////////////////////////////////////////
// libsrc/MathMatrix.h
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

#ifndef __MATHMATRIX_H__
#define __MATHMATRIX_H__

#include <set>
#include <string>

#include "MathVector.h"

#define DECLARE_EIGEN_MATRIX(matRef, varName)                              \
  Eigen::Map<Eigen::MatrixXd> varName((matRef).data.data(), (matRef).rows, \
                                      (matRef).cols);
#define DECLARE_EIGEN_CONST_MATRIX(matRef, varName)               \
  Eigen::Map<const Eigen::MatrixXd> varName((matRef).data.data(), \
                                            (matRef).rows, (matRef).cols);

// column-major
class Matrix {
 public:
  Matrix();
  Matrix(int nr, int nc);
  Matrix(const Matrix& m);
  Matrix& operator=(const Matrix& m);
  double& operator()(int i, int j) { return data[i + j * rows]; }
  double operator()(int i, int j) const { return data[i + j * rows]; }
  /**
   * Resize matrix to @param nr by @param nc
   * Existing contents is preserved
   */
  void Dimension(int nr, int nc);

  /**
   * Set all matrix elements to @param val
   */
  void Dimension(int nr, int nc, double val);
  void DimensionQuick(int nr, int nc);
  void Reserve(int nr, int nc) { data.reserve(nr * nc); }
  void Fill(double val) { std::fill(data.begin(), data.end(), val); }
  void Ones(int nr, int nc) {
    rows = nr;
    cols = nc;
    colLabel.resize(nc);
    Fill(1.0);
  }
  void Zero() { Fill(0.); }
  const std::string& GetColumnLabel(int idx) const { return colLabel[idx]; }
  void SetColumnLabel(int idx, const std::string& label) {
    colLabel[idx] = label;
  }
  double Min() const;
  double Max() const;
  void Product(const Matrix& in1, const Matrix& in2);
  void Transpose(const Matrix& old);
  /**
   * Each element is multiplied by @param s
   */
  Matrix& Multiply(double s);
  /**
   * @param rowIndexToRemove each index should be in [0, ... , rows-1]
   * @return number of row deleted
   */
  int RemoveByRowIndex(const std::vector<int>& rowIndexToRemove);

  /**
   * Stack @param m to the right
   */
  Matrix& StackRight(const Matrix& m);

  Matrix& MakeColumnMatrix(const Vector& v) {
    rows = v.Length();
    cols = 1;
    data = v.data;
    colLabel.resize(1);
    return *this;
  }
#if 0
  operator Eigen::Map<Eigen::MatrixXd>() {
    Eigen::Map<Eigen::MatrixXd> ret(data.data(), rows, cols);
    return ret;
  }
  operator Eigen::Map<const Eigen::MatrixXd>() const {
    Eigen::Map<const Eigen::MatrixXd> ret(data.data(), rows, cols);
    return ret;
  }
#endif
 public:
  int rows;
  int cols;
  std::vector<double> data;
  std::vector<std::string> colLabel;
};

#if 0
#include "Error.h"
#include "MathVector.h"

#include <stdio.h>

class ColumnExtras
   {
   private:
      bool dirty;
      int  precision, width;

      void Init();
      void Copy(ColumnExtras & c);

   public:
      String label;

      ColumnExtras()
         { Init(); }
      ColumnExtras(ColumnExtras & original)
         { Init(); Copy(original); }
      ~ColumnExtras();

      void SetLabel(const char * name);
      void SetPrecision(int p)
         {
         precision = p;
         dirty = true;
         }
      void SetWidth(int w)
         {
         width = w;
         dirty = true;
         }

      int    GetWidth();
      int    GetPrecision()
         { return precision; }

      ColumnExtras & operator = (ColumnExtras & rhs)
         { Copy(rhs); return (*this); }

      void Swap(ColumnExtras & rhs);
   };

class Matrix
   {
   public:
      String    label;
      ColumnExtras * extras;
      int       rows, cols, size, extraSize;
      Vector ** data;

   Matrix()
      { Init(); }
   Matrix(Matrix & m)
      { Init(); Copy(m); }
   Matrix(Matrix & m, const char * name)
      { Init(); Copy(m); SetLabel(name); }
   Matrix(int n, int m)
      { Init(); Dimension(n, m); }
   Matrix(const char * name)
      { Init(); SetLabel(name); }
   Matrix(const char * name, int n, int m)
      { Init(); Dimension(n, m); SetLabel(name); }
   ~Matrix();

   void Dimension(int m, int n);
   void Dimension(int m, int n, double value);
   void GrowTo(int m, int n)
      { Dimension(m > rows ? m : rows, n > cols ? n : cols); }
   void GrowTo(int m, int n, double value)
      { Dimension(m > rows ? m : rows, n > cols ? n : cols, value); }

   void SetLabel(const char * name);
   void SetColumnLabel(int n, const char * name)
      { extras[n].SetLabel(name); }
   const char * GetColumnLabel(int n)
      { return extras[n].label; }
   void SetColWidth(int n, int w)
      { extras[n].SetWidth(w); }
   void SetColPrecision(int n, int p)
      { extras[n].SetPrecision(p); }
   void CopyLabels(Matrix & m);

   void Negate();
   void Identity();
   void Zero();
   void Set(double k);

   void Copy(const Matrix & m);
   void Transpose(const Matrix & m);
   void Add(const Matrix & m);
   void AddMultiple(double k, const Matrix & m);
   void Product(const Matrix & left, const Matrix & right);

   void Add(double k);
   void Multiply(double k);

   // Reduces a matrix to row echelon form, assuming
   // values smaller than tol are zero
   void Reduce(double tol = 0.0);

   Vector & operator [] (int i)
      { assert(i < rows); return *(data[i]); }

   const Vector & operator [] (int i) const
      { assert(i < rows); return *(data[i]); }

   void DeleteRow(int r);
   void DeleteColumn(int c);

   void SwapRows(int r1, int r2)
      { Vector * temp = data[r1];
        data[r1] = data[r2];
        data[r2] = temp;
      };

   void SwapColumns(int c1, int c2);

   void MultiplyRow(int r1, double k);
   void AddRows(int r1, int r2);
   void AddRows(double k, int r1, int r2);

   // Sort according to numeric values in the first column
   void Sort();

   void Print(FILE * f, int maxRows = -1, int maxCols = -1);
   void PrintUpper(FILE * f, int maxRows = -1, int maxCols = -1, bool print_diag = false);
   void PrintLower(FILE * f, int maxRows = -1, int maxCols = -1, bool print_diag = false);
   void SetupPrint(FILE *f, int r, int c, int & column_zero, int * precision, int * width);

   void Read(FILE * f);

   Matrix & operator = (const Matrix & rhs)
      { Copy(rhs); return *this; }

   bool operator == (const Matrix & rhs) const;
   bool operator != (const Matrix & rhs) const { return !(*this == rhs); }

   Matrix & operator *= (double rhs)
      { Multiply(rhs); return *this; }
   Matrix & operator /= (double rhs)
      { Multiply(1.0/rhs); return *this; }

   // Stack a matrix to the bottom of the current matrix
   void StackBottom(const Matrix & m);

   // Stack a matrix to the left of the current matrix
   void StackLeft(const Matrix & m);

   // Swap dynamic allocation for two matrices
   void Swap(Matrix & m);

   // Functions that calculate basic summary statistics
   double Min() const;
   double Max() const;
   double Mean() const;

   // Functions that calculate summary statistics in the presence of missing data
   double SafeMin() const;
   double SafeMax() const;
   double SafeMean() const;
   int SafeCount() const;

   // Return the last row in matrix
   Vector & Last() { return *(data[rows - 1]); }

   private:
      static int alloc;
      static int CompareRows(Vector ** row1, Vector ** row2);

   void Init();
   };
#endif
#endif
