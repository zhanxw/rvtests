//: src:kbac.h
// kbac test
// Copyright 2011 Gao Wang
#ifndef _KBAC_H        
#define _KBAC_H 
#include <iostream>
#include <list>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <iomanip>

class KbacTest 
{

  public:
    KbacTest(int* nn, int* qq, double* aa, double* mafUpper, double* xdatIn, double* ydatIn, double* mafIn, int* xcol, int* ylen);
    ~KbacTest();

    void calcKbacP(double* pvalue, int* sided);

  private:
    std::vector< std::vector<double> > __xdat;
    std::vector<double> __ydat;
    std::vector<double> __observedMafs;
    double __mafLower;
    double __mafUpper;
    bool __quiet;
    double __alpha;
    double __nPermutations;
    unsigned int __adaptive;

    double m_checkAdaptivePvalue(unsigned int permcount1, unsigned int permcount2, unsigned int currentIdx, unsigned int checkPoint, unsigned int alternative) const;   
    void m_trimXdat();
};

namespace std {

//!- Dump a vector to screen
template<class T> ostream & operator<<(ostream & out, const vector<T> & vec)
  {
    if (!vec.empty()) {
      typename vector<T>::const_iterator it = vec.begin();
      out << *it;
      for (++it; it != vec.end(); ++it)
        out << " " << *it ;
    }
    return out;
  }
}
#endif
