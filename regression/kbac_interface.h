//: src:kbac_interface.h
// kbac class interface
// Copyright 2011 Gao Wang
#ifndef _KBAC_INTERFACE_H        
#define _KBAC_INTERFACE_H 
void set_up_kbac_test(int* nn, int* qq, double* aa, double* mafUpper, double* xdatIn, double* ydatIn, double* mafIn, int* xcol, int* ylen); 
void do_kbac_test(double* pvalue, int* twosided);
void clear_kbac_test();
#endif
