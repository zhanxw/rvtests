//: src:kbac_interface.cpp
// kbac class interface
// Copyright 2011 Gao Wang

#include "kbac.h"
#include "kbac_interface.h"

static KbacTest* Ktest = NULL;

void set_up_kbac_test(int* nn, int* qq, double* aa, double* mafUpper, double* xdatIn, double* ydatIn, double* mafIn, int* xcol, int* ylen) {
    Ktest = new KbacTest(nn, qq, aa, mafUpper, xdatIn, ydatIn, mafIn, xcol, ylen);
    return;
}

void do_kbac_test(double* pvalue, int* sided) {
    Ktest->calcKbacP(pvalue, sided);
    return;
}

void clear_kbac_test() {
    delete Ktest;
    return;
}
