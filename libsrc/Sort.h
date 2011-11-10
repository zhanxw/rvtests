////////////////////////////////////////////////////////////////////// 
// libsrc/Sort.h 
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
 
#ifndef __SORT_H__
#define __SORT_H__

#include "Constant.h"

#include <stddef.h>

void QuickSort(void *base, size_t nelem, size_t width,
               int (*cmp)(const void *, const void *));

void QuickSort2(void *base, void * base2, size_t nelem, size_t width,
               int (*cmp)(const void *, const void *));

void * BinarySearch(const void *key, const void *base,
               size_t nelem, size_t width,
               int (*cmp)(const void *, const void *));

#endif
 
