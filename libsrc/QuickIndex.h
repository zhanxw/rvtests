////////////////////////////////////////////////////////////////////// 
// libsrc/QuickIndex.h 
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
 
#ifndef __QUICKINDEX_H__
#define __QUICKINDEX_H__

#include "MathVector.h"
#include "StringArray.h"
#include "StringHash.h"
#include "IntArray.h"
#include "StringMap.h"

class QuickIndex : public IntArray
   {
   public:
      QuickIndex();
      QuickIndex(const IntArray & source_data)
         { Index(source_data); }
      QuickIndex(const StringArray & source_data)
         { Index(source_data); }
      QuickIndex(const Vector & source_data)
         { Index(source_data); }

      void Index(const IntArray & source_data);
      void Index(const StringArray & source_data);
      void Index(const Vector & source_data);
      void IndexCounts(const StringIntMap & source_data);
      void IndexCounts(const StringIntHash & source_data);

   private:
      const void * source;
      int    datatype;

      bool IsBefore(int i, int j);
      void Sort();
   };

#endif

 
