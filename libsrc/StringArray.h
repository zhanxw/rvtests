////////////////////////////////////////////////////////////////////// 
// libsrc/StringArray.h 
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
 
#ifndef __STRING_ARRAY_H__
#define __STRING_ARRAY_H__

#include "StringBasics.h"

class StringArray
   {
   protected:
      String ** strings;
      int size, count;

   public:
      static int alloc;
      static bool lazyMemoryManagement;

      StringArray(int startsize = 0);
      StringArray(StringArray & original);
      virtual ~StringArray();

      // Each line in a file is parsed into a separate array element
      //

      void Read(FILE * f);
      void Write(FILE * f);
      void WriteLine(FILE * f);
      void Read(const char * filename);
      void Write(const char * filename);
      void WriteLine(const char * filename);

#ifdef __ZLIB_AVAILABLE__
      void Read(IFILE & f);
#endif

      // Write all strings to the screen
      void Print();
      void PrintLine();

      void Grow(int newsize);
      void Clear();

      int Length() const { return count; }
      int Dimension(int newcount);
      int CharLength();

      String & operator [] (int i) { return *(strings[i]); }
      const String & operator [] (int i) const { return *(strings[i]); }

      // These functions divide a string into tokens and append these to the
      // array. Return value is the new array length
      //

      int AddColumns(const String & s, char ch = '\t');
      int AddColumns(const String & s, char ch, int maxColumns);
      int AddTokens(const String & s, char ch);
      int AddTokens(const String & s, const String & separators = " \t\r\n");

      int ReplaceColumns(const String & s, char ch = '\t')
         { Clear(); return AddColumns(s, ch); }
      int ReplaceTokens(const String & s, const String & separators = " \t\r\n")
         { Clear(); return AddTokens(s, separators); }

      // These functions add, insert or remove a single array element
      //

      int  Add(const String & s);
      void InsertAt(int position, const String & s);
      void Delete(int position);

      // These functions manipulate a string as a stack
      //

      String & Last() const;
      int      Push(const String & s) { return Add(s); }
      String   Pop();

      // Linear search (N/2 comparisons on average) for a single element
      // If searching is required, StringMaps are a better option
      //

      int Find(const String & s) const;
      int FastFind(const String & s) const;
      int SlowFind(const String & s) const;

      // Alphetically orders strings
      //
      void Sort();

      // Trims strings to remove whitespace
      void Trim();

      StringArray & operator = (const StringArray & rhs);

      bool operator == (const StringArray & rhs);
      bool operator != (const StringArray & rhs)
         { return !(*this == rhs); }

      void Swap(StringArray & s);

   private:
      static int ComparisonForSort(const void * a, const void * b);
   };

#endif

 
