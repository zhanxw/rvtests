////////////////////////////////////////////////////////////////////// 
// libsrc/StringMap.h 
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
 
#ifndef __STRINGMAP_H__
#define __STRINGMAP_H__

#include "StringBasics.h"

class StringMap
   {
   protected:
      ::String ** strings;
      void     ** objects;
      int         count, size;

   public:
      static int alloc;

      StringMap(int startsize = 0);
      virtual ~StringMap();

      void Grow(int newsize);
      void Clear();
      int  Length() const { return count; }

      void * Object(int i) const { return objects[i]; }
      void * Object(const ::String & key) const
         {
         int index = Find(key);
         return (index >= 0) ? objects[index] : NULL;
         }
      void * Object(const ::String & key, void * (*create_object)())
         { return objects[Find(key, create_object)]; }

      void SetObject(int i, void * object)
         { objects[i] = object; }
      void SetObject(const ::String & key, void * object)
         { Add(key, object); }

      int Add(const ::String & s, void * object = NULL);
      int Find(const ::String & s, void * (*create_object)() = NULL);
      int Find(const ::String & s) const;
      int FindStem(const ::String & stem) const;
      int FindFirstStem(const ::String & stem) const;

      StringMap & operator = (const StringMap & rhs);

      const ::String & operator [] (int i) const { return *(strings[i]); }
      ::String & operator [] (int i) { return *(strings[i]); }
      ::String & String(int i) { return *(strings[i]); }

      static void * CreateMap();

      void Delete(int index);
   };

class StringIntMap
   {
   protected:
      ::String ** strings;
      int       * integers;
      int         count, size;

   public:
      static int alloc;

      StringIntMap(int startsize = 0);
      virtual ~StringIntMap();

      void Grow(int newsize);
      void Clear();
      int  Length() const { return count; }

      int Integer(int i) const { return integers[i]; }
      int Integer(const ::String & key) const
         {
         int index = Find(key);
         return (index >= 0) ? (int) integers[index] : -1;
         }

      void SetInteger(int i, int value)
         { integers[i] = value; }
      void SetInteger(const ::String & key, int value)
         { Add(key, value); }

      int Add(const ::String & s, int i);
      int Find(const ::String & s, int defaultValue);
      int Find(const ::String & s) const;
      int FindStem(const ::String & stem) const;

      StringIntMap & operator = (const StringIntMap & rhs);

      const ::String & operator [] (int i) const { return *(strings[i]); }
      ::String & operator [] (int i) { return *(strings[i]); }
      ::String & String(int i) { return *(strings[i]); }

      static void * CreateMap();

      int IncrementCount(const ::String & key);
      int DecrementCount(const ::String & key);
      int GetCount(const ::String & key) const;
      int GetCount(int index) const { return integers[index]; }

      void Delete(int index);
   };

#endif

 
