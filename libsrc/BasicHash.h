////////////////////////////////////////////////////////////////////// 
// libsrc/BasicHash.h 
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
 
#ifndef __BASICHASH_H__
#define __BASICHASH_H__

#include <stdlib.h>

class BasicHash
   {
   protected:
      void          ** objects;
      unsigned int      * keys;
      unsigned int count, size;
      unsigned int        mask;

   public:
      BasicHash(int startsize = 32);
      virtual ~BasicHash();

      void Grow()    { SetSize(size * 2); }
      void Shrink()  { SetSize(size / 2); }

      void SetSize(int newsize);

      void Clear();

      int  Capacity() const { return size; }
      int  Entries() const  { return count; }

      void * Object(int i) const { return objects[i]; }

      void SetObject(int i, void * object)
         { objects[i] = object; }

      int Add    (int key, void * object = NULL);
      int Find   (int key);
      int Rehash (int key, int h);

      BasicHash & operator = (const BasicHash & rhs);

      void * operator [] (int i) const { return objects[i]; }

      void Delete(unsigned int index);

      bool SlotInUse(int index) { return objects[index] != NULL; }

   private:
      unsigned int Iterate(unsigned int key) const
         {
         unsigned int h = key & mask;

         while (objects[h] != NULL && keys[h] != key)
            h = (h + 1) & mask;

         return h;
         }

      unsigned int ReIterate(unsigned int key, unsigned int h) const
         {
         h = (h + 1) & mask;

         while (objects[h] != NULL && keys[h] != key)
            h = (h + 1) & mask;

         return h;
         }
   };

#endif
 
