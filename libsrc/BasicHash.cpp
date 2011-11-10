////////////////////////////////////////////////////////////////////// 
// libsrc/BasicHash.cpp 
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
 
#include "BasicHash.h"
#include "Error.h"

#include <stdio.h>

BasicHash::BasicHash(int startsize)
   {
   count = 0;
   size  = startsize;
   mask  = startsize - 1;

   // In this implementation, the size of hash tables must be a power of two
   if (startsize & mask)
      error("BasicHash: Hash table size must be a power of two.\n");

   objects = new void * [size];
   keys    = new unsigned int [size];

   for (unsigned int i = 0; i < size; i++)
      { objects[i] = NULL; }
   };

BasicHash::~BasicHash()
   {
   delete [] objects;
   delete [] keys;
   }

void BasicHash::Clear()
   {
//   printf("Clearing...\n");

   count = 0;

   if (size > 16)
      SetSize(16);

   for (unsigned int i = 0; i < size; i++)
      objects[i] = NULL;
   }

void BasicHash::SetSize(int newsize)
   {
   int newmask = newsize - 1;

   void     ** newobjects = new void * [newsize];
   unsigned int * newkeys = new unsigned int [newsize];

   for (int i = 0; i < newsize; i++)
       { newobjects[i] = NULL; }

   if (count)
      for (unsigned int i = 0; i < size; i++)
         if (objects[i] != NULL)
            {
            unsigned int key = keys[i];
            unsigned int h   = key & newmask;

            while ( newobjects[h] != NULL && newkeys[h] != h)
               h = (h + 1) & newmask;

            newkeys[h] = key;
            newobjects[h] = objects[i];
            }

   delete [] objects;
   delete [] keys;

   objects = newobjects;
   keys = newkeys;
   size = newsize;
   mask = newmask;
   }

int BasicHash::Add(int key, void * object)
   {
   if (count * 2 > size)
      Grow();

   unsigned int h = Iterate(key);

   while ((objects[h] != NULL) && (objects[h] != object))
      h = ReIterate(key, h);

   if (objects[h] == NULL)
      {
//      printf("At position %d, inserted %x\n", h, key);
      keys[h] = key;
      count++;
      }

   objects[h] = object;

   return h;
   }

int BasicHash::Find(int key)
   {
   int h = Iterate(key);

   return objects[h] == NULL ? -1 : h;
   }

int BasicHash::Rehash(int key, int h)
   {
   h = ReIterate(key, h);

   return objects[h] == NULL ? -1 : h;
   }

void BasicHash::Delete(unsigned int index)
   {
   if (index >= size || objects[index] == NULL)
      return;

   objects[index] = NULL;
   count--;

   if (count * 8 < size && size > 32)
      Shrink();
   else
      {
      // rehash the next entries until we find empty slot
      index = (index + 1) & mask;

      while (objects[index] != NULL)
         {
         if ((keys[index] & mask) != index)
            {
            unsigned int h = Iterate(keys[index]);

            while ((objects[h] != NULL) && (objects[h] != objects[index]))
               h = ReIterate(keys[index], h);

            if (h != (unsigned int) index)
               {
               keys[h] = keys[index];
               objects[h] = objects[index];
               objects[index] = NULL;
               }
            }

         index = (index + 1) & mask;
         }
      }
   }
 
