////////////////////////////////////////////////////////////////////// 
// libsrc/LongArray.cpp 
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
 
#include "LongArray.h"
#include "Hash.h"
#include "Sort.h"

#include <string.h>

int LongArray::alloc = 4;

LongArray::LongArray(int start_size)
   {
   count = start_size;
   size = (count + alloc) / alloc * alloc;
   items = new longint [size];
   }

LongArray::LongArray(LongArray & source)
   {
   count = source.count;
   size = source.size;
   items = new longint [size];

   for (int i = 0; i < count; i++)
      items[i] = source.items[i];
   }

LongArray::~LongArray()
   {
   delete [] items;
   }

void LongArray::Grow(int new_size)
   {
   if (new_size > size)
      {
      if ((new_size >> 1) >= size)
         size = (new_size + alloc) / alloc * alloc;
      else
         {
         size = alloc;
         while (size <= new_size)
            size *= 2;
         }

      longint * new_items = new longint [size];
      for (int i = 0; i < count; i++)
         new_items[i] = items[i];
      delete [] items;
      items = new_items;
      }
   }

int LongArray::Append(longint value)
   {
   Grow(count + 1);
   items[count++] = value;
   return count;
   }

void LongArray::Set(longint value)
   {
   for (int i = 0; i < count; i++)
      items[i] = value;
   }

int LongArray::Delete(int index)
   {
   count--;
   if (count - index)
      memmove(items + index, items + index + 1, sizeof(longint) * (count - index));
   return count;
   }

void LongArray::InsertAt(int index, longint value)
   {
   Grow(count + 1);
   memmove(items + index + 1, items + index, sizeof(longint) * (count - index));
   items[index] = value;
   count++;
   }

LongArray & LongArray::operator = (const LongArray & rhs)
   {
   Grow(rhs.count);
   count = rhs.count;
   for (int i = 0; i < count; i++)
      items[i] = rhs.items[i];
   return *this;
   }

int LongArray::Find(longint value) const
   {
   for (int i = 0; i < count; i++)
      if (value == items[i])
         return i;
   return -1;
   }

void LongArray::Zero()
   {
   for (int i = 0; i < count; i++)
      items[i] = 0;
   }

void LongArray::Reverse()
   {
   for (int i = 0, j = count - 1; i < j; i++, j--)
      Swap(i, j);
   }

bool LongArray::operator == (const LongArray & rhs) const
   {
   if (count != rhs.count)
      return false;

   for (int i = 0; i < rhs.count; i++)
      if (items[i] != rhs.items[i])
         return false;

   return true;
   }

bool LongArray::operator != (const LongArray & rhs) const
   {
   return !(*this == rhs);
   }

int LongArray::Hash(int initval)
   {
   return hash((unsigned char *) items, sizeof(longint) * count, initval);
   }
 
