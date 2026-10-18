//=============================================================================
//
// CArray : Class Definition
// 
// Array Template
//
// $Author: Andor $
// $Date: 6.12.03 21:22 $
// $Revision: 2 $
// 
//=============================================================================
#ifndef CARRAY_H
#define CARRAY_H

#include <stdio.h>

template <class T> class CArray
{
private:
   // Elemente
   T    *array;

   int  size;
   int  max;

public:

   // -------------------------------------------------------------------------
   // CArray constructor
   // -------------------------------------------------------------------------
   CArray(int initialSize = 100)
   {
      max   = initialSize;
      size  = 0;
      array = new T[max];
   }

   // -------------------------------------------------------------------------
   // Aufraeumen vor dem Freigeben des Objekts
   // -------------------------------------------------------------------------
   ~CArray()
   {
      if (array) delete[] array;

      max   = 0;
      size  = 0;
      array = NULL;
   }

   // -------------------------------------------------------------------------
   // löscht alle Elemente vom Array
   // -------------------------------------------------------------------------
   void Clear()
   {
      size = 0;
   }

   // -------------------------------------------------------------------------
   // liefert das Element an der gegebenen Stelle
   // -------------------------------------------------------------------------
   T &operator[](int index)
   {
      if (index < 0 || index >= size) {
         printf("Index - Array out of Bounds!!!\n");
         index = MIN(index, size -1);
         index = MAX(index, 0);
      }
      return array[index];
   }

   // -------------------------------------------------------------------------
   // liefert die Anzahl der Elemente des Heaps
   // -------------------------------------------------------------------------
   int Size()
   {
      return size;
   }

   // -------------------------------------------------------------------------
   // Fügt dem Array ein Element hinzu und liefert
   // seinen Index
   // -------------------------------------------------------------------------
   int Add(T &object)
   {
      if (size >= max) {
         // array vergrössern
         max = max + 1000;

         // reallocate
         T *tmp   = new T[max];
         memcpy(tmp, array, sizeof(T) * size);
         delete[] array;
         array = tmp;
      }
      // object an letzter stelle einfügen
      array[size++] = object;
      // position zurückliefern
      return size - 1;
   }

};


#endif