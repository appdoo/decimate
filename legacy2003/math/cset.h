//=============================================================================
//
// CSet : Class Definition
// 
// Array Template
//
// $Author: Andor $
// $Date: 11.12.03 0:16 $
// $Revision: 3 $
// 
//=============================================================================
#ifndef CSET_H
#define CSET_H

#include <stdio.h>

template <class T> class CSet
{
private:
   // Elemente
   T    *set;

   int  size;
   int  max;

public:
   typedef T* iterator;
   iterator begin() { return &set[0]; }
   iterator end()   { return &set[size]; }
   iterator rbegin() { return &set[size - 1]; }
   iterator rend()   { return &set[-1]; }

   // -------------------------------------------------------------------------
   // CArray constructor
   // -------------------------------------------------------------------------
   CSet(int initialSize = 0)
   {
      max  = initialSize;
      size = 0;
      set  = NULL;
   }

   // -------------------------------------------------------------------------
   // Aufraeumen vor dem Freigeben des Objekts
   // -------------------------------------------------------------------------
   ~CSet()
   {
      if (set) delete[] set;

      max  = 0;
      size = 0;
      set  = NULL;
   }

   // -------------------------------------------------------------------------
   // liefert das Element an der gegebenen Stelle
   // -------------------------------------------------------------------------
   T &operator[](int index)
   {
      if (index < 0 || index >= size) {
         printf("Set:[] 'Index - Array out of Bounds!!!\n'");
         index = MIN(index, size -1);
         index = MAX(index, 0);
      }
      return array[index];
   }

   // -------------------------------------------------------------------------
   // löscht das Element x vom Set
   // -------------------------------------------------------------------------
   int Erase(T &x)
   {
      // position suchen
      int pos = Contains(x);

      // check if element exists
      if (pos != -1) {
         // remove it
         for (int i=pos; i<size - 1; i++) set[i] = set[i + 1];
         size--;
      }
      return pos;
   }

   // -------------------------------------------------------------------------
   // Tauscht das Element in dem Set
   // -------------------------------------------------------------------------
   int Swap(T &o, T &n)
   {
      // position suchen
      int pos = Contains(o);

      if (pos != -1) {
         while (pos > 0 && set[pos - 1] > n) {
            set[pos] = set[pos - 1];
            pos--;
         }
         while (pos < size - 1 && set[pos + 1] < n) {
            set[pos + 1] = set[pos];
            pos++;
         }
         set[pos] = n;
      }
      return pos;
   }

   // -------------------------------------------------------------------------
   // löscht alle Elemente vom Set
   // -------------------------------------------------------------------------
   void Clear()
   {
      if (set) delete[] set;
      set  = NULL;
      size = 0;
      max  = 0;
   }

   // -------------------------------------------------------------------------
   // search for element, return index or -1
   // -------------------------------------------------------------------------
   int Contains(T &x)
   {
      if (!set) return -1;

      // binäre suche
      int anf = 0;
      int end = (size - 1);
      int m   = (size - 1) >> 1;

      while (anf <= end && x!=set[m]) {
         if (x<set[m]) end = m - 1;
         else          anf = m + 1;

         m = (anf + end) >> 1;
      }
      return (x == set[m]) ? m : -1;
   }

   // -------------------------------------------------------------------------
   // Fügt ein Element in das Set ein
   // -------------------------------------------------------------------------
   int Insert(T &x)
   {
      // check for set
      int pos = Contains(x);
      // set does only contain an object 1 time
      if (pos != -1) return pos;

      // reallocate
      if (size >= max) Reallocate(max + 10);

      // sortieren durch direktes einfügen
      // bestehendes set ist bereits sortiert
      int ug = 0; 
      int og = size - 1;
      int j  = 0;
      // bin suche
      while (ug<=og) {
         j = (ug + og) >> 1;
         if (x<set[j]) og = j - 1;
         else          ug = j + 1;
      }
      // verschieben
      for (j = size; j>ug; j--) set[j] = set[j - 1];
      // einfügen
      set[ug] = x;
      size++;
      // position zurückliefern
      return ug;
   }

   // -------------------------------------------------------------------------
   // Inserts all elements from cset into own set
   // -------------------------------------------------------------------------
   void Merge(CSet<T> &cset)
   {
      if (cset.size == 0)
         return;
      // new temporary set
      T  *nset  = new T[size + cset.size];

      // merge
      int i=0; int j=0; int n=0;

      while (i<size && j<cset.size) {
         if      (set[i] < cset.set[j]) nset[n++] = set[i++];
         else if (cset.set[j] < set[i]) nset[n++] = cset.set[j++];
         else {
            nset[n++] = set[i++];
            j++;
         }
      }
      // copy rest
      while (i < size)      nset[n++] = set[i++];
      while (j < cset.size) nset[n++] = cset.set[j++];

      if (n>=max) Reallocate(n);

      memcpy(set, nset, sizeof(T) * n);
      size = n;

      delete[] nset;
   }

   // -------------------------------------------------------------------------
   // shrink / grow of set
   // -------------------------------------------------------------------------
   void Reallocate(int minsize)
   {
      minsize = MAX(minsize, size);
      // array vergrössern
      max = minsize;

      // reallocate
      T *nset = new T[max];
      if (set) memcpy(nset, set, sizeof(T) * size);
      if (set) delete[] set;
      set = nset;
   }

   // -------------------------------------------------------------------------
   // liefert die Anzahl der Elemente des Heaps
   // -------------------------------------------------------------------------
   int Size()
   {
      return size;
   }
};


#endif