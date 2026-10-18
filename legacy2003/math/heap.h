/**
 * Class Heap
 * 
 * Implementierung eines Heap
 * (Das kleinste Element ist immer an erster Stelle)
 *
 * @author: Andor Greissl
 * @date:   22.03.2003
 */
#ifndef HEAP_H
#define HEAP_H

#include <stdio.h>

/******************************************************************************
 *
 * Heap Klasse
 *
 * Der Heap arbeitet nur mir Pointern auf die eigentlichen
 * Heap - Elemente
 *
 ******************************************************************************/

template <class T, class V> class Heap
{
   template <class X> inline void swap(X &a, X &b)
   {
      X t;
      t = a;
      a = b;
      b = t;
   }

private:
   // Heap zum sortieren und einfügen
   T *heap;
   V *values;
   int         heapSize;
   int         maxElements;

public:
	/**
	 * Erzeugung eines Heap Objektes
	 * das die angegebene Anzahl von Elementen
    * speichern kann
	 */
	Heap(int initialSize = 100)
   {
      maxElements = initialSize;
      heapSize    = 0;
      heap        = new T[maxElements];
      values      = new V[maxElements];
   }

   /**
    * Aufraeumen vor dem Freigeben des Objekts
    */
   ~Heap()
   {
      if (heap) delete[] heap;
      if (values) delete[] values;

      maxElements = 0;
      heapSize    = 0;
      heap        = NULL;
      values      = NULL;
   }

   /***************************************************************************
    * löscht alle Elemente vom Heap
    **************************************************************************/
   void clear()
   {
      heapSize = 0;
   }

   // liefert das Element an der gegebenen Stelle
   void get(int index, T *object, V *value)
   {
      // copy element to outgoing parameters
      *object = heap[index];
      *value  = values[index];
   }

   // repariert den heap
   void repair()
   {
      int size = heapSize;
      heapSize = 0;

      for (int i=0; i<size; i++)
      {
         push(heap[i], values[i]);
      }
   }

   // Funktion dient dazu um die Position
   // eines Elements im Heap zu korrigieren
   void repair(int pos)
   {
      // vorgaenger groesser
      if ((pos>0)  &&  (values[pos] < values[(pos - 1) >> 1])) 
      {
         // element nach oben wandern lassen
         upHeap(pos);
      }
      else downHeap(pos);
   }

   // liefert die Anzahl der Elemente des Heaps
   int size()
   {
      return heapSize;
   }

   // setzt das Element an der gegebenen Stelle,
   // zerstört die Heap-Eigenschaft
   void set(int index, T object, V value)
   {
      // copy element to outgoing parameters
      heap[index]   = object;
      values[index] = value;
   }

   /**
    * fügt dem Heap einen Pointer
    * auf ein HeapElement hinzu
    */
   void push(T element, V value)
   {
      if (!heapSize) {
         heap[heapSize] = element;
         values[heapSize++] = value;
      }
      else {
         if (heapSize >= maxElements) {
            // reallocate heap
            maxElements = maxElements + 1000;
            T *tmp_heap   = new T[maxElements];
            V *tmp_values = new V[maxElements];
            memcpy(tmp_heap, heap, sizeof(T) * heapSize);
            memcpy(tmp_values, values, sizeof(V) * heapSize);
            // cleanup
            delete[] heap;
            delete[] values;
            // set pointers
            heap = tmp_heap;
            values = tmp_values;
         }
         // vertex paar als letztes element in den heap einfügen
         int k     = heapSize++;
         heap[k]   = element;
         values[k] = value;

         upHeap(k);
      }
   }

   /**
    * entfernt ein Element aus dem
    * Heap
    */
   bool pop(T *object, V *value)
   {
      // check if we have any elements
      if (heapSize < 1)
         return false;
      // copy element to outgoing parameters
      *object = heap[0];
      *value  = values[0];
      // insert the last element in position 0
      heap[0] = heap[--heapSize];
      values[0] = values[heapSize];
      // element 'nach unten fallenlassen'
      downHeap(0);
      // success
      return true;
   }

private:
   /******************************************************************************
   * vertauschen des angegebenen Elements mit dem
   * Vorgänger, solange dieser kleiner ist
   *****************************************************************************/
   void upHeap(int k)
   {
      // calculate predecessor
      int v = (k - 1) >> 1;

      while ((k) && (values[k] < values[v])) {
         // change the two elements
         swap(heap[k], heap[v]);
         swap(values[k], values[v]);
         //
         k = v;
         v = (k - 1) >> 1;
      }
   }

   /******************************************************************************
   * vertauschen des angegebenen Elements mit
   * dem jeweils kleinsten Nachfolger
   *****************************************************************************/
   void downHeap(int pos)
   {
      for(;;) {
         int childPos = pos + pos + 1;

         // 1. Abbruchkriterium
         if (childPos >= heapSize)
            break;

         // Ast mit dem kleinerem Fehler wählen
         if ((childPos + 1 < heapSize) && (values[childPos + 1] < values[childPos]))
         {
            childPos++;
         }
         
         // 2. Abbruchkriterium
         if (values[pos] < values[childPos])
            break;

         // vertauschen, und index berichtigen
         swap(heap[pos], heap[childPos]);
         swap(values[pos], values[childPos]);
         //
         pos = childPos;
      }
   }
};


#endif

