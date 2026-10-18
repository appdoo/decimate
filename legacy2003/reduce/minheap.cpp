/**
 * Class MinHeap
 * 
 * Implementierung eines Minimum-Heap
 * (Das kleinste Element ist immer an erster Stelle)
 *
 * Elemente des Heap wissen über Ihre Position im Heap Bescheid
 *
 * @author: Andor Greissl
 * @date:   22.03.2003
 */
#include "minheap.h"


#define MAX(x,y) ((x)>(y) ? (x):(y))
#define MIN(x,y) ((x)<(y) ? (x):(y))

template <class T> inline void swap(T &a, T&b)
{
   T t;
   t = a;
   a = b;
   b = t;
}


/**
 * Erzeugung eines Heap Objektes
 * das die angegebene Anzahl von Elementen
 * speichern kann
 */
MinHeap::MinHeap(int initialSize)
{
   maxElements = initialSize;
   heapSize    = 0;
   minHeap     = new HeapElement*[maxElements];
}


/**
 * Freigabe eines Heap Objektes
 */
MinHeap::~MinHeap()
{
   if (minHeap) delete[] minHeap;

   maxElements = 0;
   heapSize    = 0;
   minHeap     = NULL;
}




/**
 * fügt dem Heap einen Pointer
 * auf ein HeapElement hinzu
 */
void MinHeap::add(HeapElement *element)
{
   if (!heapSize) {
      minHeap[heapSize++] = element;
      element->heapPos    = 0;
   }
   else {
      // vertex paar als letztes element in den heap einfügen
      int k            = heapSize++;
      minHeap[k]       = element;
      element->heapPos = k;

      upHeap(k);
   }
}

/**
 * entfernt ein Element aus dem
 * Heap
 */
HeapElement *MinHeap::remove(HeapElement *element)
{
   return remove(element->heapPos);
}


/**
 * entfernt ein Element aus dem
 * Heap
 */
HeapElement *MinHeap::remove(int position)
{
   // wertebereich prüfen
   if (heapSize < 1)
      return NULL;

   // element merken
   HeapElement *element = minHeap[position];

   // oberstes paar ersetzen
   minHeap[position] = minHeap[--heapSize];
   minHeap[position]->heapPos = position;

   // element 'nach unten fallenlassen'
   downHeap(position);

   return element;
}


/**
 * Funktion dient dazu um die Position
 * eines Elements im Heap zu korrigieren
 */
void MinHeap::repair(HeapElement *element)
{
   repair(element->heapPos);
}


/**
 * Funktion dient dazu um die Position
 * eines Elements im Heap zu korrigieren
 */
void MinHeap::repair(int pos)
{
   // vorgaenger groesser
   if ((pos>0)  &&  (minHeap[pos]->heapValue < minHeap[(pos - 1) >> 1]->heapValue)) 
   {
      // element nach oben wandern lassen
      upHeap(pos);
   }
   else downHeap(pos);
}


/**
 * vertauschen des angegebenen Elements mit dem
 * Vorgänger, solange dieser kleiner ist
 */
void MinHeap::upHeap(int k)
{
   // vorgaenger bestimmen
   int v = (k - 1) >> 1;

   while ((k) && (minHeap[k]->heapValue < minHeap[v]->heapValue)) {

      swap(minHeap[k], minHeap[v]);
      minHeap[k]->heapPos = k;
      minHeap[v]->heapPos = v;
      //
      k = v;
      v = (k - 1) >> 1;
   }
}


/**
 * vertauschen des angegebenen Elements mit
 * dem jeweils kleinsten Nachfolger
 */
void MinHeap::downHeap(int pos)
{
   for(;;) {
      int childPos = pos + pos + 1;

      // 1. Abbruchkriterium
      if (childPos >= heapSize)
         break;

      // Ast mit dem kleinerem Fehler wählen
      if ((childPos + 1 < heapSize) && (minHeap[childPos + 1]->heapValue < minHeap[childPos]->heapValue))
      {
         childPos++;
      }
      
      // 2. Abbruchkriterium
      if (minHeap[pos]->heapValue < minHeap[childPos]->heapValue)
         break;

      // vertauschen, und index berichtigen
      swap(minHeap[pos], minHeap[childPos]);
      minHeap[pos]->heapPos      = pos;
      minHeap[childPos]->heapPos = childPos;
      //
      pos = childPos;
   }
}


