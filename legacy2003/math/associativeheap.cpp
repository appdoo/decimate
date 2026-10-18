//=============================================================================
//
// AssociativeHeap : AssociativeHeap Implementation
//
// Elemente des Heap wissen über Ihre Position im Heap Bescheid
//
// $Author: Andor $
// $Date: 31.10.03 13:21 $
// $Revision: 1 $
// 
//=============================================================================

// local includes
#include "associativeheap.h"
#include "..\math\math3D.h"



// ----------------------------------------------------------------------------
// Erzeugung eines Heap Objektes
// das die angegebene Anzahl von Elementen
// speichern kann
// ----------------------------------------------------------------------------
AssociativeHeap::AssociativeHeap(int initialSize)
{
   minHeap = NULL;
   Init(initialSize);
}


// ----------------------------------------------------------------------------
// Set capacity of Heap
// ----------------------------------------------------------------------------
void AssociativeHeap::Init(int initialSize)
{
   if (minHeap) delete[] minHeap;
   //
   maxElements = initialSize;
   heapSize    = 0;
   minHeap     = new Heapable*[maxElements];
}


// ----------------------------------------------------------------------------
// Freigabe eines Heap Objektes
// ----------------------------------------------------------------------------
AssociativeHeap::~AssociativeHeap()
{
   if (minHeap) delete[] minHeap;

   maxElements = 0;
   heapSize    = 0;
   minHeap     = NULL;
}


// ----------------------------------------------------------------------------
// fügt dem Heap einen Pointer
// auf ein HeapElement hinzu
// ----------------------------------------------------------------------------
void AssociativeHeap::Push(Heapable *heapable)
{
   if (!heapSize) {
      minHeap[heapSize++] = heapable;
      heapable->heapPos   = 0;
   }
   else {
      // vertex paar als letztes element in den heap einfügen
      int k             = heapSize++;
      minHeap[k]        = heapable;
      heapable->heapPos = k;

      UpHeap(k);
   }
}

// ----------------------------------------------------------------------------
// entfernt das oberste Element
// aus dem Heap
// ----------------------------------------------------------------------------
Heapable *AssociativeHeap::Pop()
{
   return Remove(0);
}


// ----------------------------------------------------------------------------
// entfernt ein Element aus dem
// Heap
// ----------------------------------------------------------------------------
Heapable *AssociativeHeap::Remove(Heapable *element)
{
   return Remove(element->heapPos);
}


// ----------------------------------------------------------------------------
// entfernt ein Element aus dem
// Heap
// ----------------------------------------------------------------------------
Heapable *AssociativeHeap::Remove(int position)
{
   // wertebereich prüfen
   if (heapSize < 1)
      return NULL;

   // element merken
   Heapable *heapable = minHeap[position];

   // oberstes paar ersetzen
   minHeap[position] = minHeap[--heapSize];
   minHeap[position]->heapPos = position;

   // element 'nach unten fallenlassen'
   DownHeap(position);

   return heapable;
}


// ----------------------------------------------------------------------------
// Funktion dient dazu um die Position
// eines Elements im Heap zu korrigieren
// ----------------------------------------------------------------------------
void AssociativeHeap::Repair(Heapable *element)
{
   Repair(element->heapPos);
}


// ----------------------------------------------------------------------------
// Funktion dient dazu um die Position
// eines Elements im Heap zu korrigieren
// ----------------------------------------------------------------------------
void AssociativeHeap::Repair(int pos)
{
   // vorgaenger groesser
   if ((pos>0)  &&  (minHeap[pos]->heapValue < minHeap[(pos - 1) >> 1]->heapValue)) 
   {
      // element nach oben wandern lassen
      UpHeap(pos);
   }
   else DownHeap(pos);
}


// ----------------------------------------------------------------------------
// vertauschen des angegebenen Elements mit dem
// Vorgänger, solange dieser kleiner ist
// ----------------------------------------------------------------------------
void AssociativeHeap::UpHeap(int k)
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


// ----------------------------------------------------------------------------
// vertauschen des angegebenen Elements mit
// dem jeweils kleinsten Nachfolger
// ----------------------------------------------------------------------------
void AssociativeHeap::DownHeap(int pos)
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


