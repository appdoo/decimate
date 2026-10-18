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
#ifndef MINHEAP_H
#define MINHEAP_H

#include <stdio.h>

/******************************************************************************
 *
 * Interface fuer ein HeapElement
 *
 ******************************************************************************/
class HeapElement
{
public:
   double  heapValue;
   int    heapPos;
   int    arrPos;
};


/******************************************************************************
 *
 * Minimum Heap Klasse
 *
 * Der Heap arbeitet nur mir Pointern auf die eigentlichen
 * Heap - Elemente
 *
 ******************************************************************************/

class MinHeap
{
private:
   // Heap zum sortieren und einfügen
   HeapElement **minHeap;
   int         heapSize;

   int         maxElements;

public:
	/**
	 * Erzeugung eines Heap Objektes
	 * das die angegebene Anzahl von Elementen
    * speichern kann
	 */
	MinHeap(int initialSize = 100);


   /**
    * Aufraeumen vor dem Freigeben des Objekts
    */
   ~MinHeap();


   /**
    * fügt dem Heap einen Pointer
    * auf ein HeapElement hinzu
    */
   void add(HeapElement *element);


   /**
    * entfernt ein Element aus dem
    * Heap
    */
   HeapElement *remove(int position);


   /**
    * entfernt ein Element aus dem
    * Heap
    */
   HeapElement *remove(HeapElement *element);


   /**
    * Funktion dient dazu um die Position
    * eines Elements im Heap zu korrigieren
    */
   void repair(int position);

   
   /**
    * Funktion dient dazu um die Position
    * eines Elements im Heap zu korrigieren
    */
   void repair(HeapElement *element);


private:
   /**
    * vertauschen des angegebenen Elements mit dem
    * Vorgänger, solange dieser kleiner ist
    */
   void upHeap(int k);


   /**
    * vertauschen des angegebenen Elements mit
    * dem jeweils kleinsten Nachfolger
    */
   void downHeap(int k);
};


#endif

