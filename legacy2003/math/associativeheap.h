// ============================================================================
//
// AssociativeHeap : AssociativeHeap
//
// Elemente des Heap wissen über Ihre Position im Heap Bescheid
//
// $Author: Andor $
// $Date: 31.10.03 13:21 $
// $Revision: 1 $
// 
// ============================================================================
#ifndef ASSOCIATIVEHEAP_H
#define ASSOCIATIVEHEAP_H

#include <stdio.h>


// ============================================================================
//
// Element des Heaps
//
// ============================================================================
class Heapable
{
public:
   double  heapValue;
   int     heapPos;
};


// ============================================================================
//
// Associative Heap Class Definition
//
//
class AssociativeHeap
{
private:
   // Heap zum sortieren und einfügen
   Heapable **minHeap;
   int      heapSize;
   int      maxElements;

public:
   // construct
	AssociativeHeap(int initialSize = 100);

   // destruct
   ~AssociativeHeap();

   // Set capacity of Heap
   void Init(int initialSize);

   // adds a heap element
   void Push(Heapable *element);

   // removes topmost element
   Heapable *Pop();

   // removes element at given position
   Heapable *Remove(int position);

   // removes element
   Heapable *Remove(Heapable *element);

   // repairs heap at given position
   void Repair(int position);

   // repairs heap for given element
   void Repair(Heapable *element);

private:
   // moves elements up
   void UpHeap(int k);

   // moves elements down
   void DownHeap(int k);
};


#endif


