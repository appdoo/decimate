//=============================================================================
//
// CSparseBitSet : Class Definition
//
// BitSet Klasse
//
// $Author: Andor $
// $Date: 3.12.03 20:28 $
// $Revision: 1 $
// 
//=============================================================================
#ifndef CSPARSEBITSET_H
#define CSPARSEBITSET_H

#include "math3D.h"

class CSparseBitSet
{
   int *bits;
   int *offs;
   int length;
   int maxLength;

   // initializes the bitset arrays
   void Init(int maxLength);

   // Finds the correct bit index
   int Find(int offset);

   // Inserts a new Entry
   int Insert(int bit);

public:
   // class constructor
   CSparseBitSet();

   // class destructor
   ~CSparseBitSet();

   // Sets a bit into the BitSet
   void insert(int bit);

   // kills a bit
   void erase(int bit);
};

#endif