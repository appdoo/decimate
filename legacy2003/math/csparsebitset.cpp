//=============================================================================
//
// CSparseBitSet : Class Implementation
//
// BitSet Klasse
//
// $Author: Andor $
// $Date: 3.12.03 20:28 $
// $Revision: 1 $
// 
//=============================================================================

// Local includes
#include "csparsebitset.h"


// ----------------------------------------------------------------------------
// class constructor
// ----------------------------------------------------------------------------
CSparseBitSet::CSparseBitSet()
{
   Init(5);
}


// ----------------------------------------------------------------------------
// class destructor
// ----------------------------------------------------------------------------
CSparseBitSet::~CSparseBitSet()
{
   delete[] bits;
   delete[] offs;
}


// ----------------------------------------------------------------------------
// initializes the bitset arrays
// ----------------------------------------------------------------------------
void CSparseBitSet::Init(int maxLength)
{
   this->maxLength = maxLength;
   this->bits      = new int[maxLength];
   this->offs      = new int[maxLength];
   this->length    = 0;
}


// ----------------------------------------------------------------------------
// Inserts a new entry
// ----------------------------------------------------------------------------
int CSparseBitSet::Insert(int bit)
{
   if (length >= maxLength) {
      maxLength   += 5;
      int *tmpBits = new int[maxLength];
      int *tmpOffs = new int[maxLength];

      memcpy(tmpBits, bits, sizeof(int) * length);
      memcpy(tmpOffs, offs, sizeof(int) * length);

      delete[] bits;
      delete[] offs;

      bits = tmpBits;
      offs = tmpOffs;
   }
   offs[length] = bit >> 5;
   bits[length] = 0;
   //
   return length++;
}

// ----------------------------------------------------------------------------
// Finds the BitSet offset
// ----------------------------------------------------------------------------
int CSparseBitSet::Find(int offset)
{
   for (int i=0; i<length; i++) {
      if (offs[i] == offset) return i;
   }
   return -1;
}


// ----------------------------------------------------------------------------
// Sets a bit into the BitSet
// ----------------------------------------------------------------------------
void CSparseBitSet::insert(int bit)
{
   int offset = Find(bit >> 5);
   // insert new entry
   if (offset < 0) offset = Insert(bit);
   bits[offset] |= 1 << (bit & 0x1F);
}

// ----------------------------------------------------------------------------
// kills a bit
// ----------------------------------------------------------------------------
void CSparseBitSet::erase(int bit)
{
   int offset = Find(bit >> 5);
   if (offset > 0) {
      bits[offset] &= !(1 << bit & 0x1F);
   }
}



