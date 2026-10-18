//=============================================================================
//
// CVertexPairCollection Class Definition
//
//
// $Author: Andor $
// $Date: 26.08.03 15:48 $
// $Revision: $
// 
//=============================================================================
#ifndef CVERTEXPAIRCOLLECTION_H
#define CVERTEXPAIRCOLLECTION_H

// local includes
#include "ivertexpaircollection.h"
#include "../math/halfedge.h"

// begin namespace
namespace decimate {


// ============================================================================
// CVertexPairIterator Class Definition
//
//
class CVertexPairCollection : public IVertexPairCollection
{
   VertexPair *pairs;
   int numPairs;

   // discard unused space
   void reallocate();
   // add one entry to pairs
   void add(Halfedge *h, bool boundary);

public:
   CVertexPairCollection();
   virtual ~CVertexPairCollection();

   virtual void Initialize(int numPoints, int numFaces, Face *faces);
   virtual void Destroy();

   // element count
   virtual  int size();

   // random access
   virtual VertexPair *operator[](int index);
   virtual int find(VertexPair *vertexPair);

   // iterator
   virtual iterator begin();
   virtual iterator end();
};

} // end namespace
#endif