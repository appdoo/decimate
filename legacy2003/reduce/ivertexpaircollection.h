//=============================================================================
//
// VertexPairCollection
//
// $Author: Andor $
// $Date: 24.11.03 11:10 $
// $Revision: 2 $
// 
//=============================================================================
#ifndef IVERTEXPAIRCOLLECTION_H
#define IVERTEXPAIRCOLLECTION_H

// local includes
#include "../math/associativeheap.h"
#include "../math/math3D.h"

// begin namespace
namespace decimate {

// define a VertexPair structure for internal communication
typedef struct VertexPair : public Heapable
{
   int   p1;
   int   p2;
   Vec3  p;
   bool  boundary;
}
VertexPair;

// ============================================================================
// Interface for VertexPairCollection Classes
//
//
class IVertexPairCollection
{
public:
   typedef VertexPair *iterator;

   virtual void Initialize(int numPoints, int numFaces, Face *faces) = 0;
   virtual void Destroy() = 0;

   // element count
   virtual  int size() = 0;

   // random access
   virtual VertexPair *operator[](int index) = 0;
   virtual int find(VertexPair *vertexPair) = 0;

   // iterator
   virtual iterator begin() = 0;
   virtual iterator end() = 0;
};

} // end namespace
#endif