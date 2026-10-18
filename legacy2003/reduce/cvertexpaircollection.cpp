//=============================================================================
//
// CVertexPairCollection Class Implementation
//
//
// $Author: Andor $
// $Date: 10.12.03 17:14 $
// $Revision: 3 $
// 
//=============================================================================

// local includes
#include "cvertexpaircollection.h"

// begin namespace
namespace decimate {


// ============================================================================
// CVertexPairIterator Implementation
//
//

// ----------------------------------------------------------------------------
// initialize
// ----------------------------------------------------------------------------
CVertexPairCollection::CVertexPairCollection()
{
   pairs = NULL;
   numPairs = 0;
}


// ----------------------------------------------------------------------------
// destruct
// ----------------------------------------------------------------------------
CVertexPairCollection::~CVertexPairCollection()
{
   Destroy();
}


// ----------------------------------------------------------------------------
// discard unused space
// ----------------------------------------------------------------------------
void CVertexPairCollection::reallocate()
{
   VertexPair *tmp = new VertexPair[numPairs];
   memcpy(tmp, pairs, sizeof(VertexPair) * numPairs);
   // cleanup
   delete[] pairs;
   pairs = tmp;
}


// ----------------------------------------------------------------------------
// adds one entry
// ----------------------------------------------------------------------------
void CVertexPairCollection::add(Halfedge *h, bool boundary)
{
   // nächstes Paar in der Liste holen
   VertexPair *pair = &pairs[numPairs++];

   pair->p1       = h->end;
   pair->p2       = h->next->next->end;
   pair->boundary = boundary;
}


// ----------------------------------------------------------------------------
// setup
// ----------------------------------------------------------------------------
void CVertexPairCollection::Initialize(int numPoints, int numFaces, Face *faces)
{
   // HalfEdge Representation erzeugen
   int numHalfEdges;
   Halfedge *halfedges;
   HalfEdgeMesh::Create(numPoints, numFaces, faces, &numHalfEdges, &halfedges);

   // temp. marker erzeugen
   bool *visited = new bool[numHalfEdges];
   // initialisieren des markers
   for (int i=0; i<numHalfEdges; i++) {
      visited[i] = false;
   }
   // maximale anzahl pairs erzeugen
   pairs = new VertexPair[numHalfEdges];
   numPairs = 0;

   // iterate over each halfedge
   for (i=0; i<numHalfEdges; i++) {
      // check if we've been already there
      if (visited[i])
         continue;
      // fetch halfedge
      Halfedge *h = &halfedges[i];
      // mark halfedge, symmetric halfedge
      visited[i] = true;
      if (h->sym) visited[h->sym - halfedges] = true;
      // add halfedge to possible contraction pairs
      add(h, h->sym == NULL);
   }
   delete[] halfedges;
   delete[] visited;

   reallocate();
}


// ----------------------------------------------------------------------------
// cleanup
// ----------------------------------------------------------------------------
void CVertexPairCollection::Destroy()
{
   if (pairs) delete[] pairs;
   pairs = NULL;
   numPairs = 0;
}


// ----------------------------------------------------------------------------
// number of vertex - pairs
// ----------------------------------------------------------------------------
int CVertexPairCollection::size()
{
   return numPairs;
}


// ----------------------------------------------------------------------------
// random access
// ----------------------------------------------------------------------------
VertexPair *CVertexPairCollection::operator[](int index)
{
   return &pairs[index];
}


// ----------------------------------------------------------------------------
// returns index of pair
// ----------------------------------------------------------------------------
int CVertexPairCollection::find(VertexPair *vertexPair)
{
   return vertexPair - pairs;
}


// ----------------------------------------------------------------------------
// iterator to first element
// ----------------------------------------------------------------------------
IVertexPairCollection::iterator CVertexPairCollection::begin()
{
   return pairs;
}


// ----------------------------------------------------------------------------
// iterator to last element + 1
// ----------------------------------------------------------------------------
IVertexPairCollection::iterator CVertexPairCollection::end()
{
   return &pairs[numPairs];
}

} // end namespace
