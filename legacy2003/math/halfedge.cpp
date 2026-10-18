// ============================================================================
//
// Class HalfEdgeMesh : Implementation
//
// Erzeugen einer Halfedge Datenstruktur
//
// $Author: Andor $
// $Date: 10.12.03 17:14 $
// $Revision: 4 $
//
// ============================================================================
#include "halfedge.h"


// ----------------------------------------------------------------------------
// constructs an empty halfedge
// ----------------------------------------------------------------------------
HalfEdgeMesh::HalfEdgeMesh()
{
   numHalfEdges = 0;
   numPoints = 0;
   halfedges = NULL;
   points = NULL;
}

// ----------------------------------------------------------------------------
// deletes halfedge structure
// ----------------------------------------------------------------------------
HalfEdgeMesh::~HalfEdgeMesh()
{
   if (halfedges) delete[] halfedges;
   if (points)    delete[] points;
}

// ----------------------------------------------------------------------------
// searches for the symmetric wedge
// ----------------------------------------------------------------------------
Halfedge *HalfEdgeMesh::Symmetric(Halfedge *h, Halfedge *halfedges, CSet<int> *connect)
{
   // get the point indices of the edge
   int a = h->next->next->end;
   int b = h->end;
   // get corresponding faces for a
   CSet<int> &set = connect[a];
   // loop through the faces
   for (CSet<int>::iterator it = set.begin(); it!=set.end(); it++) {
      int edgeIndex = *it * 3;
      // get possible corresponding wedge-set
      for (int i=0; i<3; i++) {
         Halfedge *sym = &halfedges[edgeIndex++];
         if ((sym->end==a) && (sym->next->next->end==b)) {
            // return corresponding halfedge
            return sym;
         }
      }
   }
   return NULL;
}

// ----------------------------------------------------------------------------
// writes a face into a halfedge buffer
// ----------------------------------------------------------------------------
void HalfEdgeMesh::WriteFace(int face_index, Face *faces, Halfedge *halfedges)
{
   // every face has 3 halfedges
   int index = face_index * 3;
   // allocate the three wedges for the face
   Halfedge *ac = &halfedges[index++];
   Halfedge *cb = &halfedges[index++];
   Halfedge *ba = &halfedges[index++];
   // enter face info
   ac->face = face_index;
   cb->face = face_index;
   ba->face = face_index;
   // enter end point connectivity
   ac->end = faces[face_index].c;
   cb->end = faces[face_index].b;
   ba->end = faces[face_index].a;
   // enter wedge 'next' connectivity
   ac->next = cb;
   cb->next = ba;
   ba->next = ac;
}

// ----------------------------------------------------------------------------
// creates a halfedge mesh representation from a standard mesh description
// ----------------------------------------------------------------------------
void HalfEdgeMesh::Create(int numPoints, int numFaces, Face *faces,
                          int *outNumHalfEdges, Halfedge **outHalfEdges)
{
   // allocate memory for halfedges, and initialize to zero
   int numHalfEdges = numFaces * 3;
   Halfedge *halfedges = new Halfedge[numHalfEdges];
   memset(halfedges, 0, sizeof(Halfedge) * numHalfEdges);

   // create connectivity information between adjacent faces
   // an edge of a face has one neigbour or no neighbour
   // 1. create point to face connectivity information
   CSet<int> *connect = new CSet<int>[numPoints];
   // add all faces to create connectivity
   for (int i=0; i<numFaces; i++) {
      connect[faces[i].a].Insert(i);
      connect[faces[i].b].Insert(i);
      connect[faces[i].c].Insert(i);
   }

   // create winged edges without symmetric
   // halfedge connectivity
   // add wegdes in a counter - clockwise order
   for (int i=0; i<numFaces; i++) {
      WriteFace(i, faces, halfedges);
   }
   // create the symmetric connectivity information
   for (int i=0; i < numHalfEdges; i++) {
      if (!halfedges[i].sym) {
         halfedges[i].sym = Symmetric(&halfedges[i], halfedges, connect);
         if (halfedges[i].sym) halfedges[i].sym->sym = &halfedges[i];
      }
   }
   // cleanup
   delete[] connect;
   // set outgoing parameters
   *outNumHalfEdges = numHalfEdges;
   *outHalfEdges = halfedges;
}


