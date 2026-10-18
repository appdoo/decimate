// ============================================================================
//
// Class HalfEdgeMesh : Definition
//
// Erzeugen einer Halfedge Datenstruktur
//
// $Author: Andor $
// $Date: 10.12.03 17:14 $
// $Revision: 4 $
//
// ============================================================================
#ifndef HALFEDGE_H
#define HALFEDGE_H

#include "math3D.h"
#include "cset.h"

// ----------------------------------------------------------------------------
// halfedge
// ----------------------------------------------------------------------------
typedef struct Halfedge
{
   // face index, the winged edge belongs to
   int face;
   // index of vertex at the end of directed edge
   int end;
   // points to symmetric halfedge
   Halfedge *sym;
   // points to next counterclockwise edge
   Halfedge *next;
}
Halfedge;

// halfedges
typedef std::vector<Halfedge *> halfedge_array;
typedef std::set<Halfedge *> halfedge_set;
typedef std::set<Halfedge *>::iterator halfedge_iterator;

// ----------------------------------------------------------------------------
// utility functions
// ----------------------------------------------------------------------------

// first outgoing halfedge, starting with an incoming halfedge
inline Halfedge *first_outgoing(Halfedge *h)
{
   return h->next;
}

// next outgoing halfedge starting with an already outgoing halfegde
inline Halfedge *next_outgoing(Halfedge *h)
{
   return h->sym ? h->sym->next : NULL;
}

// ----------------------------------------------------------------------------
// Halfedge based mesh representation
// ----------------------------------------------------------------------------
class HalfEdgeMesh
{
   // number of edges
   int numHalfEdges;
   // halfedge pointer
   Halfedge *halfedges;
   // number of points
   int numPoints;
   // coordinate pointer
   Vec3 *points;

   // -------------------------------------------------------------------------
   // constructs an empty halfedge
   // -------------------------------------------------------------------------
   HalfEdgeMesh();

   // -------------------------------------------------------------------------
   // deletes halfedge structure
   // -------------------------------------------------------------------------
   ~HalfEdgeMesh();

   // -------------------------------------------------------------------------
   // searches for the symmetric wedge
   // -------------------------------------------------------------------------
   static Halfedge *Symmetric(Halfedge *h, Halfedge *halfedges, CSet<int> *connect);

   // -------------------------------------------------------------------------
   // writes a face into a halfedge buffer
   // -------------------------------------------------------------------------
   static void WriteFace(int face_index, Face *faces, Halfedge *halfedges);
public:
   // -------------------------------------------------------------------------
   // creates a halfedge mesh representation from a standard mesh description
   // -------------------------------------------------------------------------
   static void Create(int numPoints, int numFaces, Face *faces,
                      int *outNumHalfEdges, Halfedge **outHalfEdges);
};
#endif
