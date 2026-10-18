//=============================================================================
//
// CChartBoundary : Class Definition
//
// Buchführung über Boundary - Eigenschaften von Kanten eines
// Polygon - Meshes
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 10 $
// 
//=============================================================================

// local includes
#include "cchartboundary.h"

namespace chartification {
// ============================================================================
//
// CChartBoundary class implementation
//
//

// ----------------------------------------------------------------------------
// constructs CChartBoundary
// ----------------------------------------------------------------------------
CChartBoundary::CChartBoundary(int numPoints, int numHalfedges, Halfedge *halfedges)
{
   // reference data
   this->numHalfedges = numHalfedges;
   this->halfedges    = halfedges;
   this->numPoints    = numPoints;

   // create point degrees and boundary flags
   degrees    = new int[numPoints];
   boundaries = new bool[numHalfedges];

   // reset to an initial state
   Reset();
}

// ----------------------------------------------------------------------------
// destructs CChartBoundary
// ----------------------------------------------------------------------------
CChartBoundary::~CChartBoundary()
{
   if (degrees)    delete[] degrees;
   if (boundaries) delete[] boundaries;
}

// ----------------------------------------------------------------------------
// erases a (half)edge from boundaries
// ----------------------------------------------------------------------------
void CChartBoundary::Erase(Halfedge *h)
{
   // Kantenbasiert, Lower mappt auf 
   // die kleinere Halbkante
   h = Lower(h);

   if (boundaries[h - halfedges]) {
      // löschen der Außenkanten - Flags
      boundaries[h - halfedges] = false;
      if (h->sym) boundaries[h->sym - halfedges] = false;

      // Herabsetzen der Grade
      if (degrees[h->end] > 0)             degrees[h->end]--;
      if (degrees[h->next->next->end] > 0) degrees[h->next->next->end]--;

      // Rekursives Löschen falls Grad(Endpunkt) < 2
      if (degrees[h->end] < 2) {
         Halfedge *invalid = SearchBoundary(h);
         if (invalid) Erase(invalid);
      }

      // Rekursives Löschen falls Grad(Startpunkt) < 2
      if (degrees[h->next->next->end] < 2) {
         Halfedge *invalid = SearchBoundary(h->next->next);
         if (invalid) Erase(invalid);
      }
   }
}


// ----------------------------------------------------------------------------
// tells the programmer, wheter this (half)edge is boundary, or not
// ----------------------------------------------------------------------------
bool CChartBoundary::IsBoundary(Halfedge *h)
{
   return boundaries[Lower(h) - halfedges];
}

// ----------------------------------------------------------------------------
// Find a boundary edge around end of start edge
// ----------------------------------------------------------------------------
Halfedge *CChartBoundary::SearchBoundary(Halfedge *start)
{
   // Search first direction
   Halfedge *loop = start;

   do {
      if (IsBoundary(loop))       return loop;
      if (IsBoundary(loop->next)) return loop->next;
      //
      loop = loop->next->sym;
      //
   }  while (loop && loop!=start);
   // check exit condition
   if (loop == start) return NULL;

   // Search second direction
   loop = start;

    do {
      if (IsBoundary(loop))                   return loop;
      if (loop->sym && IsBoundary(loop->sym)) return loop->sym;
      //
      // Search second direction
      if (loop->sym) loop = loop->sym->next->next;
      else           loop = NULL;

   } while (loop && loop!=start);
   // nothing found
   return NULL;
}

// ----------------------------------------------------------------------------
// returns the halfedge with lower face index of h, h->sym
// ----------------------------------------------------------------------------
Halfedge *CChartBoundary::Lower(Halfedge *h)
{
   return (h->sym && h->sym->face < h->face) ? h->sym : h;
}


// ----------------------------------------------------------------------------
// resets the boundary structure
// ----------------------------------------------------------------------------
void CChartBoundary::Reset()
{
   bool *visited = new bool[numHalfedges];

   // initialize flags
   memset(visited,   false, sizeof(bool) * numHalfedges);
   memset(boundaries, true, sizeof(bool) * numHalfedges);
   memset(degrees,       0, sizeof(int)  * numPoints);

   // initialize data
   for (int i=0; i<numHalfedges; i++) {
      // if halfedge wasn't visited process it
      if (!visited[i]) {
         Halfedge *h = Lower(&halfedges[i]);
         // mark edge
         visited[h - halfedges] = true;
         if (h->sym) visited[h->sym - halfedges] = true;

         // increase degree of both endpoints
         degrees[h->end]++;
         degrees[h->next->next->end]++;
      }
   }
   delete[] visited;
}

// ----------------------------------------------------------------------------
// Random access to set Boundary Flag
// ----------------------------------------------------------------------------
void CChartBoundary::SetBoundaryFlag(Halfedge *h, bool boundary)
{
   // konsistenzpruefung
   if (boundary == IsBoundary(h)) return;

   // Lower half
   h = Lower(h);

   // set boundary flag
   boundaries[h - halfedges] = boundary;
   if (h->sym) boundaries[h->sym - halfedges] = boundary;

   // change degrees
   if (boundary) {
      degrees[h->end]++;
      degrees[h->next->next->end]++;
   }
   else {
      degrees[h->end]--;
      degrees[h->next->next->end]--;
   }
}


} // end namespace chartification

