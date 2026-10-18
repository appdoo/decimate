//=============================================================================
//
// Grundgerüst für einen Halfedge basierten Edge-Collapse Algorithmus
//
// $Author: Andor $
// $Date: 11.12.03 23:00 $
// $Revision: 1 $
// 
//=============================================================================

// Local Includes
#include "chedgecollapse.h"

// begin namespace
namespace decimate {

// ============================================================================
// CHalfedgeCollapse Class Implementation
//
//
// ----------------------------------------------------------------------------
// creates an chedgecollapse object
// ----------------------------------------------------------------------------
CHalfedgeCollapse::CHalfedgeCollapse(int numPoints, const Vec3 *points, int numFaces, const Face *faces, IErrorSystem *errorSystem)
{
   this->numPoints = numPoints;
   this->points    = new Vec3[numPoints];
   this->numFaces  = numFaces;
   this->faces     = new Face[numFaces];

   this->numValidPoints       = numPoints;
   this->errorSystem          = errorSystem;

   // copy data
   memcpy(this->points, points, sizeof(Vec3) * numPoints);
   memcpy(this->faces,   faces, sizeof(Face) * numFaces);

   Initialize();
}

// ----------------------------------------------------------------------------
// deletes an edgecollapse object
// ----------------------------------------------------------------------------
CHalfedgeCollapse::~CHalfedgeCollapse()
{
   errorSystem->Destroy();

   delete[] collapseEntries;
   delete[] halfedges;
   delete[] points;
   delete[] faces;
   delete[] invalidFaces;
}

// ----------------------------------------------------------------------------
// decimates the data structure
// ----------------------------------------------------------------------------
bool CHalfedgeCollapse::Decimate()
{
   return true;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void CHalfedgeCollapse::DecimateToPoints(int numDesiredPoints, int *numCoord, Vec3 **coords, int *numFaces, Face **faces)
{
}

// ----------------------------------------------------------------------------
// Setup data structures
// ----------------------------------------------------------------------------
void CHalfedgeCollapse::Initialize()
{
   // Fehlersystem initialisieren
   errorSystem->Initialize(numPoints, points, numFaces, faces);

   // Halbeckenstruktur generieren
   HalfEdgeMesh::Create(numPoints, numFaces, faces, &numHalfedges, &halfedges);

   // flags fuer ungueltige Flaechen
   invalidFaces = new bool[numFaces];
   memset(invalidFaces, false, sizeof(bool) * numFaces);

   // Heap initialisieren
   collapseEntries = new CollapseEntry[numHalfedges];
   bool *visited   = new bool[numHalfedges];
   memset(visited, false, sizeof(bool) * numHalfedges);

   for (int i=0; i<numHalfedges; i++) {
      // check for edge
      if (visited[i]) continue;

      // set flag
      visited[i] = true;
      if (halfedges[i].sym) visited[halfedges[i].sym - halfedges] = true;

      // Fehler und neuen Punkt berechnen
      errorSystem->Optimize(
         halfedges[i].next->next->end,
         halfedges[i].end,
         halfedges[i].sym == NULL,
         collapseEntries[i].heapValue,
         collapseEntries[i].p
         );

      // Paar auf heap legen
      heap.Push(&collapseEntries[i]);
   }
   delete[] visited;
}

// ----------------------------------------------------------------------------
// Update
// ----------------------------------------------------------------------------
void CHalfedgeCollapse::UpdateError(int index)
{
}

// ----------------------------------------------------------------------------
// Punkt - Referenzen ersetzen
// ----------------------------------------------------------------------------
void CHalfedgeCollapse::ReplacePoint(Halfedge *start)
{
   // p1 der halbkante
   int p1 = start->next->next->end;

      // Search first direction
   Halfedge *loop = start;

   do {
      // replace point
      loop->end = p1;
      // next halfedge
      loop = loop->next->sym;
      //
   }  while (loop && loop!=start);
   // check exit condition
   if (loop == start) return;

   // Search second direction
   loop = start;

    do {
      // replace point
      loop->end = p1;
      // Search second direction
      if (loop->sym) loop = loop->sym->next->next;
      else           loop = NULL;

   } while (loop && loop!=start);
}

// ----------------------------------------------------------------------------
// Collapse Operation for a halfedge
// ----------------------------------------------------------------------------
void CHalfedgeCollapse::Collapse(Halfedge *h)
{
   // Punkt ersetzen
   ReplacePoint(h);

   // angrenzende Fläche ausklinken
   if (h->next->sym)       h->next->sym->sym = h->next->next->sym;
   if (h->next->next->sym) h->next->next->sym->sym = h->next->sym;

   invalidFaces[h->face] = true;

   // auf symmetrische Fläche prüfen
   if (h->sym) {
      // Benachbarte Halbkante behandeln
      h = h->sym;

      // angrenzende Fläche ausklinken
      if (h->next->sym)       h->next->sym->sym = h->next->next->sym;
      if (h->next->next->sym) h->next->next->sym->sym = h->next->sym;

      invalidFaces[h->face] = true;
   }
}


} // end namespace
