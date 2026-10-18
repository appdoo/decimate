//=============================================================================
//
// Grundgerüst für einen Halfedge basierten Edge-Collapse Algorithmus
//
// $Author: Andor $
// $Date: 11.12.03 23:00 $
// $Revision: 1 $
// 
//=============================================================================

#ifndef CHALFEDGECOLLAPSE_H
#define CHALFEDGECOLLAPSE_H


// Local Includes
#include "ierrorsystem.h"
#include "..\math\halfedge.h"
#include "..\math\cset.h"

// begin namespace
namespace decimate {

struct CollapseEntry : public Heapable
{
public:
   Vec3 p;
};



// ============================================================================
// CHalfedgeCollapse Class Definition
//
//
class CHalfedgeCollapse
{
   Halfedge       *halfedges;
   int             numHalfedges;

   IErrorSystem   *errorSystem;
   AssociativeHeap heap;

   CollapseEntry  *collapseEntries;

   bool           *invalidFaces;

   int             numFaces;
   int             numPoints;

   int             numValidPoints;

   // Flächen und Punkte des Polygon-Meshes
   Face           *faces;
   Vec3           *points;


public:
   // creates an chedgecollapse object
   CHalfedgeCollapse(int numPoints, const Vec3 *points, int numFaces, const Face *faces, IErrorSystem *errorSystem);

   // deletes an edgecollapse object
   ~CHalfedgeCollapse();

   // decimates the data structure
   bool Decimate();
   void DecimateToPoints(int numDesiredPoints, int *numCoord, Vec3 **coords, int *numFaces, Face **faces);


private:
   // Setup data structures
   void Initialize();

   // Update
   void UpdateError(int index);

   // Punkt - Referenzen ersetzen
   void ReplacePoint(Halfedge *start);

   // Collapse Operation for a halfedge
   void Collapse(Halfedge *h);
};


} // end namespace
#endif



