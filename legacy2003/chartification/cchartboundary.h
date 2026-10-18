//=============================================================================
//
// CChartBoundary : Class Definition
//
// Buchführung über Boundary - Eigenschaften von Kanten eines
// Polygon - Meshes
//
// $Author: Andor $
// $Date: 11.12.03 0:16 $
// $Revision: 7 $
// 
//=============================================================================
#ifndef CCHARTBOUNDARY_H
#define CCHARTBOUNDARY_H

// local includes
#include "../math/math3D.h"
#include "../math/halfedge.h"

namespace chartification {


// ============================================================================
//
// CChartBoundary class defintion
//
//
class CChartBoundary
{
   int          numPoints;
   int          numHalfedges;
   Halfedge     *halfedges;
   int          *degrees;
   bool         *boundaries;

public:
   // constructs CChartBoundary
   CChartBoundary(int numPoints, int numHalfedges, Halfedge *halfedges);

   // destructs CChartBoundary
   ~CChartBoundary();

   // erases a (half)edge from boundaries
   void Erase(Halfedge *h);

   // tells the programmer, wheter this (half)edge is boundary, or not
   bool IsBoundary(Halfedge *h);

   // Random access to set Boundary Flag
   void SetBoundaryFlag(Halfedge *h, bool boundary);

   // resets the boundary structure
   void Reset();

private:
   // removes non extremal halfedges
   void RemoveNonExtremal(Halfedge *start);

   // Find a boundary edge around end of start edge
   Halfedge *SearchBoundary(Halfedge *start);

   // returns the halfedge with lower face index of h, h->sym
   Halfedge *Lower(Halfedge *h);

}; // end class CChartBoundary

} // end namespace chartification
#endif
