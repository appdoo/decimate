//=============================================================================
//
// CChartSet : Class Definition
//
// $Author: Andor $
// $Date: 6.12.03 1:14 $
// $Revision: 2 $
// 
//=============================================================================
#ifndef CCHARTSET_H
#define CCHARTSET_H

#include <set>

// Local includes
#include "..\math\math3D.h"
#include "..\math\halfedge.h"
#include "cchart.h"

namespace chartification {

class CChartSet
{
   std::set<CChart *> charts;

   // values for chart creation
   Vec3     *points;
   Vec3     *normals;
   Halfedge *halfedges;

public:
   // define iterator
   typedef std::set<CChart *>::iterator iterator;
   
   // creates a ChartSet
   CChartSet(Vec3 *points, Vec3 *normals, Halfedge *halfedges);

   // deletes ChartSet and all contained Charts
   ~CChartSet();

   // returns the area of this chartset
   float Area(Vec3 *points);

   // creates a new CChart
   CChart *Create(int seed);

   // Removes Invalid Charts
   bool RemoveInvalidCharts();

   // CChart iterators
   iterator Begin();
   iterator End();

   // returns the number of CCharts contained
   int Size();
};

} // end namespace
#endif