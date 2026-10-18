//=============================================================================
//
// CChartSet : Class Implementation
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 4 $
// 
//=============================================================================
//
// Local includes
#include "cchartset.h"

namespace chartification {

// ----------------------------------------------------------------------------
// creates a ChartSet
// ----------------------------------------------------------------------------
CChartSet::CChartSet(Vec3 *points, Vec3 *normals, Halfedge *halfedges)
{
   this->points        = points;
   this->normals       = normals;
   this->halfedges     = halfedges;
}


// ----------------------------------------------------------------------------
// deletes ChartSet and all contained Charts
// ----------------------------------------------------------------------------
CChartSet::~CChartSet()
{
   // delete stored CChart instances
   for (iterator it=Begin(); it!=End(); it++) {
      delete (*it);
   }
}

// ----------------------------------------------------------------------------
// returns area of this chartset
// ----------------------------------------------------------------------------
float CChartSet::Area(Vec3 *points)
{
   // initialize area to zero
   float area = 0.0f;
   
   // sum over all charts
   for (iterator it=Begin(); it!=End(); it++) {
      // add chart area
      area = area + (*it)->Area(points);
   }
   return area;
}

// ----------------------------------------------------------------------------
// creates a new CChart
// ----------------------------------------------------------------------------
CChart *CChartSet::Create(int seed)
{
   CChart *chart = new CChart(seed, normals[seed]);
   //
   charts.insert(chart);
   //
   return chart;
}


// ----------------------------------------------------------------------------
// entfernt ungültige charts aus dem chartset
// charts sind ungültig, wenn sie weniger als 2 flächen beinhalten
// ----------------------------------------------------------------------------
bool CChartSet::RemoveInvalidCharts()
{
   CChart **invalidCharts = new CChart*[charts.size()];
   int    invalidChartCount = 0;

   // search for invalid charts
   for (iterator it=Begin(); it!=End(); it++) {
      if ((*it)->ClusterCount() < 2) invalidCharts[invalidChartCount++] = (*it);
   }

   // correct chart_set
   for (int i=0; i<invalidChartCount; i++) {
      // remove chart from chart_set
      charts.erase(invalidCharts[i]);
      delete invalidCharts[i];
   }
   delete[] invalidCharts;

   // return if there were invalid charts
   return invalidChartCount > 0;
}


// ----------------------------------------------------------------------------
// CChart iterators
// ----------------------------------------------------------------------------
CChartSet::iterator CChartSet::Begin()
{
   return charts.begin();
}


// ----------------------------------------------------------------------------
// CChart iterators
// ----------------------------------------------------------------------------
CChartSet::iterator CChartSet::End()
{
   return charts.end();
}


// ----------------------------------------------------------------------------
// returns the number of CCharts contained
// ----------------------------------------------------------------------------
int CChartSet::Size()
{
   return charts.size();
}

} // end namespace