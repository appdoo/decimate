//=============================================================================
//
// CParametrizedChartSet : Class Definition
//
// $Author: Andor $
// $Date: 4.12.03 21:29 $
// $Revision: 1 $
// 
//=============================================================================
#ifndef CPARAMETRIZEDCHARTSET_H
#define CPARAMETRIZEDCHARTSET_H

#include <set>

// Local includes
#include "..\math\math3D.h"
#include "..\math\halfedge.h"
#include "..\chartification\cchartset.h"
#include "..\chartification\cchart.h"
#include "cparametrizedchart.h"
#include "iparametrizer.h"

namespace parametrization {

class CParametrizedChartSet
{
   // references to parametrzed charts
   std::set<CParametrizedChart *> pcharts;

   // values for chart parametrization
   Vec3 *points;
   Face *faces;

public:
   // define iterator
   typedef std::set<CParametrizedChart *>::iterator iterator;
   
   // creates a ChartSet
   CParametrizedChartSet(IParametrizer parametrizer, Vec3 *points, Face *faces);

   // deletes ChartSet and all contained Parametrized Charts
   ~CParametrizedChartSet();

   // parametrization run
   void Parametrize(chartification::CChartSet *charts);
   void Parametrize(chartification::CChart *chart);

   // CChart iterators
   iterator Begin();
   iterator End();

   // returns the number of CCharts contained
   int Size();
};

} // end namespace
#endif