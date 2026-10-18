//=============================================================================
//
// CParametrizedChart : Class Definition
//
// $Author: Andor $
// $Date: 5.12.03 15:02 $
// $Revision: 2 $
// 
//=============================================================================
#ifndef CPARAMETRIZEDCHART_H
#define CPARAMETRIZEDCHART_H

#include <set>

// Local includes
#include "..\math\math3D.h"
#include "..\math\halfedge.h"

namespace parametrization {

class CParametrizedChart
{
   // values for chart parametrization
   Vec3 *points;
   Face *faces;

   // vertices for texture coordinates
   Vertex *vertices;
   int    numVertices;

   // boundary point set
   int_set boundaryPoints;
   
public:
   // creates a ChartSet
   CParametrizedChart(Vec3 *points, Face *faces);

   // deletes ChartSet and all contained Parametrized Charts
   ~CParametrizedChart();

   // adds a face to this representation
   void Add(Face *face);

   // adds a boundary information
   void Add(Edge *edge);

   // access to Vertex structure
   Vertex *Vertices();

   // number of vertices generated
   int VertexCount();
};

} // end namespace
#endif