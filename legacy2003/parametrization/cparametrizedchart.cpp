//=============================================================================
//
// CParametrizedChart : Class Implementaion
//
// $Author: Andor $
// $Date: 5.12.03 15:02 $
// $Revision: 2 $
// 
//=============================================================================

// C++ includes
#include <set>

// Local includes
#include "..\math\math3D.h"
#include "..\math\halfedge.h"
#include "cparametrizedchart.h"

namespace parametrization {

// ----------------------------------------------------------------------------
// creates a ChartSet
// ----------------------------------------------------------------------------
CParametrizedChart::CParametrizedChart(Vec3 *points, Face *faces)
{
}

// ----------------------------------------------------------------------------
// deletes ChartSet and all contained Parametrized Charts
// ----------------------------------------------------------------------------
CParametrizedChart::~CParametrizedChart()
{
}


// ----------------------------------------------------------------------------
// adds a face to this representation
// ----------------------------------------------------------------------------
void CParametrizedChart::Add(Face *face)
{
}


// ----------------------------------------------------------------------------
// adds a boundary information
// ----------------------------------------------------------------------------
void CParametrizedChart::Add(Edge *edge)
{
   // only take endpoint as a boundary point
   boundaryPoints.insert(edge->b);
}


// ----------------------------------------------------------------------------
// access to Vertex structure
// ----------------------------------------------------------------------------
Vertex *CParametrizedChart::Vertices()
{
   return vertices;
}


// ----------------------------------------------------------------------------
// number of vertices generated
// ----------------------------------------------------------------------------
int CParametrizedChart::VertexCount()
{
   return numVertices;
}

} // end namespace
