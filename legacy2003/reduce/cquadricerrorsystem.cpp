//=============================================================================
//
// CQuadricErrorSystem : ErrorSystem Implementation
//
// ErrorSystem using quadric error metric
//
// $Author: Andor $
// $Date: 11.12.03 0:16 $
// $Revision: 4 $
// 
//=============================================================================

// Local includes
#include "cquadricerrorsystem.h"

// begin namespace
namespace decimate {

// ============================================================================
// CQuadricErrorSystem Class Implementation
//
//


// ----------------------------------------------------------------------------
// creates an instance
// ----------------------------------------------------------------------------
CQuadricErrorSystem::CQuadricErrorSystem()
{
}


// ----------------------------------------------------------------------------
// deletes an instance
// ----------------------------------------------------------------------------
CQuadricErrorSystem::~CQuadricErrorSystem()
{
   Destroy();
}

// ----------------------------------------------------------------------------
// initializes instance
// ----------------------------------------------------------------------------
void CQuadricErrorSystem::Initialize(int numPoints, Vec3 *points, int numFaces, Face *faces)
{
   pointRef = points;

   errorMetrics = new QError[numPoints];
   QError error;

   for (int i=0; i<numPoints; i++)
      errorMetrics[i].Init();

   for (int i=0; i<numFaces; i++) {
      // Now calculate the plane equation for each plane
      // and add its errorMetric to each point belonging to the face
      Face *f = &faces[i];

      error.SetPlane(points[f->a], points[f->b], points[f->c]);

      errorMetrics[f->a].Add(error);
      errorMetrics[f->b].Add(error);
      errorMetrics[f->c].Add(error);
   }
}


// ----------------------------------------------------------------------------
// cleanup
// ----------------------------------------------------------------------------
void CQuadricErrorSystem::Destroy()
{
   if (errorMetrics) delete[] errorMetrics;
   errorMetrics = NULL;
}


// ----------------------------------------------------------------------------
// notifies errorsystem about a edge collapse operation
// ----------------------------------------------------------------------------
void CQuadricErrorSystem::Update(int p1, int p2)
{
   errorMetrics[p1].Add(errorMetrics[p2]);
}


// ----------------------------------------------------------------------------
// calculates the optimal point and error for an edge collapse operation
// ----------------------------------------------------------------------------
void CQuadricErrorSystem::Optimize(int p1, int p2, bool boundary, double &error, Vec3 &p)
{
   QError qerror(errorMetrics[p1], errorMetrics[p2]);

   if (boundary || !qerror.Optimize(p))
   {
      p = (pointRef[p1] + pointRef[p2]) * 0.5f;
   }
   error = qerror.GetError(p);
}


} // end namespace
