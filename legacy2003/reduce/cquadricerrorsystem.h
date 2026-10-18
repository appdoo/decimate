//=============================================================================
//
// CQuadricErrorSystem : ErrorSystem Implementation
//
// ErrorSystem using quadric error metric
//
// $Author: Andor $
// $Date: 11.12.03 0:16 $
// $Revision: 3 $
// 
//=============================================================================
#ifndef CQUADRICERRORSYSTEM_H
#define CQUADRICERRORSYSTEM_H

// Local includes
#include "ierrorsystem.h"
#include "qerror.h"

// begin namespace
namespace decimate {

// ============================================================================
// CQuadricErrorSystem Class Definition
//
//
class CQuadricErrorSystem : public IErrorSystem
{
   QError *errorMetrics;
   Vec3 *pointRef;

public:
   CQuadricErrorSystem();
   virtual ~CQuadricErrorSystem();

   virtual void Initialize(int numPoints, Vec3 *points, int numFaces, Face *faces);
   virtual void Destroy();

   // notifies errorsystem about a edge collapse operation
   virtual void Update(int p1, int p2);

   // calculates the optimal point and error for an edge collapse operation
   virtual void Optimize(int p1, int p2, bool boundary, double &error, Vec3 &p);
};

} // end namespace
#endif