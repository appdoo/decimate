//=============================================================================
//
// ErrorSystem Interface
//
// Grundgerüst für die Fehlerbestimmung Für den EdgeCollapse -
// Algorithmus
//
// $Author: Andor $
// $Date: 11.12.03 0:16 $
// $Revision: 3 $
// 
//=============================================================================
#ifndef IERRORSYSTEM_H
#define IERRORSYSTEM_H

// Local includes
#include "ivertexpaircollection.h"

// begin namespace
namespace decimate {

// ============================================================================
// Interface for ErrorSystem Classes
//
//
class IErrorSystem
{
public:
   virtual void Initialize(int numPoints, Vec3 *points, int numFaces, Face *faces) = 0;
   virtual void Destroy() = 0;

   // notifies errorsystem about an edge collapse operation
   virtual void Update(int p1, int p2) = 0;

   // calculates the optimal point and error for an edge collapse operation
   virtual void Optimize(int p1, int p2, bool boundary, double &error, Vec3 &p) = 0;
};

} // end namespace
#endif