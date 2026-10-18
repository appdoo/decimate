//=============================================================================
//
// CHorizon : Class Definition
//
// Diskreter Horizont für den Atlaspacker
//
// $Author: Andor $
// $Date: 6.12.03 1:14 $
// $Revision: 1 $
// 
//=============================================================================
#ifndef CHORIZON_H
#define CHORIZON_H

// Local Includes
#include "../math/math3D.h"
#include "../math/vec.h"

namespace chartification {

class CHorizon
{
public:
   // integer horizon with given resolution
   int *horizon;
   // resolution of the horizon
   const int resolution;
   // scalefactor for incoming charts
   float scale;

   // create an Atlas with the given resolution
   CHorizon(int newresolution);
   
   // delete Atlas
   ~CHorizon();

   // access to integer horizon
   int operator[](int index);

   // reset horizon to zero
   void Reset();

   // updates the horizon
   void Update(int xoffs, int yoffs, int size, int *chart);

   // returns the discrete area of this atlas
   int Area();

   // calculates the discrete area between horizon and chart
   int Area(int xoffs, int yoffs, int size, int *chart);
};

} // end namespace

#endif
