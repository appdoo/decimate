//=============================================================================
//
// CHorizon : Class Implementation
//
// Diskreter Horizont für den Atlaspacker
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 2 $
// 
//=============================================================================

// Local Includes
#include "../math/math3D.h"
#include "../math/vec.h"
#include "chorizon.h"

namespace chartification {

// ----------------------------------------------------------------------------
// create an Atlas with the given resolution
// ----------------------------------------------------------------------------
CHorizon::CHorizon(int newresolution) : resolution(newresolution)
{
   horizon = new int[resolution];
   //
   Reset();
}

// ----------------------------------------------------------------------------
// delete Atlas
// ----------------------------------------------------------------------------
CHorizon::~CHorizon()
{
   if (horizon) delete[] horizon;
}

// ----------------------------------------------------------------------------
// access to integer horizon
// ----------------------------------------------------------------------------
int CHorizon::operator[](int index)
{
   return horizon[index];
}

// ----------------------------------------------------------------------------
// reset horizon to zero
// ----------------------------------------------------------------------------
void CHorizon::Reset()
{
   memset(horizon, 0, sizeof(int) * resolution);
}

// ----------------------------------------------------------------------------
// updates the horizon
// ----------------------------------------------------------------------------
void CHorizon::Update(int xoffs, int yoffs, int size, int *chart)
{
   // update horizon to new values
   for (int i=0; i<size; i++) {
      // check bounds
      if (xoffs + i < resolution) {
         // update horizon with a little border
         horizon[xoffs + i] = chart[i] + yoffs + 2;
      }
   }
}

// ----------------------------------------------------------------------------
// returns the discrete area of this atlas
// ----------------------------------------------------------------------------
int CHorizon::Area()
{
   return resolution * resolution;
}

// ----------------------------------------------------------------------------
// calculates the discrete area between horizon and chart
// ----------------------------------------------------------------------------
int CHorizon::Area(int x, int y, int size, int *chart)
{
   // Fläche zwischen der beiden Horizonten
   int area = 0;

   for (int i=0; i<size; i++) {
      // Prüfen ob sich Horizont des Charts komplett innerhalb
      // des Texturatlas befindet
      if ((x + i) < resolution && (horizon[x + i] < (y + chart[i]))) {
         area += y + chart[i] - horizon[x + i];
      }
      // Abbruch
      else return -1;
   }
   return area;
}

} // end namespace