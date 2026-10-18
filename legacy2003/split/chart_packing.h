/******************************************************************************
 *
 * Chart packing Algorithm
 *
 *
 *****************************************************************************/
#ifndef CHART_PACKING_H
#define CHART_PACKING_H

// Local Includes
#include "../math/math3D.h"
#include "../math/vec.h"

namespace chartpacking {

class Horizon
{
public:
   // integer horizon with given resolution
   int *horizon;
   // resolution of the horizon
   const int resolution;
   // scalefactor for incoming charts
   float scale;

   /***************************************************************************
    * create an Atlas with the given resolution
    **************************************************************************/
   Horizon(int newresolution) : resolution(newresolution)
   {
      horizon = new int[resolution];
      //
      reset();
   }

   /***************************************************************************
    * delete Atlas
    **************************************************************************/
   ~Horizon()
   {
      if (horizon) delete[] horizon;
   }

   /***************************************************************************
    * reset horizon to zero
    **************************************************************************/
   void reset()
   {
      memset(horizon, 0, sizeof(int) * resolution);
   }

   /***************************************************************************
    * updates the horizon
    **************************************************************************/
   void update(int xoffs, int yoffs, int size, int *chart)
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

   /***************************************************************************
    * returns the discrete area of this atlas
    **************************************************************************/
   int area()
   {
      return resolution * resolution;
   }

   /***************************************************************************
    * calculates the discrete area between horizon and chart
    **************************************************************************/
   int area(int xoffs, int yoffs, int size, int *chart)
   {
      // discrete area
      int area = 0;

      for (int i=0; i<size; i++) {
         // check bounds
         if ((xoffs + i) < resolution && (horizon[xoffs + i] < (yoffs + chart[i]))) {
            area += yoffs + chart[i] - horizon[xoffs + i];
         }
         // stop condition reached
         else return -1;
      }
      return area;
   }
};


/******************************************************************************
 * packs a chart into the texture - atlas
 *
 * @param horizon,     a discrete horizon
 * @param numEdges,    number of boundary edges
 * @param edges,       the boundary edges
 * @param numTextures, number of texture coordinates
 * @param textures,    texture coordinates
 *
 *****************************************************************************/
void chartpack(Horizon &horizon, int numEdges, Edge *edges, int numTextures, Tex2 *textures);


/******************************************************************************
 * include chartpacker class
 *****************************************************************************/
#include "chartpacker.h"


} // end namespace

#endif
