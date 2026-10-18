/******************************************************************************
 *
 * Chart packing algorithm
 *
 * Basisklasse
 *
 *****************************************************************************/
#ifndef CHART_PACKER_H
#define CHART_PACKER_H

#include "../math/math3d.h"
#include <vector>

class ChartPacker
{
   std::vector<float> sizes;
   std::vector<int> numTextures;
   std::vector<Tex2 *> textures;
   std::vector<int> numPaths;
   std::vector<Edge *> paths;

public:
   /***************************************************************************
    *
    * Constructs a class for chart packing
    *
    **************************************************************************/
   ChartPacker()
   {
   }

   /***************************************************************************
    * 
    **************************************************************************/
   ~ChartPacker()
   {
   }

   /***************************************************************************
    * adds a reference to texture coordinates
    **************************************************************************/
   void addRef(float size, int numTexture, Tex2 *texture, int numPath, Edge *path)
   {
      int index = findIndex(size);
      // step to position
      std::vector<float>::iterator  sit = sizes.begin();
      std::vector<int>::iterator   ntit = numTextures.begin();
      std::vector<Tex2 *>::iterator tit = textures.begin();
      std::vector<int>::iterator   npit = numPaths.begin();
      std::vector<Edge *>::iterator pit = paths.begin();

      for (int i=0; i<index; i++) {
         sit++;
         ntit++;
         tit++;
         npit++;
         pit++;
      }
      // insert values
      sizes.insert(sit, size);
      numTextures.insert(ntit, numTexture);
      textures.insert(tit, texture);
      numPaths.insert(npit, numPath);
      paths.insert(pit, path);
   }


   /***************************************************************************
    * optimizes the chart packing and returns the packing ratio
    **************************************************************************/
   float optimize(float surface_area, int discreet_size)
   {
      // start with 99 percent packing ratio
      float ratio = 99.0f;
      // try to pack charts - if no success, diminish ratio
      while (!pack_charts(surface_area, discreet_size, ratio)) ratio-=1.0f;

      return ratio;
   }

private:
   /***************************************************************************
    * finds index in list
    **************************************************************************/
   int findIndex(float size)
   {
      // set start index
      int index = 0;

      for (std::vector<float>::iterator it=sizes.begin(); it!=sizes.end(); it++) {
         // iindex found
         if (*it > size) return index;
         // increment index
         index++;
      }
      return index;
   }


   /***************************************************************************
    * try to pack the charts
    **************************************************************************/
   bool pack_charts(float surface_area, int discreet_size, float ratio)
   {
      // create a discreet horizon
      Horizon horizon(discreet_size);
      horizon.scale = ratio / (sqrtf(surface_area) * horizon.resolution);


      return false;
   }
};





#endif



