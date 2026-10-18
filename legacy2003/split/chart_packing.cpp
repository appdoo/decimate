/******************************************************************************
 *
 * Chart packing Algorithm
 *
 *
 *****************************************************************************/

// local includes
#include "../math/math3D.h"
#include "chart_packing.h"

namespace chartpacking
{

   /******************************************************************************
   * discretize using a line algorithm
   *****************************************************************************/
   void discretize(int *upper, int *lower, int u1, int v1, int u2, int v2)
   {
      // render left to right
      if (u1 > u2) {
         swap (u1, u2);
         swap (v1, v2);
      }
      int du = u2 - u1 + 1;
      int dv = v2 - v1 + 1;

      for (int u=u1; u<=u2; u++) {
         // interpolate
         int v = v1 + (dv * (u - u1)) / du;
         // apply
         if (upper[u] < v) upper[u] = v;
         if (lower[u] > v) lower[u] = v;
      }
   }

   /******************************************************************************
   * returns the lower left of the chart
   *****************************************************************************/
   Tex2 minuv(int numTextures, Tex2 *textures)
   {
      Tex2 min;

      min.u = textures[0].u;
      min.v = textures[0].v;

      for (int i=1; i<numTextures; i++) {
         min.u = MIN(min.u, textures[i].u);
         min.v = MIN(min.v, textures[i].v);
      }
      return min;
   }

   /******************************************************************************
   * returns the upper right of the chart
   *****************************************************************************/
   Tex2 maxuv(int numTextures, Tex2 *textures)
   {
      Tex2 max;

      max.u = textures[0].u;
      max.v = textures[0].v;

      for (int i=1; i<numTextures; i++) {
         max.u = MAX(max.u, textures[i].u);
         max.v = MAX(max.v, textures[i].v);
      }
      return max;
   }

   /******************************************************************************
   * returns the upper right of the chart
   *****************************************************************************/
   Tex2 miduv(int numTextures, Tex2 *textures)
   {
      Tex2 uvmid;

      Tex2 uvmin = minuv(numTextures, textures);
      Tex2 uvmax = maxuv(numTextures, textures);

      uvmid.u = (uvmax.u - uvmin.u) / 2;
      uvmid.v = (uvmax.v - uvmin.v) / 2;

      return uvmid;
   }

   /******************************************************************************
   * rotate srcCoord by given angle around midpoint and copy
   * result into dstCoord
   *****************************************************************************/
   void rotate(Tex2 *dstCoord, Tex2 *srcCoord, int numCoord, float angle, Tex2 &mid)
   {
      // convert into rad
      angle *= (float)M_PI / 180.0f;
      float sina = sin(angle);
      float cosa = cos(angle);

      for (int i=0; i<numCoord; i++) {
         // load u, v pair
         float u = srcCoord[i].u - mid.u;
         float v = srcCoord[i].v - mid.v;

         // store rotated version to dst
         dstCoord[i].u = mid.u + u * cosa - v * sina;
         dstCoord[i].v = mid.v + u * sina + v * cosa;
      }
   }

   /******************************************************************************
   * returns the oprimal chart position in discrete atlas coordinates
   *****************************************************************************/
   void minarea(int &x, int &y, int &minArea, Horizon &horizon, int width, int height, int *lower)
   {
      // create return value
      x  = 0;
      y  = 0;

      // set wasted area to maximum
      minArea = horizon.area();

      for (int i=0; i<horizon.resolution - width; i++) {
         for (int j=horizon.horizon[i]; j<horizon.resolution - height; j++) {
            // calculate area
            int area = horizon.area(i, j, width, lower);
            // check for new minimum
            if ((area!=-1) && ((area < minArea) || (j<y))) {
               minArea = area;
               x = i;
               y = j;
               // set stop condition for j
               //i = atlas.resolution;
               //break;
            }
         }
      }
   }

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
   void chartpack(Horizon &horizon, int numEdges, Edge *edges, int numTextures, Tex2 *textures)
   {
      // local copy of texture coordinates
      int numTexCoord = numEdges * 2;
      Tex2 *srcCoord  = new Tex2[numTexCoord];
      Tex2 *dstCoord  = new Tex2[numTexCoord];

      // chart top and bottom variables
      int *upper = new int[horizon.resolution];
      int *lower = new int[horizon.resolution];
      int size   = 0;

      // copy texture coordinates
      for (int i=0; i<numEdges; i++) {
         srcCoord[i+i]   = textures[edges[i].a];
         srcCoord[i+i+1] = textures[edges[i].b];
      }
      
      // uv middle for rotation
      Tex2 uvmid = miduv(numTexCoord, srcCoord);

      // rotation parameters
      int steps   = 16;
      float ainc  = 360.0f / steps;
      float angle = 0.0f;

      // best values go here
      int best_area    = horizon.area();
      float best_angle = 0.0f;
      int best_x       = 0;
      int best_y       = 0;
      int *best_upper  = new int[horizon.resolution];
      int best_size    = 0;

      // check every position
      for (int i=0; i<steps; i++) {
         // rotate and get lower left
         rotate(dstCoord, srcCoord, numTexCoord, angle, uvmid);
         Tex2 uvmin = minuv(numTexCoord, dstCoord);

         // clear
         for (int j=0; j<horizon.resolution; j++) {
            upper[j] = 0;
            lower[j] = horizon.resolution;
         }
         int width  = 0;
         int height = 0;

         // create a discret representation
         for (int j=0; j<numTexCoord; j+=2) {
            // get transformed coordinates
            int u1 = (int)((dstCoord[j].u   - uvmin.u) * horizon.scale);
            int v1 = (int)((dstCoord[j].v   - uvmin.v) * horizon.scale);
            int u2 = (int)((dstCoord[j+1].u - uvmin.u) * horizon.scale);
            int v2 = (int)((dstCoord[j+1].v - uvmin.v) * horizon.scale);

            width = MAX(width, u1 + 1);
            width = MAX(width, u2 + 1);

            discretize(upper, lower, u1, v1, u2, v2);
         }

         // get min max height of chart
         int ymin = horizon.resolution;
         int ymax = 0;
         for (int j=0; j<width; j++) {
            ymin = MIN(ymin, lower[j]);
            ymax = MAX(ymax, upper[j]);
         }
         height = ymax - ymin;

         // find best position
         int x;
         int y;
         int area;
         minarea(x, y, area, horizon, width, height, lower);

         // check for new candidate
         if (area!=-1 && area < best_area) {
            // set as new candidate
            best_area  = area;
            best_angle = angle;
            best_x     = x;
            best_y     = y;
            best_size  = width;
            memcpy(best_upper, upper, sizeof(int) * best_size);
         }

         // increment angle
         angle += ainc;
      }

      // use best one to update horizon
      horizon.update(best_x, best_y, best_size, best_upper);

      // rotate
      rotate(dstCoord, srcCoord, numTexCoord, best_angle, uvmid);
      rotate(textures, textures, numTextures, best_angle, uvmid);
      // transformation variables
      Tex2 uvmin = minuv(numTexCoord, dstCoord);

      printf("best position at %d, %d angle[%lf]\n", best_x, best_y, best_angle);

      // correct texture coordinates
      for (int i=0; i<numTextures; i++) {
         textures[i].u -= uvmin.u;
         textures[i].v -= uvmin.v;

         textures[i].u *= horizon.scale;
         textures[i].v *= horizon.scale;

         textures[i].u /= (float) horizon.resolution;
         textures[i].v /= (float) horizon.resolution;

         textures[i].u += (float) (best_x + 0.5f) / (float) (horizon.resolution - 1);
         textures[i].v += (float) (best_y + 0.5f) / (float) (horizon.resolution - 1);
      }
      // cleanup
      delete[] srcCoord;
      delete[] dstCoord;
      delete[] upper;
      delete[] lower;
      delete[] best_upper;
   }

} // end namespace