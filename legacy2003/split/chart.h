/******************************************************************************
 *
 * Chart.h
 *
 * internal representation of a chart
 *
 *****************************************************************************/
#ifndef CHART_H
#define CHART_H

#include "../math/math3D.h"
#include "../math/sparsematrix.h"
#include "lscm.h"
#include <set>
#include <map>

class Chart;
typedef std::set<Chart *> chart_set;

class Chart
{
public:
   // the seed face
   int chart_seed;
   
   // set of all faces
   int_set    faces_set;

   // maximum distance
   float max_distance;
   
   // the chart normal
   Vec3 chart_normal;
   Vec3 average_nrm;
   
   // chart boundaries
   halfedge_set chart_boundaries;
   
   // flag indicates wheter chart is linked or not
   bool linked;

   // area and perimeter information
   float chart_area;
   float chart_perimeter;

public:
   // define iterator (use set itertor)
   typedef int_set::iterator face_iterator;
   typedef halfedge_set::iterator halfedge_iterator;

   /***************************************************************************
    * Construct a chart from a face index and distance_to_feature value
    **************************************************************************/
   Chart (int face, float distance)
   {
      faces_set.insert(face);

      chart_seed      = face;
      max_distance    = distance;
      chart_area      = 0.0f;
      chart_perimeter = 0.0f;
   }

   /***************************************************************************
    * Construct a chart from a seed index and normal value
    **************************************************************************/
   Chart (int seed, Vec3 nrm)
   {
      faces_set.insert(seed);

      // set seed and initial face normal
      chart_seed      = seed;
      chart_normal    = nrm;
      average_nrm     = nrm;
      max_distance    = 0.0f;
      chart_area      = 0.5f * nrm.length();
      chart_perimeter = 0.0f;
   }

   /***************************************************************************
    * adds a face to the chart
    **************************************************************************/
   void add(int face, float distance)
   {
      faces_set.insert(face);
      max_distance = MAX(max_distance, distance);
   }

   /***************************************************************************
    * adds a face and a normal the chart
    **************************************************************************/
   void add(int face, Vec3 nrm)
   {
      faces_set.insert(face);
      average_nrm += nrm;
      // keep track of the perimeter
      chart_area += 0.5f * nrm.length();
   }

   /***************************************************************************
    * returns the area cobert by this chart
    **************************************************************************/
   float area()
   {
      return chart_area;
   }

   /***************************************************************************
    * delivers length of the halfedge
    **************************************************************************/
   float length(Halfedge *w, Vec3 *points)
   {
      return (points[w->end] - points[w->next->next->end]).length();
   }

   /***************************************************************************
    * adds a halfedge to the chart
    **************************************************************************/
   void add(Halfedge *w, Vec3 *points)
   {
      chart_boundaries.insert(w);
      chart_perimeter += length(w, points);
   }

   /***************************************************************************
    * removes a halfedge from chart_boundaries
    **************************************************************************/
   void erase(Halfedge *w, Vec3 *points)
   {
      chart_boundaries.erase(w);
      chart_perimeter -= length(w, points);
   }

   /***************************************************************************
    * checks if the chart contains the specific face
    **************************************************************************/
   bool contains(int face_index)
   {
      return faces_set.find(face_index) != faces_set.end();
   }

   /***************************************************************************
    * merges this chart with another chart
    ***************************************************************************/
   void merge (Chart *chart)
   {
      // merge face sets
      for (int_set::iterator it=chart->faces_set.begin(); it!=chart->faces_set.end(); ++it)
         faces_set.insert(*it);
      // correct the maximum
      max_distance = MAX(max_distance, chart->max_distance);
   }

   /***************************************************************************
    * return maximum distance to feature
    **************************************************************************/
   float max_dist()
   {
      return max_distance;
   }

   /***************************************************************************
    * begin of face indices
    **************************************************************************/
   face_iterator begin_face()
   {
      return faces_set.begin();
   }

   /***************************************************************************
    * end of face indices
    **************************************************************************/
   face_iterator end_face()
   {
      return faces_set.end();
   }

   /***************************************************************************
    * begin of halfedges
    **************************************************************************/
   halfedge_iterator begin_halfedge()
   {
      return chart_boundaries.begin();
   }

   /***************************************************************************
    * end of halfedges
    **************************************************************************/
   halfedge_iterator end_halfedge()
   {
      return chart_boundaries.end();
   }

   /***************************************************************************
    * returns the average normal of all faces
    **************************************************************************/
   Vec3 normal()
   {
      Vec3 nrm = chart_normal;
      nrm.norm();
      return nrm;
   }

   /***************************************************************************
    * returns the average normal of all faces
    **************************************************************************/
   Vec3 average_normal()
   {
      Vec3 nrm = average_nrm;
      return nrm.norm();
   }

   /***************************************************************************
    * returns the chart seed index
    **************************************************************************/
   int seed()
   {
      return this->chart_seed;
   }

   /***************************************************************************
    * apply a mapping from model to texture space
    *
    * algorithm Hoppe et al.
    *
    **************************************************************************/
   void parametrize_lscm(
      Vec3 *srcpoints, Face *srcfaces, int pinned1, int pinned2,
      int *numPoints, int *numFaces, Vec3 **points, Face **faces, Tex2 **texture, int *numPath, int **path
      )
   {
      // create a LSCM parametrization object
      LSCM lscm;
      lscm.init(&faces_set, &chart_boundaries, srcpoints, srcfaces, pinned1, pinned2, chart_seed);
      lscm.parametrize(numPoints, numFaces, points, faces, texture);

      lscm.getPath(numPath, path);

      float textureArea = 0.0f;
      float umin = (*texture)[0].u;
      float umax = umin;
      float vmin = (*texture)[0].v;
      float vmax = vmin;

      chart_area = 0.0f;

      for (int i=0; i<(*numFaces); i++) {
         // add face area to texture and chart area
         Tex2 &a = (*texture)[(*faces)[i].a];
         Tex2 &b = (*texture)[(*faces)[i].b];
         Tex2 &c = (*texture)[(*faces)[i].c];

         textureArea += (a.u*b.v - a.v*b.u) + (b.u*c.v - b.v*c.u) + (c.u*a.v - c.v*a.u);
         chart_area  += normal_raw(*points, &((*faces)[i])).length();

         // calculate texture middle
         umin = MIN(umin, a.u); umax = MAX(umax, a.u);
         umin = MIN(umin, b.u); umax = MAX(umax, b.u);
         umin = MIN(umin, c.u); umax = MAX(umax, c.u);

         vmin = MIN(vmin, a.v); vmax = MAX(vmax, a.v);
         vmin = MIN(vmin, b.v); vmax = MAX(vmax, b.v);
         vmin = MIN(vmin, c.v); vmax = MAX(vmax, c.v);
      }
      // translate and scale
      float umid  = (umin + umax) / 2.0f;
      float vmid  = (vmin + vmax) / 2.0f;
      float scale = sqrtf(chart_area / textureArea);

      for (int i=0; i<(*numPoints); i++) {
         (*texture)[i].u -= umid;
         (*texture)[i].v -= vmid;
         (*texture)[i].u *= scale;
         (*texture)[i].v *= scale;
      }
   }

   /***************************************************************************
    * returns the perimeter of this chart
    **************************************************************************/
   float perimeter()
   {
      return chart_perimeter;
   }
};



#endif



