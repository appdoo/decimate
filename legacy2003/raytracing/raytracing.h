/******************************************************************************
 *
 * Implementation of ray tracing, or ray casting algorithms
 *
 *
 * $Author: Andor $
 * $Date: 24.11.03 11:10 $
 * $Revision: 2 $
 * 
 *****************************************************************************/
#ifndef RAYTRACING_H
#define RAYTRACING_H

// Local Includes
#include "../math/math3D.h"
#include "../math/vec.h"

namespace rtrace {


   /***************************************************************************
    *
    * Routine zum Errechnen des Schnittpunkt - Parameters t
    *
    *        d + N * O
    * t = - -----------
    *          N * D
    *
    * Mit N: Normalenvektor der Ebene
    *     d: Parameter der Ebenengleichung N * P + d = 0
    *     O: Ursprung des Strahls
    *     D: Richtungsvektor des Strahls
    *
    **************************************************************************/
   bool evaluate_t(float &t, const Vec3 &nrm, const float d, const Vec3 &origin, const Vec3 &dir);

   /***************************************************************************
    *
    * An Efficient Ray/Polygon Intersection
    * by Didier Badouel
    * from "Graphics Gems", Academic Press, 1990
    *
    **************************************************************************/
   bool barycentricBadouel(float &alpha, float &beta, const float *p, const float *nrm, const float *p0, const float *p1, const float *p2);


   /***************************************************************************
    * Basis Routine zum Berechnen des Schnittpunkts mit einem Strahl
    *
    * Implementation of method by 
    *
    * Fast, Minimum Storage Ray/Triangle Intersection
    *          Tomas Möller, Ben Trumbore
    *
    * Paper can be found at
    * http://www.graphics.cornell.edu/pubs/1997/MT97.html
    * 
    * Version with backface culling enabled
    *
    ***************************************************************************/
   bool intersect_triangle_bfc(
      Vec3 &orig, Vec3 &dir,
      Vec3 &v0, Vec3 &v1, Vec3 &v2,
      float *t, float *u, float *v);

   /***************************************************************************
    * Basis Routine zum Berechnen des Schnittpunkts mit einem Strahl
    *
    * Implementation of method by 
    *
    * Fast, Minimum Storage Ray/Triangle Intersection
    *          Tomas Möller, Ben Trumbore
    *
    * Paper can be found at
    * http://www.graphics.cornell.edu/pubs/1997/MT97.html
    * 
    * Version with backface culling disabled
    *
    ***************************************************************************/
   bool intersect_triangle(
      Vec3 &orig, Vec3 &dir,
      Vec3 &v0, Vec3 &v1, Vec3 &v2,
      float *t, float *u, float *v);

      //////////////////////////////////////////////////////////////////////////
   // See if the ray intersects the triangle.
   // Drew's adaption from Real Time Rendering
   // Optimized based on Moller's work in http://www.acm.org/jgt/
   // http://www.ce.chalmers.se/staff/tomasm/raytri/raytri.c
   //////////////////////////////////////////////////////////////////////////
   bool
   IntersectTriangle (Vec3 &orig, Vec3 &dir,
                     Vec3 &v1, Vec3 &v2, Vec3 &v3,
                     float *t, float *u, float *v);
} // end namespace

#endif


