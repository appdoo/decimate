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

// Local Includes
#include "raytracing.h"

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
   bool evaluate_t(float &t, const Vec3 &nrm, const float d, const Vec3 &origin, const Vec3 &dir)
   {
      float div = nrm * dir;

      // pruefen ob strahl und ebene parallel
      if (fabs(div) < EPSILONF) 
         return false;
      t = - (d + (nrm * origin)) / div;

      return true;
   }

   /***************************************************************************
    *
    * An Efficient Ray/Polygon Intersection
    * by Didier Badouel
    * from "Graphics Gems", Academic Press, 1990
    *
    * http://www.acm.org/pubs/tog/GraphicsGems/gems.html#gems
    *
    **************************************************************************/
   bool barycentricBadouelf(float &alpha, float &beta, const float *p, const float *nrm, const float *p0, const float *p1, const float *p2)
   {
      int i1, i2;

      // Der folgende Teil, Berechnung von i1, i2
      // kann auch vorberechnet werden
      float x = fabs(nrm[0]);
      float y = fabs(nrm[1]);
      float z = fabs(nrm[2]);

      // find the dominant axis of faces normal vector
      if      (x > y && x > z) { i1=1; i2=2; }
      else if (y > x && y > z) { i1=0; i2=2; }
      else                     { i1=0; i2=1; }

      float u0 = p[i1]  - p0[i1]; float v0 = p[i2]  - p0[i2];
      float u1 = p1[i1] - p0[i1]; float v1 = p1[i2] - p0[i2];
      float u2 = p2[i1] - p0[i1]; float v2 = p2[i2] - p0[i2];

      if (u1 == 0.0f) {
         beta = u0/u2;
         if ((beta >= 0.0f)&&(beta <= 1.0f)) {
            alpha = (v0 - beta*v2)/v1;
            return (alpha >= 0.0f) && ((alpha+beta) <= 1.0f);
         }
      }
      else {
         beta = (v0*u1 - u0*v1)/(v2*u1 - u2*v1);
         if ((beta >= 0.0f)&&(beta <= 1.0f)) {
            alpha = (u0 - beta*u2)/u1;
            return (alpha >= 0.0f) && ((alpha+beta) <= 1.0f);
        }
      }
      return false;
   }

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
      float *t, float *u, float *v)
   {
      // find vectors for two edges sharing v0
      Vec3 e1 = v1 - v0;
      Vec3 e2 = v2 - v0;

      // begin calculation determinant - also used to calculate U parameter
      Vec3 pv = dir ^ e2;

      // if determinant is near zero, ray lies in plane of triangle
      float det = e1 * pv;

      // version with backface culling enabled
      if (det < 1e-12f) return false;

      // calculate distance from v0 to ray origin
      Vec3 tv = orig - v0;

      // calculate u parameter and test bounds
      *u = tv * pv;
      if (*u < 0.0f || *u > det)
         return false;

      // prepare to test v parameter
      Vec3 qv = tv ^ e1;

      // calculate v parameter and test bounds
      *v = dir * qv;
      if (*v < 0.0f || *u + *v > det)
         return false;

      // calculate t, scale parameters, ray intersects triangle
      *t = e2 * qv;
      float inv_det = 1.0f / det;
      *t *= inv_det;
      *u *= inv_det;
      *v *= inv_det;

      return true;
   }

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
    * Version without backface culling
    *
    ***************************************************************************/
   bool intersect_triangle(
      Vec3 &orig, Vec3 &dir,
      Vec3 &v0, Vec3 &v1, Vec3 &v2,
      float *t, float *u, float *v)
   {
      // find vectors for two edges sharing v0
      Vec3 e1 = v1 - v0;
      Vec3 e2 = v2 - v0;

      // begin calculation determinant - also used to calculate U parameter
      Vec3 pv = dir ^ e2;

      // if determinant is near zero, ray lies in plane of triangle
      float det = e1 * pv;

      // calculate distance from v0 to ray origin
      Vec3 tv = orig - v0;
      // prepare to test v parameter
      Vec3 qv = tv ^ e1;

      if (det > EPSILONF) {
         // calculate u parameter and test bounds
         *u = (tv * pv);
         if (*u < 0.0f || *u > det)
            return false;

         // calculate v parameter and test bounds
         *v = (dir * qv);
         if (*v < 0.0f || *u + *v > det)
            return false;
      }
      else if (det < -EPSILONF) {
         // calculate u parameter and test bounds
         *u = (tv * pv);
         if (*u > 0.0f || *u < det)
            return false;

         // calculate v parameter and test bounds
         *v = (dir * qv);
         if (*v > 0.0f || *u + *v < det)
            return false;
      }
      else {
         return false;
      }
      // calculate t, scale parameters, ray intersects triangle
      *t  = (e2 * qv) / det;
      *u /= det;
      *v /= det;

      return true;
   }


   
   //////////////////////////////////////////////////////////////////////////
   // See if the ray intersects the triangle.
   // Drew's adaption from Real Time Rendering
   // Optimized based on Moller's work in http://www.acm.org/jgt/
   // http://www.ce.chalmers.se/staff/tomasm/raytri/raytri.c
   //////////////////////////////////////////////////////////////////////////
   bool
   IntersectTriangle (Vec3 &orig, Vec3 &dir,
                     Vec3 &v1, Vec3 &v2, Vec3 &v3,
                     float *t, float *u, float *v)
   {
      Vec3 edge1, edge2, tvec, pvec, qvec;
      float det,inv_det;

      // find vectors for two edges sharing vert0
      edge1 = v2 -  v1;
      edge2 = v3 -  v1;
      
      // begin calculating determinant - also used to calculate U parameter
      pvec = dir ^ edge2;

      // if determinant is near zero, ray lies in plane of triangle
      det = edge1 * pvec;

      // calculate distance from vert0 to ray origin
      tvec = orig - v1;

      qvec = tvec ^ edge1;

      if (det > EPSILONF)
      {
         *u = tvec * pvec;
         if (*u < 0.0f || *u > det)
         {
            return false;
         }
               
         // Calculate V parameter and test bounds
         *v = dir * qvec;
         if (*v < 0.0f || *u + *v > det)
         {
            return false;
         }
      }
      else if (det < -EPSILONF)
      {
         // Calculate U parameter and test bounds
         *u = tvec * pvec;
         if (*u > 0.0f || *u < det)
         {
            return false;
         }
         
         // Calculate V parameter and test bounds
         *v = dir * qvec;
         if (*v > 0.0f || *u + *v < det)
         {
            return false;
         }
      }
      else
      {
         // Ray is parallell to the plane of the triangle
         return false;
      }
      inv_det = 1.0f / det;

      *t = (edge2 * qvec) * inv_det;
      (*u) *= inv_det;
      (*v) *= inv_det;

      return true;
   }
} // end namespace

