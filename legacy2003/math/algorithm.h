/******************************************************************************
 *
 * algorithm.h
 *
 * algorithms and utility functions
 *
 *****************************************************************************/
#ifndef ALGORITHM_H
#define ALGORITHM_H

// Local Includes
#include "math3D.h"
#include "vec.h"
#include "halfedge.h"


// converts two given vectors into an orthonormal basis
// For fetails of the algorithm see definition
// of the Gram-Schmidt process to define an orthonormal basis
inline void orthonormal(Vec3 &w1, Vec3 &w2)
{
   w1.norm();
   w2 = w2 - (w1 * (w2 * w1));
   w2.norm();
}

// converts two given vectors into an orthonormal basis
// For fetails of the algorithm see definition
// of the Gram-Schmidt process to define an orthonormal basis
inline Vec3 orthonormal3(Vec3 &w1, Vec3 &w2)
{
   w1.norm();
   w2 = w2 - (w1 * (w2 * w1));
   w2.norm();
   //
   return w1 ^ w2;
}

// calculate face normal in unit length
inline Vec3 normal(Face *f, Vec3 *points)
{
   Vec3 n;

   Vec3 b = points[f->b] - points[f->a];
   Vec3 c = points[f->c] - points[f->a];

   n.x = b.y * c.z - b.z * c.y;
   n.y = b.z * c.x - b.x * c.z;
   n.z = b.x * c.y - b.y * c.x;

   return n.norm();
}

// calculate face normal
inline Vec3 normal_raw(Vec3 *points, Face *f)
{
   Vec3 n;

   Vec3 b = points[f->b] - points[f->a];
   Vec3 c = points[f->c] - points[f->a];

   n.x = b.y * c.z - b.z * c.y;
   n.y = b.z * c.x - b.x * c.z;
   n.z = b.x * c.y - b.y * c.x;

   return n;
}

// calculate face - middlepoint
inline Vec3 middle(Vec3 *points, Face *f)
{
   return (((points[f->a] + points[f->b]) * 0.5f) + points[f->c]) * 0.5f;
}

// calculate length of a halfedge
inline float length(Halfedge *h, Vec3 *points)
{
   return (points[h->end] - points[h->next->next->end]).length();
}


#endif