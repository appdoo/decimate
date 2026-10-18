// ============================================================================
//
// Normal - Sampling
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 5 $
// 
// ============================================================================

// Local Includes
#include "texturesample.h"
#include "../math/algorithm.h"
#include "../raytracing/raytracing.h"
#include "../utility.h"

namespace normalmap {

// ----------------------------------------------------------------------------
// structure representing a barycentric coordinate system
// ----------------------------------------------------------------------------
typedef struct bc3
{
   float u, v, w;
} bc3;

// ----------------------------------------------------------------------------
// test the barycentric criteria to determin, if a point is inside
// ----------------------------------------------------------------------------
inline bool inside(float u, float v, float w)
{
   return !((u < 0.0 || v < 0.0 || w < 0.0) || (u > 1.0 || v > 1.0 || w > 1.0));
}

// ----------------------------------------------------------------------------
// calculates the barycentric coordinates of p0
// within the givemn triangle (p1, p2, p3)
// ----------------------------------------------------------------------------
inline bc3 barycentric(Tex2 &p0, Tex2 &p1, Tex2 &p2, Tex2 &p3)
{
   bc3 bc;

   float s =  (p2.u - p1.u) * (p3.v - p1.v) - (p3.u - p1.u) * (p2.v - p1.v);
   bc.u    = ((p2.u - p0.u) * (p3.v - p0.v) - (p3.u - p0.u) * (p2.v - p0.v)) / s;
   bc.v    = ((p3.u - p0.u) * (p1.v - p0.v) - (p1.u - p0.u) * (p3.v - p0.v)) / s;
   bc.w    = ((p1.u - p0.u) * (p2.v - p0.v) - (p2.u - p0.u) * (p1.v - p0.v)) / s;

   return bc;
}

// ----------------------------------------------------------------------------
// calculates the barycentric coordinates of p0
// within the givemn triangle (p1, p2, p3)
// ----------------------------------------------------------------------------
inline bc3 barycentric(float u, float v, Tex2 &p1, Tex2 &p2, Tex2 &p3, float bdiv)
{
   bc3 bc;

   bc.u    = ((p2.u - u) * (p3.v - v) - (p3.u - u) * (p2.v - v)) / bdiv;
   bc.v    = ((p3.u - u) * (p1.v - v) - (p1.u - u) * (p3.v - v)) / bdiv;
   bc.w    = ((p1.u - u) * (p2.v - v) - (p2.u - u) * (p1.v - v)) / bdiv;

   return bc;
}

// ----------------------------------------------------------------------------
// do a range compress on the supplied normal and store it
// into the rgb - value
// ----------------------------------------------------------------------------
inline void rangecompress(rgb &nrmrgb, Vec3 &nrm)
{
   nrm.norm();

	nrmrgb.r = (unsigned char)((nrm.x + 1.0f) * 127.0f);
	nrmrgb.g = (unsigned char)((nrm.y + 1.0f) * 127.0f);
	nrmrgb.b = (unsigned char)((nrm.z + 1.0f) * 127.0f);
}

// ----------------------------------------------------------------------------
// get the normal at the model of higher order
// ----------------------------------------------------------------------------
bool rtraceHighOrderNormal(Vec3 &nrm, Vec3 orig, Vec3 &dir, int numFaces, Face *faces, Vec3 *points, Vec3 *normals)
{
   float mint = 100000.0;
   float minu = 0.0f;
   float minv = 0.0f;
   int minface = -1;
   float u, v, t;

   // for all faces
   for (int i=0; i<numFaces; i++) {
      // reference face
      Face &f = faces[i];

      // test for intersection
      if (rtrace::intersect_triangle(orig, dir, points[f.a], points[f.b], points[f.c], &t, &u, &v)) {
         if (fabs(t) < fabs(mint)) {
            // set new minimum
            mint    = t;
            minface = i;
            minu = u;
            minv = v;
         }
      }
   }
   // found nothing ?
   if (minface<0) 
   {
      nrm = dir;
      return true;
   }

   // reference face
   Face &f = faces[minface];

   u = minu;
   v = minv;

   nrm = (normals[f.a] * (1.0f - (u + v))) + (normals[f.b] * u) + (normals[f.c] * v);

   return true;
}


// ----------------------------------------------------------------------------
// Vec3 interpolation
// ----------------------------------------------------------------------------
inline void interpolate_vec3(Vec3 &dest, const bc3 &bc, const Vec3 &v1, const Vec3 &v2, const Vec3 &v3)
{
   dest.x = bc.u * v1.x + bc.v * v2.x + bc.w * v3.x;
   dest.y = bc.u * v1.y + bc.v * v2.y + bc.w * v3.y;
   dest.z = bc.u * v1.z + bc.v * v2.z + bc.w * v3.z;
}


// ----------------------------------------------------------------------------
// samples a triangle
// ----------------------------------------------------------------------------
void mapfacert(int texSize, rgb *mem, 
             Face *f, Tex2 *tex, Vec3 *lowPoints, Vec3 *lowNormals,
             Vec3 *highPoints, Vec3 *highNormals, int numHighFaces, Face *highFaces)
{
   // load texture references
   Tex2 t1 = tex[f->a];
   Tex2 t2 = tex[f->b];
   Tex2 t3 = tex[f->c];

   // scale to fit texture size
   t1.u *= texSize - 1; t1.v *= texSize - 1;
   t2.u *= texSize - 1; t2.v *= texSize - 1;
   t3.u *= texSize - 1; t3.v *= texSize - 1;

   Vec3 &n1 = lowNormals[f->a];
   Vec3 &n2 = lowNormals[f->b];
   Vec3 &n3 = lowNormals[f->c];

   // calculate fix part of the barycentric coordinates
   float bary = (t2.u - t1.u) * (t3.v - t1.v) - (t3.u - t1.u) * (t2.v - t1.v);

   // get min, max
   int minu = (int) (MIN(t1.u, MIN(t2.u, t3.u)));
   int maxu = (int) (MAX(t1.u, MAX(t2.u, t3.u)) + 0.5f);
   int minv = (int) (MIN(t1.v, MIN(t2.v, t3.v)));
   int maxv = (int) (MAX(t1.v, MAX(t2.v, t3.v)) + 0.5f);

   minu = MAX(minu, 0);  minu = MIN(minu, texSize - 1);
   maxu = MAX(maxu, 0);  maxu = MIN(maxu, texSize - 1);
   minv = MAX(minv, 0);  minv = MIN(minv, texSize - 1);
   maxv = MAX(maxv, 0);  maxv = MIN(maxv, texSize - 1);

   // iterate over the area
   for (int x=minu; x<=maxu; x++) {
      for (int y=minv; y<=maxv; y++) {

         Vec3 n, p, snrm;

         bc3 bc = barycentric((float)x , (float)y, t1, t2, t3, bary);
         // check if we have a valid point ( 0 <= u, v, w <= 1)
         if (!inside(bc.u, bc.v, bc.w))
            continue;

         // interpolate normal, origin, tangent space
         interpolate_vec3(n, bc, n1, n2, n3);
         interpolate_vec3(p, bc, lowPoints[f->a], lowPoints[f->b], lowPoints[f->c]);
               
         if (rtraceHighOrderNormal(snrm, p, n, numHighFaces, highFaces, highPoints, highNormals)) {
            // range compressed value
            snrm.norm();
            rangecompress(mem[y * texSize + x], snrm);
         }
      }
   }
}


// ----------------------------------------------------------------------------
// sets the texture memory to the specified value
// ----------------------------------------------------------------------------
void cleartexture(int texSize, rgb *mem, unsigned char r, unsigned char g, unsigned char b)
{
   const int size = texSize * texSize;

   for (int i=0; i<size; i++) {
      mem[i].r = r;
      mem[i].g = g;
      mem[i].b = b;
   }
}


// ----------------------------------------------------------------------------
// Lässt die Grenzen eines Bildes zerfliessen
// ----------------------------------------------------------------------------
void dilate (rgb* img, int texsize)
{
   rgb *scratch = new rgb[texsize * texsize];
   rgb *alloc   = scratch;

   rgb rgb_zero;
   rgb_zero.init(128, 128, 128);

   for (int f = 0; f < 6; f++)
   {
      // Loop over the old image and create a new one.
      for (int y = 0; y < texsize; y++) {
         for (int x = 0; x < texsize; x++) {
            // check for empty pixel
            int idx = y * texsize + x;
            if (img[idx] == rgb_zero) {
               // 
               int accu[3] = {0, 0, 0};
               int acnt    = 0;
               
               // loop over near pixels
               for (int i=x-1; i<x+2; i++) {
                  for (int j=y-1; j<y+2; j++) {
                     // check for valid pixel
                     if (i==x && j==y)
                        continue;
                     // get near index
                     int xn = i % texsize;
                     int yn = j % texsize;
                     if (xn < 0) xn = texsize - 1;
                     if (yn < 0) yn = texsize - 1;
                     const int nx = yn * texsize + xn;
                     // check for filled pixel
                     if (img[nx] != rgb_zero) {
                        // accumulate
                        accu[0] += img[nx].r;
                        accu[1] += img[nx].g;
                        accu[2] += img[nx].b;
                        //
                        acnt++;
                     }
                  }
               }
               // fill pixel with average
               if (acnt > 0) {
                  scratch[idx].r = accu[0] / acnt;
                  scratch[idx].g = accu[1] / acnt;
                  scratch[idx].b = accu[2] / acnt;
               }
               else {
                  scratch[idx] = img[idx];
               }
            } // end if this pixel is empty
            else {
                  scratch[idx] = img[idx];
            }
         } // end for x
      } // end for y
      rgb *tmp = img;
      img = scratch;
      scratch = tmp;
   } // end for f
   delete[] alloc;
}

bool contains(Vec3 &p, Vec3 *nrms, float *d)
{
   bool ba = ((nrms[0] * p)- d[0]>=0) || ((nrms[1] * p)- d[1]>=0);
   bool cb = ((nrms[2] * p)- d[2]>=0) || ((nrms[3] * p)- d[3]>=0);
   bool ac = ((nrms[4] * p)- d[4]>=0) || ((nrms[5] * p)- d[5]>=0);

   return ba && cb && ac;
}

// ----------------------------------------------------------------------------
// Einschränkung der zu testeten Flächen
// ----------------------------------------------------------------------------
void getTestFaces(Face *testFaces, int *numTestFaces, Face &f, Vec3 *lowPoints, Vec3 *lowNormals, Vec3 *highPoints, int numHighFaces, Face *highFaces)
{
   *numTestFaces = 0;

   Vec3 nrms[6];
   float d[6];

   Vec3 a = lowPoints[f.a];
   Vec3 b = lowPoints[f.b];
   Vec3 c = lowPoints[f.c];

   Vec3 m = (a + b + c) / 3.0f;
   float scale = 1.5f;

   a = ((a - m) * scale) + m;
   b = ((b - m) * scale) + m;
   c = ((c - m) * scale) + m;

   Vec3 ab = a - b;
   Vec3 bc = b - c;
   Vec3 ca = c - a;

   nrms[0] = (ab ^ lowNormals[f.a]).norm();
   nrms[1] = (ab ^ lowNormals[f.b]).norm();
   nrms[2] = (bc ^ lowNormals[f.b]).norm();
   nrms[3] = (bc ^ lowNormals[f.c]).norm();
   nrms[4] = (ca ^ lowNormals[f.c]).norm();
   nrms[5] = (ca ^ lowNormals[f.a]).norm();

   d[0] = (nrms[0] * a);
   d[1] = (nrms[1] * b);
   d[2] = (nrms[2] * b);
   d[3] = (nrms[3] * c);
   d[4] = (nrms[4] * c);
   d[5] = (nrms[5] * a);

   for (int i=0; i<numHighFaces; i++) {
      Face &o = highFaces[i];
      if (contains(highPoints[o.a], nrms, d) || contains(highPoints[o.b], nrms, d) || contains(highPoints[o.c], nrms, d))
      {
         testFaces[(*numTestFaces)++] = highFaces[i];
      }
   }
}


// ----------------------------------------------------------------------------
// computes normal values for texture map
// ----------------------------------------------------------------------------
void mapnormals(
                   int texSize,
                   rgb *mem,
                   int  numLowPoints,
                   Vec3 *lowPoints,
                   Vec3 *lowNormals,
                   Tex2 *lowTex,
                   int   numLowFaces, 
                   Face *lowFaces,

                   Vec3 *highPoints,
                   Vec3 *highNormals,
                   int  numHighFaces,
                   Face *highFaces
                   )
{
   // create an up - normal vector
   cleartexture(texSize, mem, 128, 128, 128);

   Face *testFaces = new Face[numHighFaces];
   int numTestFaces = 0;

   setup_timer();
   start_timer();

   for (int i=0; i<numLowFaces; i++) {
      getTestFaces(testFaces, &numTestFaces, lowFaces[i], lowPoints, lowNormals, highPoints, numHighFaces, highFaces);
      try {
      mapfacert(texSize, mem, &lowFaces[i], lowTex, lowPoints, lowNormals, highPoints, highNormals, numTestFaces, testFaces);
      }
      catch (...)
      {
         printf("mapnormals : exception in face i\n");
      }
      printf("\rRay casting face %d/%d", i, numLowFaces);
   }
   printf("\n");
   double time = stop_timer();
   printf("Raytracing %lf milliseconds\n", time);

   delete[] testFaces;

   dilate(mem, texSize);
}

} // end namespace