//=============================================================================
//
// CAtlasPacker : Class Definition
//
// Packen von Texturstücken zu einem recjteckigen
// Texturatlas
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 7 $
// 
//=============================================================================

// C++ includes
#include <list>

// local includes
#include "catlaspacker.h"
#include "../utility.h"

namespace chartification {

// ----------------------------------------------------------------------------
// comparator for cchart compare
// ----------------------------------------------------------------------------
class Cmp_Charts
{
public:
   int operator()(CChart * &c1, CChart * &c2)
   {
      return c2->FaceCount() < c1->FaceCount();
   }
};

// ----------------------------------------------------------------------------
// constructs an Atlaspacker with given
// discreet resolution, points for vertex resolution
// ----------------------------------------------------------------------------
CAtlasPacker::CAtlasPacker(Vec3 *points, int discreetSize)
{
   this->points       = points;
   this->discreetSize = discreetSize;
   this->horizon      = new CHorizon(discreetSize);

   // create upper and lower chart horizon
   upper = new int[discreetSize * 2];
   lower = new int[discreetSize * 2];
}

// ----------------------------------------------------------------------------
// deletes an CAtlasPacker
// ----------------------------------------------------------------------------
CAtlasPacker::~CAtlasPacker()
{
   if (horizon) delete horizon;
   if (upper)   delete[] upper;
   if (lower)   delete[] lower;
}

// ----------------------------------------------------------------------------
// returns the oprimal chart position in discrete atlas coordinates
// ----------------------------------------------------------------------------
bool CAtlasPacker::FindOptimal(int &x, int &y, int &minArea, int width, int height)
{
   // create return value
   x  = -1;
   y  = -1;

   // set wasted area to maximum
   minArea = horizon->Area();

   for (int i=0; i<discreetSize - width; i++) {
      for (int j=(*horizon)[i]; j<discreetSize - height; j++) {
         // calculate area
         int area = horizon->Area(i, j, width, lower);

         // check for new minimum
         if ((area!=-1) && ((area < minArea) || (j<y))) {
            minArea = area;
            x = i;
            y = j;

            return true;
         }
      }
   }
   return (x!=-1) && (y!=-1);
}

// ----------------------------------------------------------------------------
// maximum of u, v
// ----------------------------------------------------------------------------
Tex2 CAtlasPacker::Max(int size, Vertex *vertices)
{
   // initialize
   Tex2 max;
   max.u = vertices[0].u;
   max.v = vertices[0].v;

   // search maximum
   for (int i=1; i<size; i++) {
      max.u = MAX(max.u, vertices[i].u);
      max.v = MAX(max.v, vertices[i].v);
   }
   return max;
}

// ----------------------------------------------------------------------------
// middle ov u, v
// ----------------------------------------------------------------------------
Tex2 CAtlasPacker::Mid(int size, Vertex *vertices)
{
   Tex2 mid;

   Tex2 min = Min(size, vertices);
   Tex2 max = Max(size, vertices);

   mid.u = (max.u - min.u) / 2.0f;
   mid.v = (max.v - min.v) / 2.0f;

   return mid;
}

// ----------------------------------------------------------------------------
// minimum of u, v
// ----------------------------------------------------------------------------
Tex2 CAtlasPacker::Min(int size, Vertex *vertices)
{
   // initialize
   Tex2 min;
   min.u = vertices[0].u;
   min.v = vertices[0].v;

   // search minimum
   for (int i=1; i<size; i++) {
      min.u = MIN(min.u, vertices[i].u);
      min.v = MIN(min.v, vertices[i].v);
   }
   return min;
}


// ----------------------------------------------------------------------------
// pack this chartset, return packing ratio
// ----------------------------------------------------------------------------
int CAtlasPacker::Pack(CChartSet *charts)
{
   // get chart Area
   float area = 1.0f;//TextureArea(charts);
   printf("Packing Area %lf\n", area);
   //
   // calculate scale factor
   float ratio = 4.0f;
   float scale = ratio / sqrtf(area);

   ScaleProportional(charts);

   // create a sorted list
   std::list<CChart *> list;
   //
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // insert chart into list
      list.push_back(*it);
   }
   list.sort(Cmp_Charts());

   bool success = true;
   do {
      // initialize flags
      success  = true;
      ratio   -= 0.01f;
      scale    = ratio / sqrtf(area);
      horizon->Reset();

      // pack all charts
      for (std::list<CChart *>::iterator it=list.begin(); it!=list.end(); it++) {
         show_progress(); printf(" %4.2lf%%     ", ratio * 100.0);
         // pack each chart
         success &= Pack((*it), ratio, false);
         // check for success
         if (!success) break;
      }
   } while (!success);
   printf("\r");

   horizon->Reset();

   // pack all charts
   for (std::list<CChart *>::iterator jt=list.begin(); jt!=list.end(); jt++) {
      // pack each chart
      Pack((*jt), ratio, true);
   }
   printf("packing ratio %4.2lf%%     \n\n", ratio * 100.0);
   return (int) scale;
}


// ----------------------------------------------------------------------------
// pack this chart
// ----------------------------------------------------------------------------
bool CAtlasPacker::Pack(CChart *chart, float scale, bool write)
{
   if (IsFZero(TextureArea(chart)))
      return true;

   // get parameters from chart
   Vertex *vertices = chart->Vertices();
   int  vertexCount = chart->VertexCount();
   //
   Edge      *edges = chart->Boundaries();
   int    edgeCount = chart->BoundaryCount();

   // get 2 vertices per edge to represent boundary structure
   int boundaryCount = edgeCount * 2;
   Vertex *src       = new Vertex[boundaryCount];
   Vertex *dst       = new Vertex[boundaryCount];
   //
   for (int i=0; i<edgeCount; i++) {
      // copy vertices
      src[i+i]   = vertices[edges[i].a];
      src[i+i+1] = vertices[edges[i].b];
   }
      
   // uv middle for rotation
   Tex2 mid = Mid(boundaryCount, src);

   // rotation parameters
   int steps   = 16;
   float ainc  = 360.0f / steps;
   float angle = 0.0f;

   // best values go here
   int best_area    = horizon->Area();
   float best_angle = 0.0f;
   int best_x       = -1;
   int best_y       = -1;
   int *best_upper  = new int[discreetSize * 2];
   int best_size    = 0;

   // check every position
   for (int i=0; i<steps; i++) {
      // rotate and get lower left
      Rotate(dst, src, boundaryCount, angle, mid);
      //
      Tex2 min = Min(boundaryCount, dst);

      // clear upper and lower chart horizon
      for (int j=0; j<discreetSize; j++) {
         upper[j] = 0;
         lower[j] = discreetSize - 1;
      }
      int width  = 0;
      int height = 0;

      // create a discret representation
      for (int j=0; j<boundaryCount; j+=2) {
         // get transformed coordinates
         int u1 = (int)((dst[j].u   - min.u) * scale * (discreetSize - 1));
         int v1 = (int)((dst[j].v   - min.v) * scale * (discreetSize - 1));
         int u2 = (int)((dst[j+1].u - min.u) * scale * (discreetSize - 1));
         int v2 = (int)((dst[j+1].v - min.v) * scale * (discreetSize - 1));

         width = MAX(width, u1 + 1);
         width = MAX(width, u2 + 1);

         RasterLine(u1, v1, u2, v2);
      }

      if (width > discreetSize)
         continue;

      // get min max height of chart
      int ymin = discreetSize;
      int ymax = 0;
      for (int j=0; j<width; j++) {
         ymin = MIN(ymin, lower[j]);
         ymax = MAX(ymax, upper[j]);
      }
      height = ymax - ymin;

      if (height >= discreetSize)
         continue;

      // find best position
      int x, y, area;
      if (FindOptimal(x, y, area, width, height)) {
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
      }
      // increment angle
      angle += ainc;
   }

   if (best_x!=-1 && best_y!=-1) {
      // use best one to update horizon
      horizon->Update(best_x, best_y, best_size, best_upper);
   }

   // check for success
   if (write && best_x!=-1 && best_y!=-1) {
      // rotate
      Rotate(dst, src, boundaryCount, best_angle, mid);
      Rotate(vertices, vertices, vertexCount, best_angle, mid);
      // transformation variables
      Tex2 min = Min(boundaryCount, dst);

      printf("best position at %d, %d angle[%lf]\n", best_x, best_y, best_angle);

      // correct texture coordinates
      for (int i=0; i<vertexCount; i++) {
         vertices[i].u -= min.u;
         vertices[i].v -= min.v;

         vertices[i].u *= scale;
         vertices[i].v *= scale;

         vertices[i].u += (float) (best_x + 0.5f) / (float) (discreetSize - 1);
         vertices[i].v += (float) (best_y + 0.5f) / (float) (discreetSize - 1);
      }
   }
   delete[] src;
   delete[] dst;
   delete[] best_upper;

   return (best_x!=-1 && best_y!=-1);
}

// ----------------------------------------------------------------------------
// raster a line into upper and lower - min max values
// ----------------------------------------------------------------------------
void CAtlasPacker::RasterLine(int u1, int v1, int u2, int v2)
{
   // render left to right
   if (u1 > u2) {
      swap (u1, u2);
      swap (v1, v2);
   }
   
   if ((u1 < 0 || u1 >= discreetSize) || (u2 < 0 || u2 >= discreetSize)) {
      return;
   }
   if ((v1 < 0 || v1 >= discreetSize) || (v2 < 0 || v2 >= discreetSize)) {
      printf("");
   }
   int du = u2 - u1 + 1;
   int dv = v2 - v1 + 1;

   if (du<=0 || dv<=0)
      printf("");

   for (int u=u1; u<=u2; u++) {
      if (u<0 || u>=discreetSize)
         continue;
      // interpolate
      int v = v1 + (dv * (u - u1)) / du;
      if (v < 0) 
         printf("");
      // apply
      if (upper[u] < v) upper[u] = v;
      if (lower[u] > v) lower[u] = v;
   }
}

// ----------------------------------------------------------------------------
// rotate src by given angle around midpoint and copy
// result into dst
// ----------------------------------------------------------------------------
void CAtlasPacker::Rotate(Vertex *dst, Vertex *src, int numSrc, float angle, Tex2 &mid)
{
   // convert into rad
   angle *= (float)M_PI / 180.0f;
   float sina = sin(angle);
   float cosa = cos(angle);

   for (int i=0; i<numSrc; i++) {
      // load u, v pair
      float u = src[i].u - mid.u;
      float v = src[i].v - mid.v;

      // store rotated version to dst
      dst[i].u = mid.u + u * cosa - v * sina;
      dst[i].v = mid.v + u * sina + v * cosa;
   }
}


// ----------------------------------------------------------------------------
// proportionales skalieren der Charts
// ----------------------------------------------------------------------------
void CAtlasPacker::ScaleProportional(CChartSet *charts)
{
   float c_allArea = charts->Area(points);
   float t_allArea = TextureArea(charts);

   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // Jeden Chart einzeln behandeln
      CChart *chart = *it;

      // bestimmen der Fläche pro chart
      float chartArea   = chart->Area(points);
      float textureArea = TextureArea(chart);

      if (IsFZero(textureArea) || IsFZero(chartArea))
         continue;

      float cscale = sqrtf(chartArea / (c_allArea * textureArea));
      
      Vertex *vertices    = chart->Vertices();
      int     numVertices = chart->VertexCount();

      for (int i=0; i<numVertices; i++) {
         vertices[i].u *= cscale;
         vertices[i].v *= cscale;
      }
      printf("Einzel-Texturflaeche: %lf\n", TextureArea(chart));
   }
   printf("Gesamt-Texturflaeche: %lf\n", TextureArea(charts));
   printf("Gesamt-Modellflaeche: %lf\n", charts->Area(points));
}


// ----------------------------------------------------------------------------
// calculates the texture area
// ----------------------------------------------------------------------------
float CAtlasPacker::TextureArea(CChartSet *charts)
{
   // Initialize area
   float area = 0.0f;

   for (CChartSet::iterator i=charts->Begin(); i!=charts->End(); i++) {
      // sum texture area
      area += TextureArea(*i);
   }
   return area;
}

// ----------------------------------------------------------------------------
// calculates the texture area
// ----------------------------------------------------------------------------
float CAtlasPacker::TextureArea(CChart *chart)
{
   // load parameters from chart
   Vertex *vertices = chart->Vertices();
   Face   *faces    = chart->Faces();
   int     numFaces = chart->FaceCount();

   // Initialize area
   float area = 0.0f;
   
   for (int i=0; i<numFaces; i++) {
      // current face
      Vertex &a = vertices[faces[i].a];
      Vertex &b = vertices[faces[i].b];
      Vertex &c = vertices[faces[i].c];

      // sum texture area
      area += fabs(a.u*b.v - a.v*b.u) + (b.u*c.v - b.v*c.u) + (c.u*a.v - c.v*a.u);
   }
   return area * 0.5f;
}

} // end namespace
