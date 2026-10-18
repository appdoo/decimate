//=============================================================================
//
// CAtlasPacker : Class Definition
//
// Packen von Texturstücken zu einem recjteckigen
// Texturatlas
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 4 $
// 
//=============================================================================
#ifndef CATLASPACKER_H
#define CATLASPACKER_H


// local includes
#include "cchart.h"
#include "cchartset.h"
#include "chorizon.h"

namespace chartification {

class CAtlasPacker
{
   // points reference
   Vec3 *points;
   // discreet size
   int discreetSize;
   // Horizon class
   CHorizon *horizon;

   // upper and lower
   // chart horizon
   int *upper;
   int *lower; 

public:
   // constructs an Atlaspacker with given
   // discreet resolution, points for vertex resolution
   CAtlasPacker(Vec3 *points, int discreetSize = 256);

   // deletes an CAtlasPacker
   ~CAtlasPacker();

   // pack this chartset, return packing ratio
   int Pack(CChartSet *charts);

private:
   // returns the oprimal chart position in discrete atlas coordinates
   bool FindOptimal(int &x, int &y, int &minArea, int width, int height);

   // maximum of u, v
   Tex2 Max(int size, Vertex *vertices);

   // minimum of u, v
   Tex2 Min(int size, Vertex *vertices);

   // middle ov u, v
   Tex2 Mid(int size, Vertex *vertices);

   // Pack a chart
   bool Pack(CChart *chart, float scale, bool write);

   // raster a line into upper and lower - min max values
   void RasterLine(int u1, int v1, int u2, int v2);

   // rotate src by given angle into dst
   void Rotate(Vertex *dst, Vertex *src, int numSrc, float angle, Tex2 &mid);

   // proportionales skalieren der Charts
   void ScaleProportional(CChartSet *charts);

   // calculates the texture area
   float TextureArea(CChartSet *charts);

   // calculates the texture area
   float TextureArea(CChart *chart);

};

} // end namespace
#endif
