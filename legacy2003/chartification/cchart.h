//=============================================================================
//
// CChart : Class Definition
//
// Ein Texturstück
//
// $Author: Andor $
// $Date: 7.12.03 13:34 $
// $Revision: 11 $
// 
//=============================================================================
#ifndef CCHART_H
#define CCHART_H


// Local includes
#include "..\math\math3D.h"
#include "..\math\halfedge.h"
#include "..\math\carray.h"
#include "..\math\cset.h"

namespace chartification {

class CChart
{
   // seed of chart
   int       seed;
   CSet<int> seedHistory;
   int       clusterCount;

   // variables for clustering
   Vec3     chartNormal;
   Vec3     averageNormal;

   // Containers
   CArray<Edge>   boundaries;
   CArray<Vertex> vertices;
   CArray<Face>   faces;

public:
   // accessible index
   int searchIndex;

   // creates a Chart from an initial face, its normal
   CChart(int faceIndex, Vec3 &normal);

   // destruktor
   ~CChart() {};

   // Adds a edge to the chart boundary
   void AddBoundary(int a, int b);
   // Adds a face with vertex indices
   void AddFace(int a, int b, int c);
   // adds a vertex and returns its id
   int AddVertex(int pointId);

   // calculates chart area
   float Area(Vec3 *points);

   // boundary access
   int   BoundaryCount();
   Edge *Boundaries();

   // reset everything but the seed
   void Clear();

   // number of grows
   int ClusterCount();

   // face access
   int   FaceCount();
   Face *Faces();

   // Grows the chart, and its boundary
   void Grow(Vec3 &normal);

   // returns the normal of this chart
   Vec3 Normal();

   // returns the seed index of this chart
   int Seed();

   // sets a new seed and indicates a cycle
   bool SetSeed(int seed, Vec3 &normal);

   // Vertex access
   int     VertexCount();
   Vertex *Vertices();
};

} // end namespace
#endif