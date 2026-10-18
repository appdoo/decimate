//=============================================================================
//
// CChart : Class Implementation
//
// Ein Texturstück
//
// $Author: Andor $
// $Date: 7.12.03 13:34 $
// $Revision: 12 $
// 
//=============================================================================

// Local includes
#include "..\math\algorithm.h"
#include "cchart.h"

namespace chartification {

// ----------------------------------------------------------------------------
// creates a Chart from an initial face, its normal and a reference
// to the graph structure of the mesh
// ----------------------------------------------------------------------------
CChart::CChart(int faceIndex, Vec3 &normal)
{
   seed          = faceIndex;
   chartNormal   = normal;
   averageNormal = normal;
   clusterCount  = 0;

   averageNormal.norm();
}


// ----------------------------------------------------------------------------
// Adds a edge to the chart boundary
// ----------------------------------------------------------------------------
void CChart::AddBoundary(int a, int b)
{
   Edge edge;

   edge.a = a;
   edge.b = b;

   boundaries.Add(edge);
}

// ----------------------------------------------------------------------------
// Adds a face with vertex indices
// ----------------------------------------------------------------------------
void CChart::AddFace(int a, int b, int c)
{
   Face face;

   face.a = a;
   face.b = b;
   face.c = c;

   faces.Add(face);
}

// ----------------------------------------------------------------------------
// adds a vertex and returns its id
// ----------------------------------------------------------------------------
int CChart::AddVertex(int pointId)
{
   Vertex vertex;
   //
   vertex.index = pointId;
   vertex.u     = 0.0f;
   vertex.v     = 0.0f;
   //
   return vertices.Add(vertex);
}

// ----------------------------------------------------------------------------
// calculates chart area
// ----------------------------------------------------------------------------
float CChart::Area(Vec3 *points)
{
   // initialize area
   float area = 0.0f;
   // point mapped face variable
   Face  pmFace;

   for (int i=0; i<faces.Size(); i++) {
      // map face
      pmFace.a = vertices[faces[i].a].index;
      pmFace.b = vertices[faces[i].b].index;
      pmFace.c = vertices[faces[i].c].index;

      // length of normal = area * 2
      area = area + normal_raw(points, &pmFace).length();
   }
   // return corrected area
   return area * 0.5f;
}

// ----------------------------------------------------------------------------
// boundary access
// ----------------------------------------------------------------------------
int CChart::BoundaryCount()
{
   return boundaries.Size();
}

// ----------------------------------------------------------------------------
// number of boundary edges
// ----------------------------------------------------------------------------
Edge *CChart::Boundaries()
{
   return &boundaries[0];
}

// ----------------------------------------------------------------------------
// reset everything but the seed
// ----------------------------------------------------------------------------
void CChart::Clear()
{
   boundaries.Clear();
   faces.Clear();
   vertices.Clear();
   clusterCount = 0;
}

// ----------------------------------------------------------------------------
// number of grows
// ----------------------------------------------------------------------------
int CChart::ClusterCount()
{
   return clusterCount;
}

// ----------------------------------------------------------------------------
// number of faces
// ----------------------------------------------------------------------------
int CChart::FaceCount()
{
   return faces.Size();
}

// ----------------------------------------------------------------------------
// face access
// ----------------------------------------------------------------------------
Face *CChart::Faces()
{
   return &faces[0];
}

// ----------------------------------------------------------------------------
// Grows the chart, and its boundary
// ----------------------------------------------------------------------------
void CChart::Grow(Vec3 &normal)
{
   chartNormal  += normal;
   averageNormal = chartNormal;
   //
   clusterCount++;
   //
   averageNormal.norm();
}

// ----------------------------------------------------------------------------
// returns the chart normal
// ----------------------------------------------------------------------------
Vec3 CChart::Normal()
{
   return averageNormal;
}

// ----------------------------------------------------------------------------
// returns the seed index of this chart
// ----------------------------------------------------------------------------
int CChart::Seed()
{
   return seed;
}

// ----------------------------------------------------------------------------
// returns the seed index of this chart
// ----------------------------------------------------------------------------
bool CChart::SetSeed(int seed, Vec3 &normal)
{
   this->seed          = seed;
   this->chartNormal   = normal;
   this->averageNormal = normal;
   //
   averageNormal.norm();

   if (seedHistory.Contains(seed) != -1)
      return true;

   // insert seed into history
   seedHistory.Insert(seed);
   return false;
}

// ----------------------------------------------------------------------------
// number of vertices
// ----------------------------------------------------------------------------
int CChart::VertexCount()
{
   return vertices.Size();
}

// ----------------------------------------------------------------------------
// vertex access
// ----------------------------------------------------------------------------
Vertex *CChart::Vertices()
{
   return &vertices[0];
}
} // end namespace