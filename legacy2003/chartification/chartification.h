//=============================================================================
//
// Chartification : Class Definition
//
// Aufteilung in TexturStücke mittels Chart-Growing
//
// $Author: Andor $
// $Date: 11.12.03 23:00 $
// $Revision: 18 $
// 
//=============================================================================
#ifndef CCHARTIFICATION_H
#define CCHARTIFICATION_H

// Local includes
#include "..\math\math3D.h"
#include "..\math\halfedge.h"
#include "..\math\heap.h"
#include "cchart.h"
#include "cchartboundary.h"
#include "cchartset.h"

// c++ includes
#include <vector>

namespace chartification {

class CChartification
{
   // original data structure reference
   int  numPoints;
   Vec3 *points;
   int  numFaces;
   Face *faces;

   // face normals
   Vec3 *faceNormals;
   Vec3 *faceMidpoints;

   // halfedge graph of the mesh
   int      numHalfedges;
   Halfedge *halfedges;

   // keep track of last face assigned
   int      lastface;

   // face to chart mapping
   CChart **chartMap;

   // Boundary structure
   CChartBoundary *boundary;

   // temporary info for vertex structure creation
   Vertex *vertices;
   Face   *vFaces;
   //
   int    numVertices;
   int    numVFaces;

public:
   // creates an CChartification object
   CChartification(int numPoints, Vec3 *points, int numFaces, Face *faces);

   // deletes an CChartification object
   ~CChartification();

   // starts the chartification process
   void Chartify(int numDesiredCharts, CChartSet *outCharts);

   // Visualization
   // allocates enough memory for the operations
   void Allocate(CChartSet  *charts, int *numPoints, Vec3 **points, rgbf **rgb_values, int *numFaces, Face **faces, int *numEdges, Edge **edges);
   void Visualization(int numDesiredCharts, int *numPoints, Vec3 **points, rgbf **rgb_values, int *numFaces, Face **faces, int *numEdges, Edge **edges);
   void TextureVisualization(int numDesiredCharts, int *numPoints, Vec3 **points, rgbf **rgb_values, int *numFaces, Face **faces, int *numEdges, Edge **edges);
   void Colorize(CChartSet  *charts, bool tex, int *numPoints, Vec3 **points, rgbf **rgb_values, int *numFaces, Face **faces, int *numEdges, Edge **edges);

   // creates a chart visualization
   void CreateVertices(int numDesiredCharts, int *outNumVertices, Vertex **outVertices, int *outNumFaces, Face **outFaces);

private:
   // region growing
   void ExpandCharts(CChartSet *charts);

   // repair heap
   void RepairHeap(CChart *chart, Heap<Halfedge *, double> *heap);

   // straighten the boundary
   bool Straighten(CChartBoundary *boundary);

   // returns the cost for moving from one face to the other
   double EdgeCost(Halfedge *h);
   double Distance(Halfedge *h);

   // adds a face and boundary information to chart
   void AddInfoToChart(CChart *chart, Halfedge *h);

   // creates the vertex info structure
   void CreateVertexStructure(CChartSet *charts, CChartBoundary *boundary);

   // searches is vertex is already allocated
   int SearchVertexId(Halfedge *h, bool *visited);

   // cuts along the boundary path
   void CutBoundaryPath(CChartBoundary *boundary);

   // Gets Faces from Chartset representation
   void GetFaces(CChartSet *charts, int *numFaces, Face **faces);

   // Gets Vertices from Chartset representation
   void GetVertices(CChartSet *charts, int *numVertex, Vertex **vertices);

   // seed search
   bool ComputeSeeds(CChartSet *charts, CChartBoundary *boundary);
   bool ApplySeeds(CChartSet *charts);
   void InitializeSeedSearch(Heap<Halfedge *, double> *heap, bool *visited, CChartBoundary *boundary);

}; // end CChartification

} // end namespace
#endif