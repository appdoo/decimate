//=============================================================================
//
// Grundgerüst für einen Edge-Collapse Algorithmus
//
// $Author: Andor $
// $Date: 7.12.03 13:06 $
// $Revision: 4 $
// 
//=============================================================================

#ifndef EDGECOLLAPSE_H
#define EDGECOLLAPSE_H


// Local Includes
#include "ivertexpaircollection.h"
#include "ierrorsystem.h"
#include "..\math\cset.h"

// begin namespace
namespace decimate {

// ============================================================================
// EdgeCollapse Class Definition
//
//
class CEdgeCollapse
{
   IErrorSystem *errorSystem;
   IVertexPairCollection *vertexPairCollection;

   AssociativeHeap heap;

   int numFaces;
   int numPoints;

   int numValidPoints;

   // Flächen und Punkte des Polygon-Meshes
   Face *faces;
   Vec3 *points;

   // Indexlisten für Flächen und Paare
   CSet<int> *faceIndices;
   CSet<int> *pairIndices;

   // sets für ungültige paare, flächen und punkte
   CSet<int> invalidFaces;
   CSet<int> invalidPairs;

   bool *invalidPoints;

public:
   // creates an edgecollapse object
   CEdgeCollapse(int numPoints, const Vec3 *points, int numFaces, const Face *faces, IErrorSystem *errorSystem, IVertexPairCollection *vertexPairCollection);

   // deletes an edgecollapse object
   ~CEdgeCollapse();

   // decimates the data structure
   bool Decimate();
   void DecimateToPoints(int numDesiredPoints, int *numCoord, Vec3 **coords, int *numFaces, Face **faces);


private:
   // Setup data structures
   void Initialize();

   // Compress data structure
   void Compress();

   // Update
   void UpdateError(int index);

   // remove functions
   void RemoveInvalidPairs();
   void RemoveInvalidFaces();
   void RemoveInvalidFace(int index);
   void RemovePoint(int index);

   // remapping
   void RemapFaceIndices(CSet<int> &indices, int new_value, int old_value);

   // index functions
   void MarkInvalid(int index);

   void UpdateFaceIndices(int p1, int p2);
   void UpdateFaceIndices(int index, int p1, int p2);

   void UpdatePairIndices(int p1, int p2);
   void UpdatePairIndices(int index, int p1, int p2);
};


} // end namespace
#endif



