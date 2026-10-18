/******************************************************************************
 *
 * MRPolyMeshQE
 *
 * MRPolyMesh Implementation
 * Fehlerberechnung mittels Quadric Error Metric
 *
 ******************************************************************************/

// Local Includes
#include "mrpolymeshqe.h"


/*****************************************************************************
 * Hinzufügen eines Paares
 *****************************************************************************/
void MRPolyMeshQE::addPair(halfedge *h, bool boundary)
{
   // nächstes Paar in der Liste holen
   VertexPairQE *pair = &vertexPairs[numVertexPairs++];

   pair->v1       = h->end;
   pair->v2       = h->next->next->end;
   pair->boundary = boundary;
}

/**************************************************************************
 * Erzeugung aller möglichen Vertex-Paare
 **************************************************************************/
void MRPolyMeshQE::createVertexPairs(int numPoints, vec3 *points, int numFaces, face *faces)
{
   // HalfEdge Representation erzeugen
   int numHalfEdges;
   halfedge *halfedges;
   HalfEdgeMesh::create(numPoints, numFaces, faces, &numHalfEdges, &halfedges);

   // temp. marker erzeugen
   bool *visited = new bool[numHalfEdges];
   // initialisieren des markers
   for (int i=0; i<numHalfEdges; i++) {
      visited[i] = false;
   }
   // maximale anzahl pairs erzeugen
   vertexPairs = new VertexPairQE[numHalfEdges];
   numVertexPairs = 0;

   // iterate over each halfedge
   for (i=0; i<numHalfEdges; i++) {
      // check if we've been already there
      if (visited[i])
         continue;
      // fetch halfedge
      halfedge *h = &halfedges[i];
      // mark halfedge, symmetric halfedge
      visited[i] = true;
      if (h->sym) visited[h->sym - halfedges] = true;
      // add halfedge to possible contraction pairs
      addPair(h, h->sym == NULL);
   }

   // pairs reallokieren
   VertexPairQE *tmp = new VertexPairQE[numVertexPairs];
   memcpy(tmp, vertexPairs, sizeof(VertexPairQE) * numVertexPairs);
   // cleanup
   delete[] vertexPairs;
   delete[] halfedges;
   delete[] visited;
   //
   vertexPairs = tmp;
}

/**************************************************************************
 * Fehler Metrik
 **************************************************************************/
void MRPolyMeshQE::initializeErrorMetric(int numPoints, vec3 *points, int numFaces, face *faces)
{
   // Fehler matrizen erzeugen
   errorMetrics = new QError[numPoints];
   QError error;

   // Gleichungen initialisieren
   for (int i=0; i<numPoints; i++)
      errorMetrics[i].init();

   for (int i=0; i<numFaces; i++) {
      // Gleichung für eine Fläche bestimmen
      face *f = &faces[i];
      error.setPlane(points[f->a], points[f->b], points[f->c]);

      // Gleichungen für jeden Punkt aufaddieren
      errorMetrics[f->a].add(error);
      errorMetrics[f->b].add(error);
      errorMetrics[f->c].add(error);
   }
}
