/******************************************************************************
 *
 * MRPolyMeshQE
 *
 * MRPolyMesh Implementation
 * Fehlerberechnung mittels Quadric Error Metric
 *
 ******************************************************************************/
#ifndef MRPOLYMESHQE_H
#define MRPOLYMESHQE_H

// Local Includes
#include "mrpolymesh.h"
#include "qerror.h"
#include "../math/halfedge.h"

// Vertex Paar für Quadric Error Metric
typedef struct sVertexPairQE : public VertexPair
{
   vec3 np;
   bool boundary;
} VertexPairQE;


class MRPolyMeshQE : public MRPolyMesh
{
private:
   // Liste der VertexPaare
   VertexPairQE *vertexPairs;
   int numVertexPairs;

   // Error Metrix
   QError *errorMetrics;
   int numErrorMetrics;

public:
   /**
    * Erstellen einer PolyMesh Instanz
    */
   MRPolyMeshQE(int numPoints, vec3 *points, int numFaces, face *faces)
   {
      initialize(numPoints, points, numFaces, faces);
   }

   /**
    * Löschen der PolyMesh Instanz
    */
   virtual ~MRPolyMeshQE()
   {
      if (vertexPairs) delete[] vertexPairs;
      if (errorMetrics) delete[] errorMetrics;
   }

protected:
   /*****************************************************************************
    * Hinzufügen eines Paares
    * 
    *****************************************************************************/
   void addPair(halfedge *h, bool boundary);

   /**************************************************************************
    * Erzeugung aller möglichen Vertex-Paare
    **************************************************************************/
   virtual void createVertexPairs(int numPoints, vec3 *points, int numFaces, face *faces);

   /**************************************************************************
    * Iteratoren über die Vertex-Paare
    **************************************************************************/
   virtual VertexPairIterator beginVertexPair()
   {
      return &vertexPairs[0];
   }

   virtual VertexPairIterator endVertexPair()
   {
      return &vertexPairs[numVertexPairs];
   }

   /**************************************************************************
    * Eintragen von neuen Coordinaten in den
    * referenzierten Punkt
    **************************************************************************/
   virtual void updateCoordinates(vec3 &v, VertexPair &vertexPair)
   {
      // typecast ist innerhalb der Klasse sicher
      VertexPairQE *vertexPairQE = (VertexPairQE *)(&vertexPair);
      v = vertexPairQE->np;
   }

   /**************************************************************************
    * Zugriff auf Vertex-Paare für die Basis-Klasse
    **************************************************************************/
   virtual int pairIndex(VertexPair *vertexPair)
   {
      return ((VertexPairQE *) vertexPair) - &vertexPairs[0];
   }
   
   virtual VertexPair *pair(int index)
   {
      return (VertexPair *)&vertexPairs[index];
   }

   /**************************************************************************
    * Fehler Metrik
    **************************************************************************/
   virtual void initializeErrorMetric(int numPoints, vec3 *points, int numFaces, face *faces);

   /**************************************************************************
    * Liefert den Fehler für ein
    * Vertex-Paar
    **************************************************************************/
   virtual double error(VertexPair *vertexPair)
   {
      QError error;
      VertexPairQE *pair = (VertexPairQE *)(vertexPair);

      error.init(errorMetrics[pair->v1], errorMetrics[pair->v2]);

      if (!error.optimize(pair->np))
      {
         // Punkt mitteln
         pair->np = (getPoint(pair->v1) + getPoint(pair->v2)) * 0.5f;
      }
      return error.getError(pair->np);
   }
};

#endif