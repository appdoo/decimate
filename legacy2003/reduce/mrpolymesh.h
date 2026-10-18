/******************************************************************************
 *
 * MRPolyMesh
 *
 * Interface for a MultiResolution capable Polygon Mesh
 *
 ******************************************************************************/

#ifndef MRPOLYMESH
#define MRPOLYMESH

#include "..\math\math3D.h"
#include "associativeheap.h"

/**
 *
 * v1, v2 - Point Indices
 */
typedef struct sVertexPair : public Heapable
{
   int v1, v2;
} 
VertexPair;

typedef VertexPair *VertexPairIterator;

class MRPolyMesh
{
private:
   AssociativeHeap heap;

   int numFaces;
   int numPoints;

   int numDstPoints;

   // Flächen und Punkte des Polygon-Meshes
   face *faces;
   vec3 *points;

   // Indexlisten für Flächen und Paare
   int_set *faceIndices;
   int_set *pairIndices;

   // sets für ungültige paare, flächen und punkte
   int_set invalidFaces;
   int_set invalidPairs;
   bool *invalidPoints;

public:
   /**
    * Erstellen einer PolyMesh Instanz
    */
   MRPolyMesh();

   /**
    * Löschen der PolyMesh Instanz
    */
   virtual ~MRPolyMesh();

   /**
    * Reduziert das Modell um einen Punkt, indem zwei Punkte
    * durch einen Punkt ersetzt werden
    */
   double decimate();

   /**
    * initialisiert die Instanz
    */
   void initialize(int numPoints, vec3 *points, int numFaces, face *faces);

private:
   /**
    * Erstellt die IndexListen
    * Punkt - Fläche,indem jede Fläche in die
    * Indexliste Ihrer drei Punkte a, b, c
    * eingetragen wird
    */
   void initializeFaceIndices();

   /**
    * Erstellt die IndexListen
    * Punkt - Paare, indem jedes Paar in die
    * Indexliste Ihrer zwei Punkte v1, v2
    * eingetragen wird
    */
   void initializePairIndices();

   /**
    * Erstellt die Abhängigkeiten Punkt - Fläche
    *                         bzw Punkt - VertexPaar
    */
   void initializeDependencies();

   /**
    * Initialisierung des Heaps mit allen möglichen Vertex-Paaren
    */
   void initializeHeap();

   /**
    * ersetzt in durch index bestimmten Fläche
    * den index von v2 durch den Index von v1
    */
   inline void updateFaceIndices(int index, int v1, int v2)
   {
      face &f = faces[index];

      if ((f.a==v1 || f.b==v1 || f.c==v1) && (f.a==v2 || f.b==v2 || f.c==v2)) {
         // set to invalid state
         invalidFaces.insert(index);
      }
      // change index
      if (f.a == v2) f.a = v1;
      if (f.b == v2) f.b = v1;
      if (f.c == v2) f.c = v1;
   }

   /**
    * ersetzt in allen von faceIndices angegebenen Flächen
    * den index von v2 durch den Index von v1
    */
   inline void updateFaceIndices(int_set *facesIndices, int v1, int v2)
   {
      for (int_iterator it=faceIndices->begin(); it!=faceIndices->end(); it++) {
         updateFaceIndices(*it, v1, v2);
      }
   }

   /**
    * ersetzt in allen an v1 und v2 angrenzenden Flächen
    * den index von v2 durch den Index von v1
    */
   void updateFaceIndices(int v1, int v2);

   /**
    * ersetzt in dem durch index bestimmten Paar
    * den index von v2 durch den Index von v1
    */
   inline void updatePairIndices(int index, int v1, int v2)
   {
      VertexPair *p = pair(index);

      if ((p->v1==v1 || p->v2==v1) && (p->v1==v2 || p->v2==v2)) {
         // set to invalid state
         invalidPairs.insert(index);
      }
      // change index
      if (p->v1 == v2) p->v1 = v1;
      if (p->v2 == v2) p->v2 = v1;
   }

   /**
    * ersetzt in allen von faceIndices angegebenen Paaren
    * den index von v2 durch den Index von v1
    */
   inline void updatePairIndices(int_set *pairIndices, int v1, int v2)
   {
      for (int_iterator it=pairIndices->begin(); it!=pairIndices->end(); it++) {
         updatePairIndices(*it, v1, v2);
      }
   }

   /**
    * ersetzt in allen an v1 und v2 angrenzenden Paaren
    * den index von v2 durch den Index von v1
    */
   void updatePairIndices(int v1, int v2);

   /**
    * fügt die von v2 referenzierten Flächen
    * zur indexliste von v1 hinzu
    */
   inline void addFaceIndices(int v1, int v2)
   {
      addIntSets(&faceIndices[v1], &faceIndices[v2]);
   }

   /**
    * fügt die von v2 referenzierten Flächen
    * zur indexliste von v1 hinzu
    */
   inline void addPairIndices(int v1, int v2)
   {
      addIntSets(&pairIndices[v1], &pairIndices[v2]);
   }

   /**
    * fügt alle Elemente von src zu dest hinzu
    */
   inline void addIntSets(int_set *dest, int_set *src)
   {
      for (int_iterator it=src->begin(); it!=src->end(); ++it)
         dest->insert(*it);
   }

   /**
    * Markiert einen Punkt als ungültig und
    * reduziert die Anzahl der gültigen Punkte
    */
   void markInvalid(int index);

   /**
    * Entfernen der ungültigen Paare,
    * Flächen und Punkte
    */
   void removeInvalidFace(int index);
   void removeInvalidPairs();
   void removeInvalidFaces();
   void removeInvalidPoints();

   /**
    * ermittelt den Fehler für alle an durch index
    * bestimmten Punkt neu
    */
   void updateError(int index);

protected:
   /**
    * accessor für punktkoordinaten
    */
   inline const vec3 &getPoint(int index)
   {
      return points[index];
   }

   /**
    * Erzeugung aller möglichen Vertex-Paare
    **/
   virtual void createVertexPairs(int numPoints, vec3 *points, int numFaces, face *faces) = 0;

   /**
    * Iteratoren über die Vertex-Paare
    **/
   virtual VertexPairIterator beginVertexPair() = 0;
   virtual VertexPairIterator endVertexPair() = 0;

   /**
    * Eintragen von neuen Coordinaten in den
    * referenzierten Punkt
    */
   virtual void updateCoordinates(vec3 &v, VertexPair &vertexPair) = 0;

   /**
    * Zugriff auf Vertex-Paare für die Basis-Klasse
    */
   virtual int pairIndex(VertexPair *vertexPair) = 0;
   virtual VertexPair *pair(int index) = 0;

   /**
    * Fehler Metrik
    */
   virtual void initializeErrorMetric(int numPoints, vec3 *points, int numFaces, face *faces) = 0;
   virtual double error(VertexPair *vertexPair) = 0;
};

#endif






