/******************************************************************************
 *
 * MRPolyMesh
 *
 * Implementierung der Basisklasse MRPolyMesh
 *
 * Diese Klasse implementiert den Edge-Collapse Algorithmus
 * zur Vereinfachung von Polygon-Meshes, in den Abgeleiteten Klassen
 * müssen die Verfahren zur Fehlerbestimmung und die Generierung der
 * VertexPaare erstellt werden
 *
 ******************************************************************************/
#include "mrpolymesh.h"

/**
 * Erstellen einer PolyMesh Instanz
 */
MRPolyMesh::MRPolyMesh()
{
   this->numPoints = 0;
   this->numFaces = 0;
   this->numDstPoints = 0;

   this->points = NULL;
   this->faces  = NULL;
}

/**
 * Löschen der PolyMesh Instanz
 */
MRPolyMesh::~MRPolyMesh()
{
   if (points) delete[] points;
   if (faces) delete[] faces;
   if (invalidPoints) delete[] invalidPoints;

   if (faceIndices) delete[] faceIndices;
   if (pairIndices) delete[] pairIndices;
}

/**
 * initialisiert die Instanz
 */
void MRPolyMesh::initialize(int numPoints, vec3 *points, int numFaces, face *faces)
{
   this->numPoints = numPoints;
   this->numFaces = numFaces;
   this->numDstPoints = numPoints;

   this->points = new vec3[numPoints];
   this->faces  = new face[numFaces];

   memcpy(this->points, points, sizeof(vec3) * numPoints);
   memcpy(this->faces, faces, sizeof(face) * numFaces);

   invalidPoints = new bool[numPoints];
   for (int i=0; i<numPoints; i++) invalidPoints[i] = false;

   createVertexPairs(numPoints, points, numFaces, faces);
   //
   initializeErrorMetric(numPoints, points, numFaces, faces);
   initializeDependencies();
   initializeHeap();
}

/**
 * Reduziert das Modell um einen Punkt, indem zwei Punkte
 * durch einen Punkt ersetzt werden
 */
double MRPolyMesh::decimate()
{
   VertexPair *vertexPair = (VertexPair *) heap.pop();

   if (vertexPair)
   {
      // Eintragen der neuen Koordinaten
      updateCoordinates(points[vertexPair->v1], *vertexPair);

      // in allen Dreiecken v2-index durch v1-index ersetzen
      updateFaceIndices(vertexPair->v1, vertexPair->v2);
      addFaceIndices(vertexPair->v1, vertexPair->v2);

      // in allen Paaren v2-index durch v1-index ersetzen
      updatePairIndices(vertexPair->v1, vertexPair->v2);
      addPairIndices(vertexPair->v1, vertexPair->v2);

      // v2 als ungültig markieren
      markInvalid(vertexPair->v2);

      removeInvalidFaces();
      removeInvalidPairs();

      // Fehler für alle an v1 angrenzenden Paare neu berechnen
      updateError(vertexPair->v1);

      return vertexPair->heapValue;
   }
   return -1.0f;
}

/**
 * Erstellt die IndexListen
 * Punkt - Fläche,indem jede Fläche in die
 * Indexliste Ihrer drei Punkte a, b, c
 * eingetragen wird
 */
void MRPolyMesh::initializeFaceIndices()
{
   faceIndices = new int_set[numPoints];

   for (int i=0; i<numFaces; i++) {
      faceIndices[faces[i].a].insert(i);
      faceIndices[faces[i].b].insert(i);
      faceIndices[faces[i].c].insert(i);
   }
}

/**
 * Erstellt die IndexListen
 * Punkt - Paare, indem jedes Paar in die
 * Indexliste Ihrer zwei Punkte v1, v2
 * eingetragen wird
 */
void MRPolyMesh::initializePairIndices()
{
   pairIndices = new int_set[numPoints];

   for (VertexPairIterator it=beginVertexPair(); it!=endVertexPair(); it++) {
      pairIndices[it->v1].insert(pairIndex(it));
      pairIndices[it->v2].insert(pairIndex(it));
   }
}

/**
 * Erstellt die Abhängigkeiten Punkt - Fläche
 *                         bzw Punkt - VertexPaar
 */
void MRPolyMesh::initializeDependencies()
{
   initializeFaceIndices();
   initializePairIndices();
}

/**
 * Initialisierung des Heaps mit allen möglichen Vertex-Paaren
 */
void MRPolyMesh::initializeHeap()
{
   for (VertexPairIterator it=beginVertexPair(); it!=endVertexPair(); it++) {
      // vertexpaar und heap vorbelegen
      VertexPair *vertexPair = &(*it);
      //
      vertexPair->arrPos = pairIndex(vertexPair);
      vertexPair->heapValue = error(vertexPair);
      //
      heap.push(vertexPair);
   }
}

/**
 * ersetzt in allen an v1 und v2 angrenzenden Flächen
 * den index von v2 durch den Index von v1
 */
void MRPolyMesh::updateFaceIndices(int v1, int v2)
{
   updateFaceIndices(&faceIndices[v1], v1, v2);
   updateFaceIndices(&faceIndices[v2], v1, v2);
}

/**
 * ersetzt in allen an v1 und v2 angrenzenden Paaren
 * den index von v2 durch den Index von v1
 */
void MRPolyMesh::updatePairIndices(int v1, int v2)
{
   updatePairIndices(&pairIndices[v1], v1, v2);
}

/**
 * Entfernen der ungültigen Paare vom Heap
 */
void MRPolyMesh::removeInvalidPairs()
{
   for (int_iterator it=invalidPairs.begin(); it!=invalidPairs.end(); it++) {
      VertexPair *vertexPair = (VertexPair *) pair(*it);
      //
      pairIndices[vertexPair->v1].erase(*it);
      pairIndices[vertexPair->v2].erase(*it);
      //
      heap.remove(vertexPair);
   }
   invalidPairs.clear();
}

/**
 * Markiert einen Punkt als ungültig und
 * reduziert die Anzahl der gültigen Punkte
 */
void MRPolyMesh::markInvalid(int index)
{
   invalidPoints[index] = true;
   numDstPoints--;
}

/**
 * Entfernt eine ungültige Fläche, die Liste der Flächen
 * wird gepackt - d.h. die ungültige Fläche innerhalb der Liste wird
 * durch die letzte gültige Flächen ersetzt, die Indexlisten werden
 * entsprechend korrigiert
 */
void MRPolyMesh::removeInvalidFace(int index)
{
   // FaceIndex von allen vertex strukturen entfernen
   face *fp = &faces[index];
   //
   faceIndices[fp->a].erase(index);
   faceIndices[fp->b].erase(index);
   faceIndices[fp->c].erase(index);

   // punkte auf gültigkeit überprüfen
   if (faceIndices[fp->a].size() == 0) markInvalid(fp->a);
   if (faceIndices[fp->b].size() == 0) markInvalid(fp->b);
   if (faceIndices[fp->c].size() == 0) markInvalid(fp->c);

   numFaces--;

   // letzte face-struktur an die freigewordene stelle kopieren
   if (index < numFaces) {
      // copy face info to free place
      faces[index] = faces[numFaces];
      fp           = &faces[index];

      // index in facelisten korrigieren
      faceIndices[fp->a].erase(numFaces);
      faceIndices[fp->b].erase(numFaces);
      faceIndices[fp->c].erase(numFaces);
      //
      faceIndices[fp->a].insert(index);
      faceIndices[fp->b].insert(index);
      faceIndices[fp->c].insert(index);
   }
}

/**
 * Entfernt die ungültigen Flächen, die Liste der Flächen
 * wird gepackt - d.h. ungültige Flächen innerhalb der Liste werden
 * durch gültige Flächen am Ende der Liste ersetzt
 */
void MRPolyMesh::removeInvalidFaces()
{
   for (int_set::reverse_iterator rit=invalidFaces.rbegin(); rit!=invalidFaces.rend(); rit++) {
      removeInvalidFace(*rit);
   }
   invalidFaces.clear();
}

/**
 * ermittelt den Fehler für alle an durch index
 * bestimmten Punkt neu
 */
void MRPolyMesh::updateError(int index)
{
   int_set *indices = &pairIndices[index];

   for (int_iterator it=indices->begin(); it!=indices->end(); it++) {
      // Fehler für paar neu bestimmen
      VertexPair *vertexPair = pair(*it);
      vertexPair->heapValue = error(vertexPair);
      //
      heap.repair(vertexPair);
   }
}

