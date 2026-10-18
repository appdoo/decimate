/**
 * class to reduce polygonal models
 * by keeping the appearance
 *
 * class is capable of producin a reduced model,
 * for a source model that has one texture map
 *
 * @author: Andor Greissl
 * @date:   22.03.2003
 */
#include "reduction.h"
#include "../math/halfedge.h"
#include <set>
#include <math.h>
#include <time.h>

/*
template <class T> inline void swap(T &a, T&b)
{
   T t;
   t = a;
   a = b;
   b = t;
}
*/
void inline replace(int_set *iset, int oldValue, int newValue)
{
   iset->erase(oldValue);
   iset->insert(newValue);
}



/*****************************************************************************
 *
 * Konstruktor für die Reduktion Klasse
 * 
 *****************************************************************************/
MeshSimplification::MeshSimplification(
                     vec3 *points,
                     int  numPoints,
                     face *faces,
                     int  numFaces
                    )
{
	// clear all
	this->points   = NULL;
	this->faces    = NULL;
	this->pairs    = NULL;
   this->quadrics = NULL;
   this->minHeap  = NULL;
   //
   this->numPoints   = numPoints;
   this->numFaces    = numFaces;
   this->numPairs    = 0;
   //
	// copy face data
	if (faces) {
		this->faces = (face *) memcpy(new face[numFaces], faces, numFaces * sizeof(face));
	}
   // create the connectivity structures
   vertices = new cvtx[numPoints];
   for (int i=0; i<numPoints; i++) {
      // copy geometric info and initialize sets
      vertices[i].move(points[i]);
      vertices[i].faceIndices = new int_set;
      vertices[i].p2pIndices  = new int_set;
   }

   // enter face connectivity info
   for (i=0; i<numFaces; i++) {
      vertices[faces[i].a].faceIndices->insert(i);
      vertices[faces[i].b].faceIndices->insert(i);
      vertices[faces[i].c].faceIndices->insert(i);
   }
   printf("numFaces: %d ", i);

   // create point connectivity info
   createPairs();
   createPairIndices();

   // Debug info
   printf("Moegliche Vertex-Paare: %d\n", numPairs);

   // Fehler matrizen erzeugen
   quadrics = new QError[numPoints];
   QError     planeError;

   for (i=0; i<numPoints; i++)
      quadrics[i].init();

   for (i=0; i<numFaces; i++) {
      face *f = &faces[i];
      planeError.setPlane(points[f->a], points[f->b], points[f->c]);

      quadrics[f->a].add(planeError);
      quadrics[f->b].add(planeError);
      quadrics[f->c].add(planeError);
   }

   // Heap füllen
   minHeap  = new MinHeap(numPairs);

   for (i=0; i<numPairs; i++) {
      cvtxpair &pair = pairs[i];

      getQuadricError(pair);
      minHeap->add(&pair);
   }
}



/*****************************************************************************
 * Hinzufügen eines Paares
 * 
 *****************************************************************************/
void MeshSimplification::addPair(halfedge *h, bool boundary)
{
   // nächstes Paar in der Liste holen
   cvtxpair *pair = &pairs[numPairs++];

   pair->v1     = h->end;
   pair->v2     = h->next->next->end;
   pair->v3     = h->next->end;
   pair->boundary = boundary;
   pair->arrPos = numPairs - 1;
}



/*****************************************************************************
 * Implementierung eines Algorithmus
 * zum Auffinden aller Kanten innerhalb des
 * Meshes
 *
 * Algorithmus:
 *  Einführung eines Markers pro Punkt
 *
 * Für jeden Punkt P
 *     - Füge alle von P ausgehenden Kanten, die zu nicht markierten
 *       Punkten führen, in die Liste aller Kanten ein
 *     
 *     - Markiere den Punkt
 *
 *****************************************************************************/
void MeshSimplification::createPairs()
{
   // get a halfedge representation of the mesh
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
   pairs    = new cvtxpair[numHalfEdges];
   numPairs = 0;

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
   cvtxpair *tmp = new cvtxpair[numPairs];
   memcpy(tmp, pairs, sizeof(cvtxpair) * numPairs);
   delete[] pairs;
   //
   pairs = tmp;
}
/*
void MeshSimplification::createPairs()
{
   // temp. marker erzeugen
   bool *visited = new bool[numPoints];
   // initialisieren des markers
   for (int i=0; i<numPoints; i++) {
      visited[i] = false;
   }
   // maximale anzahl pairs erzeugen
   pairs    = new cvtxpair[numFaces * 3];
   numPairs = 0;
   int_set pset;
   int_set::iterator it;

   // Suche über alle Punkte
   for (i=0; i<numPoints; i++) {
      // for all points
      cvtx &vtxRef = vertices[i];
      pset.clear();

      for (it=vtxRef.faceIndices->begin(); it!=vtxRef.faceIndices->end(); it++) {
         // get next face
         face &faceRef = faces[*it];

         if (faceRef.a!=i  &&  !visited[faceRef.a]) pset.insert(faceRef.a);
         if (faceRef.b!=i  &&  !visited[faceRef.b]) pset.insert(faceRef.b);
         if (faceRef.c!=i  &&  !visited[faceRef.c]) pset.insert(faceRef.c);
      }
      
      for (it=pset.begin(); it!=pset.end(); it++) {
         addPair(i, *it);
      }
      visited[i] = true;
   }
   delete[] visited;

   // pairs reallokieren
   cvtxpair *tmp = new cvtxpair[numPairs];
   memcpy(tmp, pairs, sizeof(cvtxpair) * numPairs);
   delete[] pairs;

   pairs = tmp;
}
*/


/*****************************************************************************
 * Indexlisten für Paare erzeugen
 *****************************************************************************/
void MeshSimplification::createPairIndices()
{
   // Indexliste erzeugen
   for (int i=0; i<numPairs; i++) {
      cvtx &v1 = vertices[pairs[i].v1];
      cvtx &v2 = vertices[pairs[i].v2];

      v1.p2pIndices->insert(i);
      v2.p2pIndices->insert(i);
   }

}



/*****************************************************************************
 * Aufraeumen vor dem Freigeben des Objekts
 *****************************************************************************/
MeshSimplification::~MeshSimplification()
{
   // connectivity info löschen
   if (vertices) 
   {
      for (int i=0; i<numPoints; i++) 
      {
         if (vertices[i].faceIndices) delete vertices[i].faceIndices;
         if (vertices[i].p2pIndices)  delete vertices[i].p2pIndices;
      }
      delete[] vertices;
   }
	// Speicher freigeben
	if (points)   delete[] points;
	if (faces)    delete[] faces;
	if (pairs)    delete[] pairs;
   if (quadrics) delete[] quadrics;
   if (minHeap)  delete   minHeap;
}


/*****************************************************************************
 * funktion um den Fehler anhand der quadrierten
 * Abstaende von den Flaechen zu ermitteln
 *****************************************************************************/
void MeshSimplification::getQuadricError(cvtxpair &vp)
{
   vec3 &v  = vp.v;
   cvtx &v1 = vertices[vp.v1];
   cvtx &v2 = vertices[vp.v2];
   cvtx &v3 = vertices[vp.v3];

   QError qerror = quadrics[vp.v1];
   qerror.add(quadrics[vp.v2]);

   // check if pair lies on a border
/*   if (vp.boundary) {
      vec3 a(v1.x, v1.y, v1.z);
      vec3 b(v2.x, v2.y, v2.z);
      vec3 c(v3.x, v3.y, v3.z);
      vec3 dist = (c - a) * 5000.0;
      //
      c = ((c - a) ^ (b - a)) + a;
      //vec3 nrm = ((b - a) ^ (c - a)) * 5.0f;
      QError error;
      error.setPlane(a , b, c);
      for (int i=0; i<1000; i++)
         qerror.add(error);
   }*/
   // try to find an optimal point to replace the pair
   if (!qerror.optimize(v))
   {
      if (vp.boundary) printf("\\");
      v.x = v1.x + (v2.x - v1.x) * 0.5f;
      v.y = v1.y + (v2.y - v1.y) * 0.5f;
      v.z = v1.z + (v2.z - v1.z) * 0.5f;
   }

   vp.heapValue = qerror.getError(v);
}


/*****************************************************************************
 * methode zum reduzieren bis zu einer gewissen 
 * polygonanzahl
 *****************************************************************************/
void MeshSimplification::decimate(int border, int *numCoord, vec3 **coords, int *numFaces, face **faces)
{
   // print debug information
   checkStructure();
   printf("\n Progress [");
   // variables for progress information
   double overall = numPoints - border;
   int    mark    = 0;

   // decimate til the desired number of points is reached
   while (numPoints  > border) {
      if (decimate() == -1) {
         printf("Error while contracting....\n");
         break;
      }
      // print progress information
      double now = overall - (numPoints - border);
      now = (now * 10.0) / overall;
      if ((int) now > mark) {
         mark++;
         printf(".");
      }
   }
   // print debug information
   printf("] done\n\n");
   // Do some checks on the created structure
   checkStructure();
   
   // copy coordinate information
   *numCoord = numPoints;
   *coords = new vec3[numPoints];
   for (int i=0; i<numPoints; i++) {
      (*coords)[i].x = vertices[i].x;
      (*coords)[i].y = vertices[i].y;
      (*coords)[i].z = vertices[i].z;
   }
   // copy face structure
   *numFaces = this->numFaces;
   *faces = new face[this->numFaces];
   memcpy(*faces, this->faces, sizeof(face) * this->numFaces);
}



/*****************************************************************************
 * decimate
 *
 * die Grundfunktionalität, um die Anzahl der Polygone zu reduzieren
 * 2 Punkte werden aus dem Modell entfernt und durch einen neuen Punkt ersetzt
 *
 * (v1, v2) -> v
 *
 * Dabei muss folgendes geschehen:
 *
 * entfernen von den zwei Punkten, und Einfügen eines neuen Punktes
 * danach muss die komplette Struktur korrigiert werden, d.h.
 *
 * entfernen aller fehlerhaften Flächen und Ecken
 *
 * zur effizienten Implementierung eines solchen Algorithmus
 * muss eine Liste über freie Punkte, Ecken, und Flächen geführt werden
 * dies ist möglich, da die Anzahl der Punkte ja nicht ansteigt
 *****************************************************************************/
int MeshSimplification::decimate()
{
   // vertex-paar vom Heap entfernen
   cvtxpair *vp = (cvtxpair *) minHeap->remove(0);
   // check for error
   if (!vp) {
      printf("Keine Vertex - Paare mehr...\n");
      return -1;
   }

   // Dateninhalte kopieren
   int v1 = vp->v1;
   int v2 = vp->v2;

   // 1. Paar aus den Indexlisten löschen
   cvtx &vtx1 = vertices[v1];
   cvtx &vtx2 = vertices[v2];
   vtx1.p2pIndices->erase(vp->arrPos);
   vtx2.p2pIndices->erase(vp->arrPos);


   // 2. neue Koordinaten eintragen
   vertices[v1].move(vp->v);

   int_set::iterator it;

   // 3. Flaechen und Paar indices addieren
   for (it=vtx2.faceIndices->begin(); it!=vtx2.faceIndices->end(); ++it)
      vtx1.faceIndices->insert(*it);

   // 4. Index reparieren und liste der degenerierte Flächen erstellen
   invalidFaces.clear();

   for (it=vtx1.faceIndices->begin(); it!=vtx1.faceIndices->end(); ++it) {
      face &f   = faces[*it];

      // Flaeche entfernen
      if ((f.a == v1 || f.b == v1 || f.c == v1) && (f.a == v2 || f.b == v2 || f.c == v2))
      {
         invalidFaces.insert(*it);
      }
      // correct face indices
      if (f.a == v2) f.a = v1;
      if (f.b == v2) f.b = v1;
      if (f.c == v2) f.c = v1;
   }
   // Fehlermetrik korrigieren
   quadrics[v1].add(quadrics[v2]);

   // Paare korrigieren
   for (it=vtx2.p2pIndices->begin(); it!=vtx2.p2pIndices->end(); ++it)
     vtx1.p2pIndices->insert(*it);

   invalidPairs.clear();
   connections.clear();

   for (it=vtx1.p2pIndices->begin(); it!=vtx1.p2pIndices->end(); ++it)
   {
      cvtxpair *pair = &pairs[*it];
      if (vp->boundary)
         pair->boundary = true;
      // index in den paarungen umstellen
      if (pair->v1 == v2) pair->v1 = v1;
      if (pair->v2 == v2) pair->v2 = v1;

      if (pair->v1 != v1)
         swap(pair->v1, pair->v2);

      if (pair->v1 == pair->v2 || (pair->v1!=v1 && pair->v2!=v1))
         printf("invalid pair\n");

      if (connections.find(pair->v2)!=connections.end())
      {
         invalidPairs.insert(*it);
      }
      else
      {
         connections.insert(pair->v2);

         getQuadricError(*pair);
         minHeap->repair(pair);
      }
   }
   try{
   // remove invalid faces, process highest index first
   for (int_set::reverse_iterator rit=invalidFaces.rbegin(); rit!=invalidFaces.rend(); rit++) {
         removeFace(*rit);
   }
   // remove invalid pairs
   for (it=invalidPairs.begin(); it!=invalidPairs.end(); it++) {
      removePair(*it);
   }

   // freien vertex löschen
   vtx2.faceIndices->clear();
   vtx2.p2pIndices->clear();
   removeVertex(v2);
   }
   catch (...)
   {
      printf("error in MeshSimplification::contract\n");
   }
   // neue Anzahl der Punkte zurückliefern
   return numPoints;
}



/*****************************************************************************
 * löschen einer Vertex struktur,
 * die struktur wird gelöscht und in die liste der 
 * 
 * und der freie platz durch den letzten vertex in der liste
 * ersetzt, die Indexlisten werden korrigiert
 *****************************************************************************/
void MeshSimplification::removeVertex(int index)
{
   // Freigabe der indices
   if (vertices[index].faceIndices) delete vertices[index].faceIndices;
   if (vertices[index].p2pIndices)  delete vertices[index].p2pIndices;
   memset(&vertices[index], 0, sizeof(cvtx));
   numPoints--;

   // nachrücken des obersten vertex
   if (index < numPoints)
   {
      vertices[index] = vertices[numPoints];
      cvtx &vtxRef    = vertices[index];

      // Fehlermetrik korrigieren
      quadrics[index] = quadrics[numPoints];
      int_set::iterator it;

      for (it=vtxRef.faceIndices->begin(); it!=vtxRef.faceIndices->end(); ++it)
      {
         face &faceRef = faces[*it];

         if (faceRef.a == numPoints) faceRef.a = index;
         if (faceRef.b == numPoints) faceRef.b = index;
         if (faceRef.c == numPoints) faceRef.c = index;
      }

      for (it=vtxRef.p2pIndices->begin(); it!=vtxRef.p2pIndices->end(); ++it)
      {
         cvtxpair &pairRef = pairs[*it];

         if (pairRef.v1 == numPoints) pairRef.v1 = index;
         if (pairRef.v2 == numPoints) pairRef.v2 = index;
      }
   }
   // löschen der alten struktur
   memset(&vertices[numPoints], 0, sizeof(cvtx));
}



/*****************************************************************************
 * prüfen der erzeugten Abhängigkeiten
 *****************************************************************************/
void MeshSimplification::checkStructure()
{
   // print debug info
   printf("checking...");
   // create a reference table for points
   bool *referenced = new bool[numPoints];
   for (int i=0; i<numPoints; i++) {
      referenced[i] = false;
   }
   // check all faces
   for (int i=0; i<numFaces; i++) {
      // get current face reference
      face &f = faces[i];

      if ((f.a >= numPoints  ||  f.b >= numPoints  ||  f.c >= numPoints) ||
          (f.a < 0  ||  f.b < 0  ||  f.c < 0) ||
          (f.a==f.b || f.a==f.c || f.b==f.c)
          )
      {
         printf("removeVertex failed on face %d\n", i);
         printf("  Face.a %d\n", f.a);
         printf("  Face.b %d\n", f.b);
         printf("  Face.c %d\n\n", f.c);

      }
      else {
         // enter point reference
         referenced[f.a] = true;
         referenced[f.b] = true;
         referenced[f.c] = true;
      }
   }
   int unreferenced = 0;
   // print unreferenced points info
   for (int i=0; i<numPoints; i++) {
      if (!referenced[i]) {
         printf("\runreferenced point %d removing from mesh", i);
         removeVertex(i);
         unreferenced++;
      }
   }
   if (unreferenced > 0)
      printf("\rRemoved %d unreferenced points\n", unreferenced);
   // cleanup
   delete[] referenced;
}



/*****************************************************************************
 * löscht eine cvtxpair struktur aus dem pairs[] array,
 * der frei gewordene Platz wird durch das letzte Element der
 * Liste ersetzt
 * die indexlisten werden entsprechend korrigiert
 *****************************************************************************/
void MeshSimplification::removePair(int pairIndex)
{
   cvtx &v1 = vertices[pairs[pairIndex].v1];
   cvtx &v2 = vertices[pairs[pairIndex].v2];

   v1.p2pIndices->erase(pairIndex);
   v2.p2pIndices->erase(pairIndex);

   // Struktur vom Heap löschen
   minHeap->remove(&pairs[pairIndex]);
}



/******************************************************************************
 * löscht eine face struktur aus dem faces[] array,
 * der frei gewordene Platz wird durch das letzte Element der
 * Liste ersetzt
 * die indexlisten werden entsprechend korrigiert
/******************************************************************************/
void MeshSimplification::removeFace(int index)
{
   // FaceIndex von allen vertex strukturen entfernen
   face *fp = &faces[index];
   //
   vertices[fp->a].faceIndices->erase(index);
   vertices[fp->b].faceIndices->erase(index);
   vertices[fp->c].faceIndices->erase(index);

   // TODO: Debug Info
   if (vertices[fp->a].faceIndices->size() == 0) printf("punkt verloren\n");
   if (vertices[fp->b].faceIndices->size() == 0) printf("punkt verloren\n");
   if (vertices[fp->c].faceIndices->size() == 0) printf("punkt verloren\n");

   numFaces--;
   //invalid_faces[numFaces] = true;

   // letzte face-struktur an die freigewordene stelle kopieren
   if (index < numFaces) {
      // copy face info to free place
      faces[index] = faces[numFaces];
      fp             = &faces[index];

      // index in facelisten korrigieren
      vertices[fp->a].faceIndices->erase(numFaces);
      vertices[fp->b].faceIndices->erase(numFaces);
      vertices[fp->c].faceIndices->erase(numFaces);
      //
      vertices[fp->a].faceIndices->insert(index);
      vertices[fp->b].faceIndices->insert(index);
      vertices[fp->c].faceIndices->insert(index);
   }
   memset(&faces[numFaces], 0, sizeof(face));
}