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
#include <set>




/**
 * create a reduction class
 * for the given parameters
 *
 * 
 */
Reduction::Reduction(vec3 *points, int numPoints, face *faces, int numFaces, tex2 *textures, int numTextures)
{
	// clear all
	this->points   = NULL;
	this->faces    = NULL;
	this->pairs    = NULL;
	this->textures = NULL;
   //
   this->numPoints   = 0;
   this->numFaces    = 0;
   this->numPairs    = 0;
   this->numTextures = 0;
		
	// points data
	if (points)
	{
		this->points = (vec3 *) memcpy(new vec3[numPoints], points, numPoints * sizeof(vec3));
		this->numPoints = numPoints;
	}

	// face data
	if (faces)
	{
		this->faces = (face *) memcpy(new face[numFaces], faces, numFaces * sizeof(face));
		this->numFaces = numFaces;
	}

	// texture data
	if (textures)
	{
		this->textures    = (tex2 *) memcpy(new tex2[numTextures], textures, numTextures * sizeof(tex2));
		this->numTextures = numTextures;
	}


   // create the connectivity structures
   vertices = new cvtx[numPoints];
   memset(vertices, 0, sizeof(cvtx) * numPoints);

   // step 1. get memory usage per vertex, edge and face lists
   for (int i=0; i<numFaces; i++) {
      vertices[faces[i].a].numFaces++;
      vertices[faces[i].b].numFaces++;
      vertices[faces[i].c].numFaces++;
   }

   // step 2: allocate memory
   for (i=0; i<numPoints; i++) {
      vertices[i].faceIndices = new int[vertices[i].numFaces];
   }

   // step 3 enter face references
   int *cnt = new int[numPoints];
   memset(cnt, 0, sizeof(int) * numPoints);

   for (i=0; i<numFaces; i++) {
      vertices[faces[i].a].faceIndices[cnt[faces[i].a]++] = i;
      vertices[faces[i].b].faceIndices[cnt[faces[i].b]++] = i;
      vertices[faces[i].c].faceIndices[cnt[faces[i].c]++] = i;
   }
   delete[] cnt;

   std::set<int> p2p;
   numPairs = 0;

   // create point connectivity info
   for (i=0; i<numPoints; i++) {

      // information zurücksetzen
      p2p.clear();
      cvtx &v = vertices[i];

      for (int j=0; j<v.numFaces; j++) {
         // jedes Dreieck bearbeiten
         face &f = faces[v.faceIndices[j]];

         if (f.a > i) p2p.insert(f.a);
         if (f.b > i) p2p.insert(f.b);
         if (f.c > i) p2p.insert(f.c);
      }

      // informationen kopieren
      v.numP2P = p2p.size();
      v.p2pIndices = new int[v.numP2P];
      numPairs += v.numP2P;

      std::set<int>::const_iterator it;
      int *pIndices = v.p2pIndices;
      //
      for (it=p2p.begin(); it!=p2p.end(); it++) {
         *pIndices++ = *it;
      }
   }
   // Debug info
   printf("Moegliche Vertex-Paare: %d\n", numPairs);

   // Paare erzeugen
   pairs = new scvtxpair[numPairs];
   memset(pairs, 0, sizeof(scvtxpair) * numPairs);
   int index = 0;

   for (i=0; i<numPoints; i++) {
      scvtx &v = vertices[i];

      for (int j=0; j<v.numP2P; j++) {
         scvtxpair &vp = pairs[index++];

         vp.v1 = &vertices[i];
         vp.v2 = &vertices[v.p2pIndices[j]];
      }
   }
}


/**
 * Aufraeumen vor dem Freigeben des Objekts
 */
Reduction::~Reduction()
{
   if (vertices) {
      for (int i=0; i<numPoints; i++) {
         if (vertices[i].faceIndices) delete[] vertices[i].faceIndices;
         if (vertices[i].p2pIndices)  delete[] vertices[i].p2pIndices;
      }
      delete[] vertices;
   }

	// Speicher freigeben
	if (points)   delete[] points;
	if (faces)    delete[] faces;
	if (pairs)    delete[] pairs;
	if (textures) delete[] textures;

   // clear all
   vertices = NULL;
	points   = NULL;
	faces    = NULL;
	pairs    = NULL;
	textures = NULL;
   //
   numPoints   = 0;
   numFaces    = 0;
   numPairs    = 0;
   numTextures = 0;
}


/**
 * contract
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
 */
void Reduction::contract(int i1, int i2, cvtx &v)
{
   // lokale Kopie der ursp. punkte erstellen
   cvtx &v1 = vertices[i1];
   cvtx &v2 = vertices[i2];

   int numFaces = 0;
   int *faceIndices = new int[v1.numFaces + v2.numFaces];

   // index i2 in allen flaechen und ecken auf i1 setzen
   for (int i=0; i<v2.numFaces; i++) {

      face &f = faces[v2.faceIndices[i]];

      // pruefen ob das polygon degeneriert
      if (f.a==i1 || f.b==i1 || f.c==i1) {
         removeFace(f);
      }
      else {
         // index umstellen
         if (f.a == i2) f.a = i1;
         if (f.b == i2) f.b = i1;
         if (f.c == i2) f.c = i1;

         faceIndices[numFaces++] = v2.faceIndices[i];
      }
   }

   // derselbe vorgang für v1
   for (i=0; i<v1.numFaces; i++) {

      face &f = faces[v1.faceIndices[i]];

      // pruefen ob das polygon degeneriert
      if (f.a==i2 || f.b==i2 || f.c==i2) {
         removeFace(f);
      }
      else {
         // index umstellen
         if (f.a == i2) f.a = i1;
         if (f.b == i2) f.b = i1;
         if (f.c == i2) f.c = i1;

         faceIndices[numFaces++] = v1.faceIndices[i];
      }
   }
   v.numFaces    = numFaces;
   v.faceIndices = new int[numFaces];
   memcpy(v.faceIndices, faceIndices, sizeof(int) * numFaces);
   delete[] faceIndices;

   removeVertex(v2);
   replaceVertex(v1, v);
}


/**
 * speicherfreigabe für einen Vertex
 */
void Reduction::deleteVertex(cvtx &v)
{
   if (v.p2pIndices)  delete[] v.p2pIndices;
   if (v.faceIndices) delete[] v.faceIndices;
   //
   if (v.quadric)     delete[] v.quadric;

   memset(&v, 0, sizeof(cvtx));
}


/**
 * ersetzen einer Vertex struktur,
 * die alte struktur wird gelöscht und eine
 * binäre kopie wird an diese stelle geschrieben
 */
void Reduction::replaceVertex(cvtx &dst, cvtx &src)
{
   // löschen der alten struktur
   deleteVertex(dst);

   // kopieren der neuen struktur
   memcpy(&dst, &src, sizeof(cvtx));
}


/**
 * löschen einer Vertex struktur,
 * die struktur wird gelöscht und in die liste der 
 * freien vertex - strukturen eingehängt
 */
void Reduction::removeVertex(cvtx &v)
{
   // löschen der alten struktur
   deleteVertex(v);

   // kopieren der neuen struktur
   
}


/**
 * löschen einer Face struktur,
 * die struktur wird gelöscht und in die liste der 
 * freien face - strukturen eingehängt
 */
void Reduction::removeFace(face &f)
{
   f.a = -1;
   f.b = -1;
   f.c = -1;
}