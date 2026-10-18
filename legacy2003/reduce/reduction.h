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
#ifndef REDUCTION_H
#define REDUCTION_H


#include "..\math\math3D.h"
#include "..\math\halfedge.h"
#include "qerror.h"
#include "minheap.h"
#include <set>


/******************************************************************************
 *
 * Zum vereinfachen eines 3D - Modells ist eine wesentlich komplexere
 * Struktur notwendig als zum simplen darstellen
 *
 * so kann es durchaus bei einem Modell vorkommen, dass ein und derselbe 
 * Punkt (d.h. dieselben Koordinatenwerte mehrfach vorhanden sind,
 * um sie z.b. mir unterschiedlichen Texturen zu verknüpfen
 *
 * beim entfernen, bzw. ersetzen eines Punktes muss dies natürlich 
 * berücksichtigt werden
 *
 ******************************************************************************/

class cvtx
{
public:
   float x;
   float y;
   float z;

   std::set<int> *faceIndices;
   std::set<int> *p2pIndices;

   void move(vec3 &v)
   {
     x = v.x;
     y = v.y;
     z = v.z;
   }
};


/******************************************************************************
 *
 * Struktur für eine moegliche Kontraktion von zwei Punkten
 *
 ******************************************************************************/
typedef struct scvtxpair : public HeapElement
{
   // verweis auf die beiden punkte
   int   v1;
   int   v2;
   int   v3;
   
   // flag, ob es sich bei diesem Paar um eine 'outer Edge' handelt
   bool boundary;

   // ersatz - punkt
   vec3  v;
} cvtxpair;


class MeshSimplification
{
private:
	// Liste mit den Dreiecken
	face *faces;
	int  numFaces;
   // maximale anzahl flaechen pro punkt
   int  maxFacesPerPoint;
	// Liste der möglichen Kontraktionen
	cvtxpair *pairs;
	int      numPairs;
	// Liste der original - punkte
	vec3 *points;
	int  numPoints;
   // Liste der Fehlermatrizen
   QError *quadrics;
   // vertex structure mit mehr informationen
   cvtx *vertices;
   // Heap zum sortieren und einfügen
   MinHeap *minHeap;
   // marker vor invalid faces
   int_set invalidFaces;
   int_set invalidPairs;
   int_set connections;


public:
	/*****************************************************************************
	 * Erzeugt eine Object vom Typ MeshSimplification
	 ****************************************************************************/
	MeshSimplification(vec3 *points, int numPoints, face *faces, int numFaces);

   /*****************************************************************************
    * funktion um den Fehler anhand der quadrierten
    * Abstaende von den Flaechen zu ermitteln
    ****************************************************************************/
   void getQuadricError(cvtxpair &vp);

   /*****************************************************************************
    * überprüft eine struktur aus punkten und flächen auf ihre Gültigkeit
    ****************************************************************************/
   void checkStructure();

   /*****************************************************************************
    * Aufraeumen vor dem Freigeben des Objekts
    ****************************************************************************/
   ~MeshSimplification();

   /*****************************************************************************
    * methode zum reduzieren bis zu einer gewissen 
    * polygonanzahl
    *****************************************************************************/
   void decimate(int border, int *numCoord, vec3 **coords, int *numFaces, face **faces);
   
   /***************************************************************************
    * Hinzufügen eines Paares
    **************************************************************************/
   void addPair(halfedge *h, bool boundary);

   /***************************************************************************
    * Entfernt ein paar aus dem Heap und löscht die
    * referenzen in den punkten
    **************************************************************************/
   void removePair(int pairIndex);

   /**
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
    */
   void createPairs();

   /**
    * Indexlisten für Paare erzeugen
    *
    */
   void createPairIndices();

   /***************************************************************************
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
    **************************************************************************/
   int decimate();

   /***************************************************************************
    * löschen einer Vertex struktur,
    * die struktur wird gelöscht und in die liste der 
    * freien vertex - strukturen eingehängt
    **************************************************************************/
   void removeVertex(int index);
   
   /***************************************************************************
    * löschen einer Face struktur,
    * die struktur wird gelöscht und in die liste der 
    * freien face - strukturen eingehängt
    **************************************************************************/
   void removeFace(int faceIndex);
};


#endif

