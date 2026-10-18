//=============================================================================
//
// Grundgerüst für einen Edge-Collapse Algorithmus
//
// $Author: Andor $
// $Date: 11.12.03 23:00 $
// $Revision: 9 $
// 
//=============================================================================

#include "edgecollapse.h"
#include "../utility.h"

// begin namespace
namespace decimate {

// ----------------------------------------------------------------------------
// creates an edgecollapse object
// ----------------------------------------------------------------------------
CEdgeCollapse::CEdgeCollapse(int numPoints, const Vec3 *points, int numFaces, const Face *faces, IErrorSystem *errorSystem, IVertexPairCollection *vertexPairCollection)
{
   this->numPoints = numPoints;
   this->points    = new Vec3[numPoints];
   this->numFaces  = numFaces;
   this->faces     = new Face[numFaces];

   this->numValidPoints       = numPoints;
   this->errorSystem          = errorSystem;
   this->vertexPairCollection = vertexPairCollection;

   // copy data
   memcpy(this->points, points, sizeof(Vec3) * numPoints);
   memcpy(this->faces,   faces, sizeof(Face) * numFaces);

   // initialize subsystems
   errorSystem->Initialize(numPoints, this->points, numFaces, this->faces);
   vertexPairCollection->Initialize(numPoints, numFaces, this->faces);

   setup_timer();
   start_timer();

   // initialize edge collapse
   Initialize();
}


// ----------------------------------------------------------------------------
// deletes an edgecollapse object
// ----------------------------------------------------------------------------
CEdgeCollapse::~CEdgeCollapse()
{
   errorSystem->Destroy();
   vertexPairCollection->Destroy();

   delete[] points;
   delete[] faces;
   delete[] invalidPoints;
}


// ----------------------------------------------------------------------------
// Setup data structures
// ----------------------------------------------------------------------------
void CEdgeCollapse::Initialize()
{
   pairIndices   = new CSet<int>[numPoints];
   faceIndices   = new CSet<int>[numPoints];
   invalidPoints = new bool[numPoints];

   heap.Init(vertexPairCollection->size());

   for (IVertexPairCollection::iterator it = vertexPairCollection->begin(); it!=vertexPairCollection->end(); it++) {
      // get VertexPait
      VertexPair *vp = it;
      // compute error for possible contraction and
      // add it to priority queque
      errorSystem->Optimize(vp->p1, vp->p2,vp->boundary, vp->heapValue, vp->p);
      heap.Push(it);

      // setup point to pair index list
      int index = vertexPairCollection->find(it);
      //
      pairIndices[it->p1].Insert(index);
      pairIndices[it->p2].Insert(index);
   }

   for (int i=0; i<numFaces; i++) {
      // setup point to face reference
      faceIndices[faces[i].a].Insert(i);
      faceIndices[faces[i].b].Insert(i);
      faceIndices[faces[i].c].Insert(i);
   }

   // validate structure
   for (int i=0; i<numPoints; i++) {
      invalidPoints[i] = false;
      if (faceIndices[i].Size() == 0) MarkInvalid(i);
   }
}


// ----------------------------------------------------------------------------
// decimates the model
// ----------------------------------------------------------------------------
void CEdgeCollapse::DecimateToPoints(int numDesiredPoints, int *numCoord, Vec3 **coords, int *numFaces, Face **faces)
{
   // get number of points to do
   double difference = numValidPoints - numDesiredPoints;
   int done       = 100;

   // Decimation loop
   while ((numDesiredPoints < numValidPoints) && Decimate()) {
      // show some progress to user
      if ((((numValidPoints - numDesiredPoints) * 100.0 )/ difference) < done) {
         show_progress(); printf("% d           ", numValidPoints);
         done = (int)(((numValidPoints - numDesiredPoints) * 100.0 )/ difference);
      }
   }
   printf("\r Time elapsed %lf\n", stop_timer());

   // correct point indices
   Compress();

   // copy coordinate information
   *numCoord = numValidPoints;
   *coords = new Vec3[numValidPoints];
   memcpy(*coords, points, sizeof(Vec3) * numValidPoints);

   // copy face structure
   *numFaces = this->numFaces;
   *faces = new Face[this->numFaces];
   memcpy(*faces, this->faces, sizeof(Face) * this->numFaces);
}


// ----------------------------------------------------------------------------
// remap the point indices within a face
// ----------------------------------------------------------------------------
void CEdgeCollapse::RemapFaceIndices(CSet<int> &indices, int new_value, int old_value)
{
   for (CSet<int>::iterator it=indices.begin(); it!=indices.end(); it++) {
      // change point indices
      Face &fref = faces[*it];

      if (fref.a == old_value) fref.a = new_value;
      if (fref.b == old_value) fref.b = new_value;
      if (fref.c == old_value) fref.c = new_value;
   }
}


// ----------------------------------------------------------------------------
// compress data structure
// ----------------------------------------------------------------------------
void CEdgeCollapse::Compress()
{
   int counter = 0;

   for (int i=0; i<numPoints; i++) {
      // serach for next valid point
      if (invalidPoints[i])
         continue;
      // remap point
      if (counter!=i) {
         points[counter] = points[i];
         RemapFaceIndices(faceIndices[i], counter, i);
      }
      counter++;
   }
}


// ----------------------------------------------------------------------------
// decimates the model
// ----------------------------------------------------------------------------
bool CEdgeCollapse::Decimate()
{
   VertexPair *vertexPair = (VertexPair *) heap.Pop();

   if (vertexPair) {
      // get index of pair and remove it
      int index = vertexPairCollection->find(vertexPair);

      pairIndices[vertexPair->p1].Erase(index);
      pairIndices[vertexPair->p2].Erase(index);

      // Eintragen der neuen Koordinaten
      points[vertexPair->p1] = vertexPair->p;

      // Fehlersystem benachrichtigen
      errorSystem->Update(vertexPair->p1, vertexPair->p2);

      // p2 als ungültig markieren
      MarkInvalid(vertexPair->p2);

      // in allen Dreiecken p2-index durch p1-index ersetzen
      UpdateFaceIndices(vertexPair->p1, vertexPair->p2);

      // in allen Paaren p2-index durch p1-index ersetzen
      UpdatePairIndices(vertexPair->p1, vertexPair->p2);

      RemoveInvalidFaces();
      RemoveInvalidPairs();

      // Fehler für alle an p1 angrenzenden Paare neu berechnen
      UpdateError(vertexPair->p1);

      return true;
   }
   return false;
}


// ----------------------------------------------------------------------------
// marks a point as invalid and reduces 
// number of valid points
// ----------------------------------------------------------------------------
void CEdgeCollapse::MarkInvalid(int index)
{
   if (!invalidPoints[index]) {
      invalidPoints[index] = true;
      numValidPoints--;
   }
#ifdef DEBUG
   else
      printf("CEdgeCollapse::MarkInvalid for already marked point\n");
#endif
}


// ----------------------------------------------------------------------------
// removes a point wich is not invalid yet
// ----------------------------------------------------------------------------
void CEdgeCollapse::RemovePoint(int index)
{
   // Assert invalid state
   if (invalidPoints[index])
      return;

   // remove point and pairs
   invalidPairs.Merge(pairIndices[index]);
   pairIndices[index].Clear();
   faceIndices[index].Clear();
   MarkInvalid(index);
}


// ----------------------------------------------------------------------------
// Entfernt eine ungültige Fläche, die Liste der Flächen
// wird gepackt - d.h. die ungültige Fläche innerhalb der Liste wird
// durch die letzte gültige Flächen ersetzt, die Indexlisten werden
// entsprechend korrigiert
// ----------------------------------------------------------------------------
void CEdgeCollapse::RemoveInvalidFace(int index)
{
   // FaceIndex von allen vertex strukturen entfernen
   Face *fp = &faces[index];
   //
   faceIndices[fp->a].Erase(index);
   faceIndices[fp->b].Erase(index);
   faceIndices[fp->c].Erase(index);

   // punkte auf gültigkeit überprüfen
   if (faceIndices[fp->a].Size() == 0) RemovePoint(fp->a);
   if (faceIndices[fp->b].Size() == 0) RemovePoint(fp->b);
   if (faceIndices[fp->c].Size() == 0) RemovePoint(fp->c);

   numFaces--;

   // letzte face-struktur an die freigewordene stelle kopieren
   if (index < numFaces) {
      // copy face info to free place
      faces[index] = faces[numFaces];
      fp           = &faces[index];

      // index in facelisten korrigieren
      faceIndices[fp->a].Swap(numFaces, index);
      faceIndices[fp->b].Swap(numFaces, index);
      faceIndices[fp->c].Swap(numFaces, index);
/*
      // index in facelisten korrigieren
      faceIndices[fp->a].Erase(numFaces);
      faceIndices[fp->b].Erase(numFaces);
      faceIndices[fp->c].Erase(numFaces);
      //
      faceIndices[fp->a].Insert(index);
      faceIndices[fp->b].Insert(index);
      faceIndices[fp->c].Insert(index);*/
   }
}


// ----------------------------------------------------------------------------
// Entfernt die ungültigen Flächen, die Liste der Flächen
// wird gepackt - d.h. ungültige Flächen innerhalb der Liste werden
// durch gültige Flächen am Ende der Liste ersetzt
// ----------------------------------------------------------------------------
void CEdgeCollapse::RemoveInvalidFaces()
{
   for (CSet<int>::iterator rit=invalidFaces.rbegin(); rit!=invalidFaces.rend(); rit--) {
      RemoveInvalidFace(*rit);
   }
   invalidFaces.Clear();
}


// ----------------------------------------------------------------------------
// ermittelt den Fehler für alle an durch index
// bestimmten Punkt neu
// ----------------------------------------------------------------------------
void CEdgeCollapse::UpdateError(int index)
{
   CSet<int> &indices = pairIndices[index];

   for (CSet<int>::iterator it=indices.begin(); it!=indices.end(); it++) {
      // Fehler für paar neu bestimmen
      VertexPair *vp = (*vertexPairCollection)[*it];
      errorSystem->Optimize(vp->p1, vp->p2, vp->boundary, vp->heapValue, vp->p);
      //
      heap.Repair(vp);
   }
}


// ----------------------------------------------------------------------------
// removes invalid pairs from heap and index lists
// ----------------------------------------------------------------------------
void CEdgeCollapse::RemoveInvalidPairs()
{
   for (CSet<int>::iterator it=invalidPairs.begin(); it!=invalidPairs.end(); it++) {
      VertexPair *vertexPair = (*vertexPairCollection)[*it];
      //
      pairIndices[vertexPair->p1].Erase(*it);
      pairIndices[vertexPair->p2].Erase(*it);
      //
      heap.Remove(vertexPair);
   }
   invalidPairs.Clear();
}


// ----------------------------------------------------------------------------
// ersetzt in allen von faceIndices angegebenen Flächen
// den index von v2 durch den Index von v1
// ----------------------------------------------------------------------------
void CEdgeCollapse::UpdateFaceIndices(int p1, int p2)
{
   CSet<int> &indices = faceIndices[p2];
   for (CSet<int>::iterator it=indices.begin(); it!=indices.end(); it++) {
      UpdateFaceIndices(*it, p1, p2);
   }
   faceIndices[p1].Merge(faceIndices[p2]);
   faceIndices[p2].Clear();
}


// ----------------------------------------------------------------------------
// ersetzt in durch index bestimmten Fläche
// den index von p2 durch den Index von p1
// ----------------------------------------------------------------------------
void CEdgeCollapse::UpdateFaceIndices(int index, int p1, int p2)
{
   Face &f = faces[index];

   if ((f.a==p1 || f.b==p1 || f.c==p1) && (f.a==p2 || f.b==p2 || f.c==p2)) {
      // set to invalid state
      invalidFaces.Insert(index);
   }
   // change index
   if (f.a == p2) f.a = p1;
   if (f.b == p2) f.b = p1;
   if (f.c == p2) f.c = p1;
}


// ----------------------------------------------------------------------------
// ersetzt in allen von faceIndices angegebenen Paaren
// den index von p2 durch den Index von p1
// ----------------------------------------------------------------------------
void CEdgeCollapse::UpdatePairIndices(int p1, int p2)
{
   for (CSet<int>::iterator it=pairIndices[p2].begin(); it!=pairIndices[p2].end(); it++) {
      UpdatePairIndices(*it, p1, p2);
   }
   pairIndices[p1].Merge(pairIndices[p2]);
   pairIndices[p2].Clear();
}


// ----------------------------------------------------------------------------
// ersetzt in dem durch index bestimmten Paar
// den index von p2 durch den Index von p1
// ----------------------------------------------------------------------------
void CEdgeCollapse::UpdatePairIndices(int index, int p1, int p2)
{
   VertexPair *p = (*vertexPairCollection)[index];

   if ((p->p1==p1 || p->p2==p1) && (p->p1==p2 || p->p2==p2)) {
      // set to invalid state
      invalidPairs.Insert(index);
   }
   // change index
   if (p->p1 == p2) p->p1 = p1;
   if (p->p2 == p2) p->p2 = p1;

   // check for idents, this is possible since after a
   // collapse operation two identical pairs can exist
   if (p->p1 == p->p2) {
      invalidPairs.Insert(index);
   }
}

} // end namespace
