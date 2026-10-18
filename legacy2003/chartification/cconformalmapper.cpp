//=============================================================================
//
// CConformalMapper : Class Implementation
//
// Abbildung von Charts in das TexturkoordinatenSystem
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 13 $
// 
//=============================================================================
#include "cconformalmapper.h"

// c++ includes
#include <map>

namespace chartification {

// ----------------------------------------------------------------------------
// Constructs a CConformalMapper class
// ----------------------------------------------------------------------------
CConformalMapper::CConformalMapper(Vec3 *points)
{
   this->points = points;
}

// ----------------------------------------------------------------------------
// Deletes a CConformalMapper
// ----------------------------------------------------------------------------
CConformalMapper::~CConformalMapper() {}

// ----------------------------------------------------------------------------
// copies a result set
// ----------------------------------------------------------------------------
void CopyResult(Vertex *vertices, double *vec, int offset, int num)
{
   for (int i=0; i<num; i++)
   {
      vertices[offset + i].u = (float) vec[i];
      vertices[offset + i].v = (float) vec[num + i];
   }
}

// ----------------------------------------------------------------------------
// calculates chart geometric diameter
// maximum geometric distance between two poits
// ----------------------------------------------------------------------------
float CConformalMapper::Diameter(int *p1, int *p2, CChart *chart)
{
   // get parameter from chart
   Vertex *vertices = chart->Vertices();
   Edge   *edges    = chart->Boundaries();
   int    size      = chart->BoundaryCount();

   // initialize maximum to zero
   float max = 0.0f;
   // initialize p1, p2
   *p1 = 0;
   *p1 = 1;

   for (int i=0; i<size; i++) {
      // get startpoint (map vertex to point id)
      Vec3 &start = points[vertices[edges[i].a].index];
      for (int j=i + 1; j<size; j++) {
         // get endpoint (map vertex to point id)
         Vec3 &end   = points[vertices[edges[j].b].index];

         // calculate distance
         float dist = (end - start).length();

         // check for new maximum
         if (dist > max) {
            // set new maximum
            *p1 = i;//edges[i].a;
            *p2 = j;//edges[j].b;
            max = dist;
         }
      }
   }
   // return maximum distance
   return max;
}

// ----------------------------------------------------------------------------
// returns the triangle gradient
//
// optimization //TODO: eliminate known zero multiplies
// ----------------------------------------------------------------------------
void CConformalMapper::Gradient(
   // outgoing parameters
   Complex *w1, Complex *w2, Complex *w3, float *dT, 
   // incoming parameters
   Face *f, Vertex *vertices
   )
{
   // local references and vectors
   Vec3 &a = points[vertices[f->a].index];
   Vec3 &b = points[vertices[f->b].index];
   Vec3 &c = points[vertices[f->c].index];

   Vec3 ba = b - a;
   Vec3 ca = c - a;

   // Use Gram-Schmidt process to
   // define an orthonormal basis
   Vec3 e1 = ba;
   Vec3 e2 = ca;
   orthonormal(e1, e2);

   // get triangle coordinates in
   // orthonormal 2D basis
   float x1 = 0.0f;
   float y1 = 0.0f;
   float x2 = (ba) * e1;
   float y2 = (ba) * e2;
   float x3 = (ca) * e1;
   float y3 = (ca) * e2;

   // setup gradients
   w1->r = x3-x2;  w1->i = y3-y2;
   w2->r = x1-x3;  w2->i = y1-y3;
   w3->r = x2-x1;  w3->i = y2-y1;

   // get area dT
   *dT = (x1*y2 - y1*x2) + (x2*y3 - y2*x3) + (x3*y1 - y3*x1);
}

// ----------------------------------------------------------------------------
// chooses the pinned coordinates and put the vertex info at the end
// of the data structure
//
// @return number of pinned coordinates
//
// ----------------------------------------------------------------------------
int CConformalMapper::InitPinned(CChart *chart)
{
   // variables for index of pinned vertices
   int p1, p2;
   // calculate p1, p2 by geometric diameter
   float dist = Diameter(&p1, &p2, chart);

   // put both vertices at the end
   Vertex *vertices = chart->Vertices();
   int    size      = chart->VertexCount();

   // put at end of list
   Vertex ve1 = vertices[size - 2];
   Vertex ve2 = vertices[size - 1];
   //
   Vertex vp1 = vertices[p1];
   Vertex vp2 = vertices[p2];
   // swap
   vertices[size - 2] = vp1;
   vertices[size - 1] = vp2;
   //
   vertices[p1]       = ve1;
   vertices[p2]       = ve2;

   // bootstrap uv coordinates
   vertices[size - 2].u = 0.0f;
   vertices[size - 2].v = 0.0f;
   vertices[size - 1].u = 1.0f;
   vertices[size - 1].v = 1.0f;

   // build a relocation map
   std::map<int, int> relocation;
   //
   relocation[p1]       = size - 2;
   relocation[p2]       = size - 1;
   relocation[size - 2] = p1;
   relocation[size - 1] = p2;

   // relocate face indices
   Face *faces = chart->Faces();
   for (int i=0; i<chart->FaceCount(); i++) {
      // get current face
      Face &face = faces[i];
      // relocate indices
      if (relocation.find(face.a) != relocation.end()) face.a = relocation[face.a];
      if (relocation.find(face.b) != relocation.end()) face.b = relocation[face.b];
      if (relocation.find(face.c) != relocation.end()) face.c = relocation[face.c];
   }

   // relocate edge indices
   Edge *edges = chart->Boundaries();
   for (int i=0; i<chart->BoundaryCount(); i++) {
      // get current face
      Edge &edge = edges[i];
      // relocate indices
      if (relocation.find(edge.a) != relocation.end()) edge.a = relocation[edge.a];
      if (relocation.find(edge.b) != relocation.end()) edge.b = relocation[edge.b];
   }
   return 2;
}

// ----------------------------------------------------------------------------
// inserts a complex value into a matrix
//
// complex matrix is decomposed as follows
//
// AC = MC (M Complex)
//
// A = | Mr -Mi |  Mr stands of real part of MC
//     | Mi  Mr |  Mi is imaginary part of MC
//
// @param smat<float>, pointer to a sparse matrix
// @param i          , row to insert
// @param j          , column to insert
// @param rows2      , number of rows / 2
// @param columns2   , number of columns / 2
//
// ----------------------------------------------------------------------------
void CConformalMapper::Insert(SMCR *mat, Complex *c, int i, int j, int rows2, int columns2)
{
   // validate indices
   if ((i>=rows2) || (j>=columns2) || (i<0) || (j<0))
      return;
   // insert upper left
   mat->insert(i, j, c->r);
   // insert upper right
   mat->insert(i, j + columns2, -c->i);
   // insert bottom left
   mat->insert(i + rows2, j, c->i);
   // insert bottom right
   mat->insert(i + rows2, j + columns2, c->r);
}

// ----------------------------------------------------------------------------
// Maps a whole chartset
// ----------------------------------------------------------------------------
void CConformalMapper::Parametrize(CChartSet *charts)
{
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // Map current chart
      Parametrize(*it);
   }
}

// ----------------------------------------------------------------------------
// chart parametrization using least squares
//
// set up matrices f - free
//                 p - pinned
//                 1 - real
//                 2 - imaginary
// ----------------------------------------------------------------------------
void CConformalMapper::Parametrize(CChart *chart)
{
   printf ("Least Squares Conformal Maps\n");

   // Definition der lokalen variablen
   int      numPoints, numFaces, numFree, numPinned, iter, bsize;
   double  *b, *x, *pinned;
   double   err;
   Vertex  *vertices;

   // Initialize with pinned vertices
   numPinned = InitPinned(chart);

   // get parameters from chart
   numPoints = chart->VertexCount();
   numFaces  = chart->FaceCount();
   numFree   = numPoints - numPinned;
   vertices  = chart->Vertices();

   // create our sparse matrices
   SMCR am(numFaces * 2, numFree * 2);
   SMCR bm(numFaces * 2, numPinned * 2);

   // Matrizen initialisieren
   InitSMCR(am, bm, numPinned, chart);


   // Vektor mit den gepinnten Koordinaten aufbauen
   pinned = new double[numPinned * 2];
   for (int i=0; i<numPinned; i++) {
      pinned[            i] = vertices[numFree + i].u;
      pinned[numPinned + i] = vertices[numFree + i].v;
   }

   // Vektor b erzeugen
   bsize = bm.rows;
   b = new double[bsize];

   // b = - bm * pinned
   bm.vecmul(b, pinned, bsize);
   for (int i=0; i<bsize; i++) b[i] = -b[i];

   // matrix am transponieren
   am.transpose();
   
   // b = am^T * b
   am.vecmul(bsize, b);
   
   printf ("matrix mul %d ", am.rows);
   setup_timer();
   start_timer();
   
   // am = A^T * A
   am.transform_atas();

   printf(": %lf milliseconds\n", stop_timer());
   printf ("matrix solve %d ", am.rows);
   start_timer();

   x = am.cgiters(b, 1e-12, 15000, &iter, &err);

  	// print elapsed time in milliseconds
   printf(": %lf milliseconds\n", stop_timer());
   printf("iter=%4d err=%14.12f\n\n", iter, err);

   // set texture coordinates
   CopyResult(vertices, x, 0, numFree);

   // scale texture to be equal to chart dimensions
   //ScaleTexture(chart);

   // cleanup vectors
   delete[] x;
   delete[] b;
   delete[] pinned;
}

// ----------------------------------------------------------------------------
// Matrizen a und B initialisieren
// ----------------------------------------------------------------------------
void CConformalMapper::InitSMCR(SMCR &am, SMCR &bm, int numPinned, CChart *chart)
{
   // lokale Variablen
   int      numPoints, numFaces, numFree;
   Complex  w1, w2, w3;
   float    dT, sqrtdT;
   Vertex  *vertices;
   Face    *faces;

   // Werte des Charts laden
   numPoints = chart->VertexCount();
   numFaces  = chart->FaceCount();
   numFree   = numPoints - numPinned;
   vertices  = chart->Vertices();
   faces     = chart->Faces();


   // iterate over all faces to fill matrices
   for (int i=0; i<numFaces; i++) {
      // get triangle
      Face *face = &faces[i];

      // load gradients
      Gradient(&w1, &w2, &w3, &dT, face, vertices);

      // calculate what we need
      sqrtdT = sqrtf(dT);

      // divide gradients by sqrt of area
      w1 /= sqrtdT;
      w2 /= sqrtdT;
      w3 /= sqrtdT;

      // fill martix A
      Insert(&am, &w1, i, face->a, numFaces, numFree);
      Insert(&am, &w2, i, face->b, numFaces, numFree);
      Insert(&am, &w3, i, face->c, numFaces, numFree);

      // fill matrix B
      Insert(&bm, &w1, i, face->a - numFree, numFaces, numPinned);
      Insert(&bm, &w2, i, face->b - numFree, numFaces, numPinned);
      Insert(&bm, &w3, i, face->c - numFree, numFaces, numPinned);
   }
}

} // end namespace
