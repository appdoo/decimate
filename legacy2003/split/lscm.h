/******************************************************************************
 *
 * Implementation of the least squares conformal maps algorithm
 *
 * [1] 
 *     Bruno Lévy Sylvain Petitjean Nicolas Ray Jérome Maillot
 *
 *
 *****************************************************************************/
#ifndef LSCM_H
#define LSCM_H

// Local Includes
#include "../math/math3D.h"
#include "../math/halfedge.h"
#include "../math/sparsematrix.h"
#include "../math/smcr.h"
#include "../math/algorithm.h"
#include "..\utility.h"


class LSCM
{
public:
   // faces belonging to this LSCM object
   int numFaces;
   Face *faces;

   // array for pinned vertices
   int   numPinned;
   int   *pinned;
   float *upinned;
   float *vpinned;

   // points
   int numPoints;
   Vec3 *points;

   // path
   int numPath;
   int *path;

   // typedef for relocation table
   typedef std::map<int, int> mapping;


   /***************************************************************************
    * constructs an empty LSCM object
    **************************************************************************/
   LSCM()
   {
      // faces belonging to this LSCM object
      numFaces  = 0;
      faces     = NULL;

      // pinned vertices
      numPinned = 0;
      pinned    = NULL;
      upinned   = NULL;
      vpinned   = NULL;

      // points
      numPoints = 0;
      points    = NULL;

      // path
      path      = NULL;
      numPath   = 0;
   }

   /***************************************************************************
    * cleanup and exit
    **************************************************************************/
   ~LSCM()
   {
      // cleanup
      if (faces)   delete[] faces;
      if (pinned)  delete[] pinned;
      if (upinned) delete[] upinned;
      if (vpinned) delete[] vpinned;
      if (points)  delete[] points;
      if (path)    delete[] path;
   }

   /***************************************************************************
    * initializes the LSCM object with the needed data structures
    **************************************************************************/
   void init(int_set *face_set, halfedge_set *boundaries, Vec3 *points, Face *faces, int pinned1, int pinned2, int seed)
   {
      // build face structure
      init_faces(face_set, faces);
      // build the boundary order
      init_pinned(pinned1, pinned2, seed, points, faces, boundaries);
      // convert path
      init_path(boundaries);

      // build a temporary relocation table
      mapping *relocation = build_relocation();
      // copy point structure
      copy_points(relocation, points);
      // relocate face - indices
      relocate_faces(relocation);
      // relocate boundary - point indices
      relocate_pinned(relocation);
      // relocate path
      relocate_path(relocation);
      // cleanup
      delete relocation;
   }

   /***************************************************************************
    * builds a relocation table for point indices
    **************************************************************************/
   mapping *build_relocation()
   {
      // create a relocation table using
      // pinned vertices and faces
      mapping *relocation = new mapping();
      int index = 0;
      // enter pinned vertex indices
      for (int i=0; i<numPinned; i++) {
         (*relocation)[pinned[i]] = index++;
      }
      // enter face indices
      for (int i=0; i<numFaces; i++) {
         // get face
         Face *f = &faces[i];
         // create relocated point - pairs
         if (relocation->find(f->a) == relocation->end()) (*relocation)[f->a] = index++;
         if (relocation->find(f->b) == relocation->end()) (*relocation)[f->b] = index++;
         if (relocation->find(f->c) == relocation->end()) (*relocation)[f->c] = index++;
      }
      // return relocation table
      return relocation;
   }

   /***************************************************************************
    * relocate function
    **************************************************************************/
   int relocate(mapping *relocation, int index)
   {
      return numPoints - (*relocation)[index] - 1;
   }
   
   /***************************************************************************
    * relocate face - point indices by relocation table
    **************************************************************************/
   void relocate_faces(mapping *relocation)
   {
      for (int i=0; i<numFaces; i++) {
         faces[i].a = relocate(relocation, faces[i].a);
         faces[i].b = relocate(relocation, faces[i].b);
         faces[i].c = relocate(relocation, faces[i].c);
      }
   }

   /***************************************************************************
    * relocate boundary - point indices by relocation table
    **************************************************************************/
   void relocate_pinned(mapping *relocation)
   {
      for (int i=0; i<numPinned; i++) {
         pinned[i] = relocate(relocation, pinned[i]);;
      }
   }

   /***************************************************************************
    * relocate path indices by relocation table
    **************************************************************************/
   void relocate_path(mapping *relocation)
   {
      for (int i=0; i<numPath; i++) {
         path[i] = relocate(relocation, path[i]);;
      }
   }

   /***************************************************************************
    * copies points into LSCM object using the relocation table
    * copy is done in reverse order, to bring the pinned vertices
    * to the end of the vertex array
    **************************************************************************/
   void copy_points(mapping *relocation, Vec3 *srcpoints)
   {
      // create new array of points
      numPoints = relocation->size();
      points = new Vec3[numPoints];
      // copy points, putting pinned vertices to the end
      for (mapping::iterator it=relocation->begin(); it!=relocation->end(); it++) {
         points[numPoints - (it->second) - 1] = srcpoints[it->first];
      }
   }

   /***************************************************************************
    * initializes the face set of this LSCM object
    **************************************************************************/
   void init_faces(int_set *face_set, Face *faces)
   {
      // create private face structure
      this->faces = new Face[face_set->size()];
      numFaces = 0;
      // copy face structure
      for (int_iterator it=face_set->begin(); it!=face_set->end(); it++)
         this->faces[numFaces++] = faces[*it];
   }

   /***************************************************************************
    * keep track of chart boundaries
    **************************************************************************/
   void init_path(halfedge_set *boundaries)
   {
      // allocate path
      numPath = boundaries->size() * 2;
      path    = new int[numPath];
      // local variable
      int index = 0;

      for (halfedge_set::iterator it=boundaries->begin(); it!=boundaries->end(); it++) {
         path[index++] = (*it)->next->next->end;
         path[index++] = (*it)->end;
      }
   }

   /***************************************************************************
    * initializes the set of pinned vertices for this
    * lscm - object
    *
    * //TODO:
    * as opposed in [1] we use only those two vertices, 
    * maximizing the distance between them
    *
    * for performance reasons take one point from the seed face
    * (the first face) and the point from chart boundary with maximum distance
    **************************************************************************/
   void init_pinned(int pinned1, int pinned2, int seed, Vec3 *points, Face *faces, halfedge_set *boundaries)
   {
      numPinned = 2;
      pinned  = new int[numPinned];
      upinned = new float[numPinned];
      vpinned = new float[numPinned];

      pinned[0] = pinned1;
      pinned[1] = pinned2;

      upinned[0] = 0.0f;
      vpinned[0] = 0.0f;
      upinned[1] = 1.0f;
      vpinned[1] = 1.0f;
   }

   /***************************************************************************
    * returns the triangle gradient value as a complex value
    *
    * optimization //TODO: eliminate known zero multiplies
    *
    **************************************************************************/
   void gradient_lscm(
      // outgoing parameters
      Complex *w1, Complex *w2, Complex *w3, float *dT, 
      // incoming parameters
      Face *f, Vec3 *points
      )
   {
      // local references and vectors
      Vec3 &a = points[f->a];
      Vec3 &b = points[f->b];
      Vec3 &c = points[f->c];

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

   /***************************************************************************
    * inserts a complex value into a matrix
    *
    * complex matrix is decomposed as follows
    *
    * AC = MC (M Complex)
    *
    * A = | Mr -Mi |  Mr stands of real part of MC
    *     | Mi  Mr |  Mi is imaginary part of MC
    *
    * @param SMat<float>, pointer to a sparse matrix
    * @param i          , row to insert
    * @param j          , column to insert
    * @param rows2      , number of rows / 2
    * @param columns2   , number of columns / 2
    *
    **************************************************************************/
   void insert(SMat *mat, Complex *c, int i, int j, int rows2, int columns2)
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

   /***************************************************************************
    * chart parametrization using least squares
    **************************************************************************/
   void parametrize(int *dstNumPoints, int *dstNumFaces, Vec3 **dstPoints, Face **dstFaces, Tex2 **dstTexture)
   {
      // set up matrices f - free
      //                 p - pinned
      //                 1 - real
      //                 2 - imaginary
      int numFree = numPoints - numPinned;
      // create our sparse matrices
      SMat am(numFaces * 2, numFree * 2);
      SMat bm(numFaces * 2, numPinned * 2);

      // placeholder for gradient generation
      Complex w1, w2, w3;
      float dT;
      // iterate over all faces to fill matrices
      for (int i=0; i<numFaces; i++) {
         // get triangle
         Face *f = &faces[i];
         // load gradients
         gradient_lscm(&w1, &w2, &w3, &dT, f, points);
         // calculate what we need
         float sqrtdT = sqrtf(dT);
         // divide gradients by sqrt of area
         w1 /= sqrtdT;
         w2 /= sqrtdT;
         w3 /= sqrtdT;
         // fill A martix
         insert(&am, &w1, i, f->a, numFaces, numFree);
         insert(&am, &w2, i, f->b, numFaces, numFree);
         insert(&am, &w3, i, f->c, numFaces, numFree);
         // fill B matrix
         insert(&bm, &w1, i, f->a - numFree, numFaces, numPinned);
         insert(&bm, &w2, i, f->b - numFree, numFaces, numPinned);
         insert(&bm, &w3, i, f->c - numFree, numFaces, numPinned);
      }

      // create pinned vertices vector
      SMat b(numPinned * 2, 1);
      for (int i=0; i<numPinned; i++) {
         b.insert(i, 0, upinned[i]);
         b.insert(numPinned + i, 0, vpinned[i]);
      }
      SMat b2 = bm * b;
      b2.neg();
      
      double *bf = new double[b2.m];
      for (int i=0; i<b2.m; i++) {
         bf[i] = b2.rows[i] ? b2.rows[i]->value : 0.0;
      }
      am.transpose();
      SMCR sma(am);
      sma.vecmul(b2.m, bf);

      printf ("matrix mul %d ", sma.rows);

      setup_timer();
      start_timer();
      sma.transform_atas();
	   double time = stop_timer();
      printf(": %lf milliseconds\n", time);

      printf ("matrix solve %d ", sma.rows);

      start_timer();
      int iter;
      double err;
      double *x = sma.cgiters(bf, 1e-8, 15000, &iter, &err);
      //double *x = smata.bcgiter(bf, 1, 1e-12, 15000, &iter, &err);
	   time = stop_timer();

  	   // print elapsed time in milliseconds
      printf(": %lf milliseconds\n", time);
      printf("iter=%4d err=%14.12f\n\n", iter, err);


      // copy outgoing stuff
      *dstNumPoints = this->numPoints;
      *dstNumFaces  = this->numFaces;
      //
      *dstPoints  = new Vec3[this->numPoints];
      *dstTexture = new Tex2[this->numPoints];
      *dstFaces   = new Face[this->numFaces];
      //
      memcpy(*dstPoints, this->points, sizeof(Vec3) * this->numPoints);
      memcpy(*dstFaces, this->faces, sizeof(Face) * this->numFaces);
      //
      for (int i=0; i<numFree; i++)
      {
         (*dstTexture)[i].u = (float) x[i];
         (*dstTexture)[i].v = (float) x[numFree + i];
      }
      for (int i=0; i<numPinned; i++)
      {
         (*dstTexture)[numFree + i].u = (float) b.getvalue(i, 0);
         (*dstTexture)[numFree + i].v = (float) b.getvalue(numPinned + i, 0);
      }

      // cleanup
      delete[] x;
      delete[] bf;
   }

   /***************************************************************************
    * returns the relocated point paths
    **************************************************************************/
   void getPath(int *numPath, int **path)
   {
      *numPath = this->numPath;
      *path    = new int[this->numPath];

      memcpy(*path, this->path, sizeof(int) * this->numPath);
   }
};

#endif