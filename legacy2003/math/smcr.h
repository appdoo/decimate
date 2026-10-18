// ============================================================================
//
// Class SMCR : Sparse Matrix, compressed row
//
// Definition and Implementation
//
// $Author: Andor $
// $Date: 12.12.03 2:27 $
// $Revision: 7 $
//
// ============================================================================

#include "sparsematrix.h"

#ifndef SMCR_H
#define SMCR_H

class SMCR
{
public:
   // Anzahl Zeilen und Spalten
   int rows;
   int cols;

   // Zeiger auf werte und indices
   real **values;
   int **indices;

   // Anzahl Elemente pro Zeile
   int *sizes;

   // temporäre Arrays für einfügeoperationen
   real *tmpvalues;
   int *tmpindices;

   // global variable für bcgiter lösungsverfahren
   real *sa;

   // -------------------------------------------------------------------------
   // constructs a Sparse row compressed matrix
   // -------------------------------------------------------------------------
   SMCR(int rows, int cols)
   {
      init (rows, cols);
   }

   // -------------------------------------------------------------------------
   // constructs a Sparse row compressed matrix
   // -------------------------------------------------------------------------
   SMCR(SMat &mat)
   {
      init (mat.m, mat.n);
      // 
      for (int i=0; i<mat.m; i++) {
         if (mat.rows[i]) insert_row(mat.rows[i]);
      }
   }

   // -------------------------------------------------------------------------
   // constructs a Sparse row compressed matrix
   // -------------------------------------------------------------------------
   void init(int rows, int cols)
   {
      this->rows = rows;
      this->cols = cols;
      //
      // allocate spaces
      values     = new real *[rows];
      indices    = new int *[rows];
      sizes      = new int[rows];
      tmpvalues  = new real[MAX(cols * 2, rows * 2)];
      tmpindices = new int[MAX(cols * 2, rows * 2)];
      // clear 
      memset(values, 0, sizeof(real *) * rows);
      memset(indices, 0, sizeof(int *) * rows);
      memset(sizes, 0, sizeof(int) * rows);
   }

   // -------------------------------------------------------------------------
   // deletes a Sparse row compressed matrix
   // -------------------------------------------------------------------------
   ~SMCR()
   {
      // free compressed row data
      for (int i=0; i<rows; i++) {
         if (values[i]) delete[] values[i];
         if (indices[i]) delete[] indices[i];
      }
      // free memory
      delete[] values;
      delete[] indices;
      delete[] sizes;
      delete[] tmpvalues;
      delete[] tmpindices;
   }

   // -------------------------------------------------------------------------
   // inserts a Sparse row
   // -------------------------------------------------------------------------
   void insert_row(Sparse *row)
   {
      int tmpsize = 0;
      Sparse *loop = row;
      // copy row data
      while (loop) {
         tmpvalues[tmpsize] = loop->value;
         tmpindices[tmpsize] = loop->col;
         tmpsize++;
         //
         loop = loop->next_col;
      }
      insert_tmp(row->row, tmpsize);
   }

   // -------------------------------------------------------------------------
   // search for element, return index or -1
   // -------------------------------------------------------------------------
   int ContainsCol(int *indices, int size, int &x)
   {
      // konsistenzprüfung
      if (size < 1) return -1;

      // binäre suche
      int anf = 0;
      int end = (size - 1);
      int m   = (size - 1) >> 1;

      while (anf <= end && x!=indices[m]) {
         if (x<indices[m]) end = m - 1;
         else              anf = m + 1;

         m = (anf + end) >> 1;
      }
      return (x == indices[m]) ? m : -1;
   }

   // -------------------------------------------------------------------------
   // Fügt ein Element in das Set ein
   // -------------------------------------------------------------------------
   int Insert(int row, int col, real &x)
   {
      // lokale variablen laden
      int   size     = sizes[row];
      real *rvalues  = values[row];
      int  *rindices = indices[row];

      // check for set
      int pos = ContainsCol(rindices, size, col);
      // set does only contain an object 1 time
      if (pos != -1) return pos;

      int  *tmpidx = new int[size + 1];
      real *tmpval = new real[size + 1];

      // sortieren durch direktes einfügen
      // bestehendes set ist bereits sortiert
      int ug = 0; 
      int og = size - 1;
      int j  = 0;
      // bin suche
      while (ug<=og) {
         j = (ug + og) >> 1;
         if (col<rindices[j]) og = j - 1;
         else                 ug = j + 1;
      }

      // oberen teil verschieben
      for (j = size; j>ug; j--) {
         tmpidx[j] = rindices[j - 1];
         tmpval[j] = rvalues[j - 1];
      }

      // element einfügen
      tmpidx[ug] = col;
      tmpval[ug] = x;

      // unteren teil kopieren
      for (j = 0; j<ug; j++) {
         tmpidx[j] = rindices[j];
         tmpval[j] = rvalues[j];
      }
      size++;

      // alte arrays löschen
      if (values[row])  delete[] values[row];
      if (indices[row]) delete[] indices[row];
      
      sizes[row]   = size;
      values[row]  = tmpval;
      indices[row] = tmpidx;

      // position zurückliefern
      return ug;
   }

   // -------------------------------------------------------------------------
   // inserts an element
   // -------------------------------------------------------------------------
   void insert(int row, int column, real value)
   {
      // konsistenz prüfung
      if (row<0 || column<0 || row >= rows || column>=cols || value==0.0)
         return;

      // ungeprüfter insert
      Insert(row, column, value);
   }

   // -------------------------------------------------------------------------
   // removes an element
   // -------------------------------------------------------------------------
   real remove(int row, int column)
   {
      int tmpsize        = 0;
      int size           = sizes[row];
      real *row_values  = values[row];
      int *row_indices   = indices[row];
      real retvalue     = 0.0f;
      // copy row data
      for (int i=0; i<size; i++) {
         if (column==row_indices[i]) {
             retvalue = row_values[i];
         }
         else {
            tmpvalues[tmpsize] = row_values[i];
            tmpindices[tmpsize] = row_indices[i];
            tmpsize++;
         }
      }
      insert_tmp(row, tmpsize);
      // return removed value
      return retvalue;
   }

   // -------------------------------------------------------------------------
   // removes an element
   // -------------------------------------------------------------------------
   real get(int row, int column)
   {
      int size           = sizes[row];
      real *row_values  = values[row];
      int *row_indices   = indices[row];
      // search index
      for (int i=0; i<size; i++) {
         if (column==row_indices[i]) {
             return row_values[i];
         }
      }
      return 0.0f;
   }

   // -------------------------------------------------------------------------
   // copy tmp values into matrix
   // -------------------------------------------------------------------------
   void insert_tmp(int row, int size)
   {
      if (size <= sizes[row]) {
         memcpy(values[row], tmpvalues, sizeof(real) * size);
         memcpy(indices[row], tmpindices, sizeof(int) * size);
         sizes[row] = size;
      }
      else if (size > 0) {
         // clear old stuff
         if (values[row]) delete[] values[row];
         if (indices[row]) delete[] indices[row];
         // allocate new stuf
         values[row] = new real[size];
         indices[row] = new int[size];
         // copy elements
         memcpy(values[row], tmpvalues, sizeof(real) * size);
         memcpy(indices[row], tmpindices, sizeof(int) * size);
         sizes[row] = size;
      }
   }

   // -------------------------------------------------------------------------
   // subtract row src from dst row, applying a factor
   // -------------------------------------------------------------------------
   void sub(int dst, int src, real factor)
   {
      // set tmp size
      int tmpsize  = 0;
      int srcsize  = sizes[src];
      int dstsize  = sizes[dst];
      int srcindex = 0;
      int dstindex = 0;
      // set start value - index
      int *srcindices = indices[src];
      int *dstindices = indices[dst];
      real *srcvalues = values[src];
      real *dstvalues = values[dst];

      while (srcindex < srcsize && dstindex<dstsize) {
         if (srcindices[srcindex] < dstindices[dstindex]) {
            // enter src - value
            tmpvalues[tmpsize] = -factor * srcvalues[srcindex];
            tmpindices[tmpsize] = srcindices[srcindex];
            tmpsize++;
            srcindex++;
         }
         else if (dstindices[dstindex] < srcindices[srcindex]) {
            // enter dst value
            tmpvalues[tmpsize] = dstvalues[dstindex];
            tmpindices[tmpsize] = dstindices[dstindex];
            tmpsize++;
            dstindex++;
         }
         else {
            // enter tmp value
            tmpvalues[tmpsize] = dstvalues[dstindex-1] - factor * srcvalues[srcindex-1];
            tmpindices[tmpsize] = dstindices[dstindex];
            if (tmpvalues[tmpsize]!=0.0f) tmpsize++;
            srcindex++;
            dstindex++;
         }
      }
      // copy src - remainders
      while (srcindex < srcsize) {
         tmpvalues[tmpsize] = -factor * srcvalues[srcindex];
         tmpindices[tmpsize] = srcindices[srcindex];
         tmpsize++;
         srcindex++;
      }
      // copy dst - remainders
      while (dstindex<dstsize) {
         // enter dst value
         tmpvalues[tmpsize] = dstvalues[dstindex];
         tmpindices[tmpsize] = dstindices[dstindex];
         tmpsize++;
         dstindex++;
      }
      // copy into table
      insert_tmp(dst, tmpsize);
   }

   // -------------------------------------------------------------------------
   // subtract row src from dst row, applying a factor
   // -------------------------------------------------------------------------
   void eliminate(int dst, int src)
   {
      int srcsize  = sizes[src];
      int dstsize  = sizes[dst];
      // set tmp pointers
      real *tmpvp = tmpvalues;
      int *tmpip = tmpindices;
      // set start value - index
      int *srcindices = indices[src];
      int *dstindices = indices[dst];
      int *srcend = &srcindices[srcsize];
      int *dstend = &dstindices[dstsize];
      real *srcvalues = values[src];
      real *dstvalues = values[dst];
      // calculate factor
      real factor = *dstvalues / *srcvalues;
      // step ahead
      srcindices++;
      dstindices++;
      srcvalues++;
      dstvalues++;

      while (srcindices < srcend && dstindices<dstend) {
         if (*srcindices < *dstindices) {
            // enter src - value
            *tmpvp++ = -factor * *srcvalues++;
            *tmpip++ = *srcindices++;
         }
         else if (*dstindices < *srcindices) {
            // enter dst value
            *tmpvp++ = *dstvalues++;
            *tmpip++ = *dstindices++;
         }
         else {
            real value = *dstvalues++ - factor * *srcvalues++;
            if (value) {
               // enter tmp value
               *tmpvp++  = value;
               *tmpip++ = *dstindices;
            }
            dstindices++;
            srcindices++;
         }
      }
      // copy src - remainders
      while (srcindices < srcend) {
         *tmpvp++ = -factor * *srcvalues++;
         *tmpip++ = *srcindices++;
      }
      // copy dst - remainders
      while (dstindices<dstend) {
         // enter dst value
         *tmpvp++ = *dstvalues++;
         *tmpip++ = *dstindices++;
      }
      // copy into table
      insert_tmp(dst, tmpvp - tmpvalues);
   }

   // -------------------------------------------------------------------------
   // transponiert die matrix
   // -------------------------------------------------------------------------
   void transpose()
   {
      // Elemente anlegen
      int   *csizes   = new int[cols];
      int  **cindices = new int*[cols];
      real **cvalues  = new real*[cols];
      memset(csizes,   0, sizeof(int)    * cols);
      memset(cindices, 0, sizeof(int *)  * cols);
      memset(cvalues,  0, sizeof(real *) * cols);

      // Anzahl Elemente der neuen Matrix bestimmen
      for (int i=0; i<rows; i++) {
         int *rindices = indices[i];
         int  size     = sizes[i];

         for (int j=0; j<size; j++)
            csizes[*rindices++]++;
      }

      // Arrays für indices und werte anlegen
      for (int i=0; i<cols; i++) {
         if (csizes[i] > 0) {
            cindices[i] = new int[csizes[i]];
            cvalues[i]  = new real[csizes[i]];
         }
      }

      // Zaehler zurücksetzen
      memset(csizes, 0, sizeof(int) * cols);

      // Werte kopieren, die Matrix bleibt dabei
      // richtig sortiert !!!!
      for (int i=0; i<rows; i++) {
         int  *rindices = indices[i];
         real *rvalues  = values[i];
         int   size     = sizes[i];

         for (int j=0; j<size; j++) {
            int col = *rindices++;

            int   idx        = csizes[col]++;
            int  *rowindices = cindices[col];
            real *rowvalues  = cvalues[col];
            
            rowindices[idx]  = i;
            rowvalues[idx]   = *rvalues++;
         }
         // Index und wert - listen löschen
         if (values[i])  delete[] values[i];
         if (indices[i]) delete[] indices[i];
      }

      // alte arrays freigeben
      delete[] values;
      delete[] indices;
      delete[] sizes;

      // neue werte setzen
      sizes   = csizes;
      values  = cvalues;
      indices = cindices;

      // Anzahl spalten, zeilen tauschen
      swap(rows, cols);
   }

   // -------------------------------------------------------------------------
   // divide a row by a factor
   // -------------------------------------------------------------------------
   void divide(int row, real factor)
   {
      int size = sizes[row];
      real *row_values = values[row];
      //
      for (int i=0; i<size; i++)
         row_values[i] /= factor;
   }

   // -------------------------------------------------------------------------
   // search for the maximum element in this column
   // -------------------------------------------------------------------------
   int pivot_search(int column, int start, int end)
   {
      real max = 0.0f;
      int pivot = -1;
      // search pivot
      for (int i=start; i<end; i++) {
         // compare value to max value
         if (sizes[i] && indices[i][0]==column && fabs(values[i][0])>max) {
            // set new pivot
            max = fabs(values[i][0]);
            pivot = i;
            // TODO:
            // changed pivot search....
            return pivot;
         }
      }
      return pivot;
   }

   // -------------------------------------------------------------------------
   // inverts this matrix
   // -------------------------------------------------------------------------
   bool invert()
   {
      // add an identity matrix to the right half (gauss)
      for (int i=0; i<rows; i++) {
         insert(i, i + cols, 1.0f);
      }

      // iterate over all rows
      for (int i=0; i<rows-1; i++) {
         // eliminate all values from subsequent rows
         for (int j=i+1; j<rows; j++) {
            // eliminate value from row
            if (sizes[j] && (indices[j][0] == i)) {
               eliminate(j, i);
            }
         }
         printf("%d von %d\n", i, rows);
      }
/*
      // substitute back
      for (int i=m-1; i>=0; i--) {
         // select source - line
         int src = i;
         // traverse back
         for (int j=i-1; j>=0; j--) {
            // search corresponding element
            Sparse **element = &rows[j];
            while ((*element) && ((*element)->col<src->col))
               element = &((*element)->next_col);
            // backsubstitute
            if ((*element) && ((*element)->col==src->col))
               subrow(element, src, (*element)->value/src->value);
         }
         // divide row by first element
         divide(src, src->value);
         printf("%d von %d\n", i, m);
      }
      for (int i=0; i<inv.m; i++)
      {
         Sparse *loop = inv.rows[i];
         while (loop) {
            loop->col -= n;
            loop = loop->next_col;
         }
      }
      */
      return true;
   }

   // -------------------------------------------------------------------------
   // apply matrix so that vec = A * vec
   // -------------------------------------------------------------------------
   void vecmul(int dim, double *vec)
   {
      // allocate tmp vector
      double *v = new double[dim];

      for (int i=0; i<rows; i++) {
         v[i] = mulrowbycolumn(i, vec);
      }
      // copy result
      memcpy(vec, v, sizeof(real) * dim);
      delete[] v;
   }

   // -------------------------------------------------------------------------
   // apply matrix so that r = A * x
   // -------------------------------------------------------------------------
   inline void vecmul(double *r, double *x, int n)
   {
      for (int i=0; i<n; i++) r[i] = mulrowbycolumn(i, x);
   }

   // -------------------------------------------------------------------------
   // apply matrix so that r = AT * x
   // -------------------------------------------------------------------------
   inline void vectmul(double *r, double *x, int n)
   {
      // start with resetting the result
      memset(r, 0, sizeof(double) * n);
      // apply a transposed matrix, vector multiply
      for (int i=0; i<n; i++) {
         // get size, index, value pointer
         int size   = sizes[i];
         int *inp   = indices[i];
         real *vap = values[i];
         // apply multiply
         for (int j=0; j<size; j++) {
            r[*inp++] += *vap++ * x[i];
         }
      }
   }

   // -------------------------------------------------------------------------
   // apply matrix so that r = A * x
   // with A stored in symmetric format
   // -------------------------------------------------------------------------
   inline void vecmuls(double *r, double *x, int n)
   {
      // local variables 
      double *valuep;
      int    *indexp;
      int    i, j, size;

      // start with resetting the result
      memset(r, 0, sizeof(double) * n);
      
      // apply a transposed matrix, vector multiply
      for (i=0; i<n; i++) {
         // get size, index, value pointer
         size   = sizes[i];
         indexp = indices[i];
         valuep = values[i];

         // multiply diagonal element
         r[i] += *valuep++ * x[*indexp++];
         
         // multiply remaining elements
         for (j=1; j<size; j++) {
            r[i]       += *valuep * x[*indexp];
            r[*indexp] += *valuep * x[i];
            indexp++;
            valuep++;
         }
      }
   }

   // -------------------------------------------------------------------------
   // row by row multiply
   // -------------------------------------------------------------------------
   real mulrowbyrow(int i, int j)
   {
      // set start value - index
      int *iindices = indices[i];
      int *jindices = indices[j];
      int *iend = &iindices[sizes[i]];
      int *jend = &jindices[sizes[j]];
      real *ivalues = values[i];
      real *jvalues = values[j];
      real sum = 0.0f;

      while (iindices < iend && jindices<jend) {
         if (*iindices < *jindices) {
            ivalues++;
            iindices++;
         }
         else if (*jindices < *iindices) {
            *jvalues++;
            *jindices++;
         }
         else {
            sum += *jvalues++ * *ivalues++;
            jindices++;
            iindices++;
         }
      }
      return sum;
   }

   // -------------------------------------------------------------------------
   // transform matrix to A(transposed) * A
   // instead of multiplying matrices as row by column, we
   // can apply row by row and calculate only one half of the matrix
   // (symmetrie)
   // -------------------------------------------------------------------------
   void transform_ata(SMCR &ata)
   {
      for (int i=0; i<rows; i++) {
         int tmpsize = 0;

         for (int j=0; j<rows; j++) {
            real sum = mulrowbyrow(j, i);
            if (sum!=0.0f) {
               ata.tmpvalues[tmpsize] = sum;
               ata.tmpindices[tmpsize] = j;
               tmpsize++;
            }
         }
         ata.insert_tmp(i, tmpsize);
      }
   }

   // -------------------------------------------------------------------------
   // transform matrix to A(transposed) * A
   // instead of multiplying matrices as row by column, we
   // can apply row by row and calculate only one half of the matrix
   // (symmetrie)
   // -------------------------------------------------------------------------
   void transform_atas()
   {
      // Elemente anlegen
      int   *nsizes   = new int[rows];
      int  **nindices = new int*[rows];
      real **nvalues  = new real*[rows];
      memset(nsizes,   0, sizeof(int)    * rows);
      memset(nindices, 0, sizeof(int *)  * rows);
      memset(nvalues,  0, sizeof(real *) * rows);

      for (int i=0; i<rows; i++) {
         int tmpsize = 0;

         // Indices und Werte berechnen
         for (int j=i; j<rows; j++) {
            real sum = mulrowbyrow(j, i);
            if (sum!=0.0f) {
               tmpvalues[tmpsize]  = sum;
               tmpindices[tmpsize] = j;
               tmpsize++;
            }
         }

         // Indices und Werte kopieren
         if (tmpsize > 0) {
            nsizes[i]   = tmpsize;
            nindices[i] = new int[tmpsize];
            nvalues[i]  = new real[tmpsize];

            memcpy(nindices[i], tmpindices, sizeof(int)  * tmpsize);
            memcpy(nvalues[i],  tmpvalues,  sizeof(real) * tmpsize);
         }
      }

      // alte arrays löschen
      for (int i=0; i<rows; i++) {
         if (values[i]) delete[] values[i];
         if (indices[i]) delete[] indices[i];
      }
      // free memory
      delete[] values;
      delete[] indices;
      delete[] sizes;
      
      // neue werte setzen
      sizes   = nsizes;
      values  = nvalues;
      indices = nindices;

      // cols = rows nach transformation
      cols = rows;
   }

   // -------------------------------------------------------------------------
   // multiplies a row by a column
   // -------------------------------------------------------------------------
   inline double mulrowbycolumn(int row, double *col)
   {
      double sum          = 0.0;
      int   size          = sizes[row];
      double *row_values  = values[row];
      int   *row_indices  = indices[row];
      //
      for (int i=0; i<size; i++) {
         sum += *row_values++ * col[*row_indices++];
      }
      return sum;
   }

   // -------------------------------------------------------------------------
   // checks matrix symmetrie
   // -------------------------------------------------------------------------
   bool check_symmetric()
   {
      if (rows!=cols)
         return false;

      for (int i=0; i<rows; i++) {
         for (int j=0; j<i; j++) {
            if (get(i,j) != get (j,i))
               return false;
         }
      }
      return true;
   }

   // -------------------------------------------------------------------------
   // solve A * x = b
   // -------------------------------------------------------------------------
   double *solve_iterative(double *b, double tol, int itmax, int *iter, double *err)
   {
      double diff2, res, aii;
      double *x;
      int n;

      // init variables
      x = new double[rows];
      n = rows;

      // transformation of a
      for (int i=0; i<rows; i++) {
         // get diagonal element
         aii  = remove(i, i);
         x[i] = b[i] /= aii;
         divide(i, aii);
      }

      // set start condition
      *iter    = 0;
      // main loop
      while (++(*iter) < itmax) {
         // reset error
         *err = 1.0f;
         // apply gauss seidel
         for (int k=0; k<n; k++) {
            //real newvalue
            res   = b[k] - mulrowbycolumn(k, x);
            diff2 = x[k] - res;
            *err += diff2 * diff2;
            x[k]  = res;
         }
         // calculate error an check stop condition
         *err = sqrt(*err);
         if (*err < tol) break;
      }
      return x;
   }

   // -------------------------------------------------------------------------
   // calculates the (L2 - Norm)^2 for a vector
   // -------------------------------------------------------------------------
   inline static double l2norm2(int dim, double *x)
   {
      double result = 0.0;
      for (int i=0; i<dim; i++) result += x[i] * x[i];
      return result;
   }

   // -------------------------------------------------------------------------
   // calculates the L2 - Norm for a vector
   // -------------------------------------------------------------------------
   inline static double l2norm(int dim, double *x)
   {
      double l2 = l2norm2(dim, x);
      //
      return (ismatherr()) ? 0.0 : sqrt(l2);
   }

   // -------------------------------------------------------------------------
   // conjugate gradient solve of A * x = b
   // -------------------------------------------------------------------------
   double *cgiter(double *b, double tol, int itmax, int *iter, double *err)
   {
      // define local variables
      double alpha, bnorm, beta, rho, rhoold = 0.0;
      double *x, *r, *p, *q;
      int n = rows;

      // alloctae temporary vectors
      x = new double[n];
      r = new double[n];
      p = new double[n];
      q = new double[n];

      // initialize x
      for (int i=0; i<n; i++) x[i] = b[i]/get(i, i);

      // compute initial residual r = b - A * x
      for (int i=0; i<n; i++)
         r[i] = b[i] - mulrowbycolumn(i, x);
      
      // set iteration conditions
      *iter = 0;
      bnorm = l2norm(n, b);
 
      // get rho als L2Norm ^ 2
      rho = l2norm2(n, r);      
 
      while (++(*iter) < itmax) {

          // check for first step
         if (*iter == 1) {
            // p = r
            for (int i=0; i<n; i++) p[i] = r[i];
         } 
         else {
            beta = rho / rhoold;
            // p = r + beta * p
            for (int i=0; i<n; i++) p[i] = r[i] + beta * p[i];
         }
         // q = A * p
         // alpha = rho / p * q
         alpha = 0.0;
         for (int i=0; i<n; i++) {
            q[i] = mulrowbycolumn(i, p);
            alpha += p[i] * q[i];
         }
         alpha = rho / alpha;
         // x = x + alpha * p
         // r = r - alpha * q
         for (int i=0; i<n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
         }
         rhoold = rho;

         // get new rho als L2Norm ^ 2
         rho = l2norm2(n, r);

         // calculate error and check stopping criterion.
         *err = sqrt(rho) / bnorm;
         if (*err <= tol) break;
      }
      delete[] r;
      delete[] p;
      delete[] q;

      return x;
   }

   // -------------------------------------------------------------------------
   // bestimmt das abbruchkriterium für den conjugate gradient algorithmus
   // -------------------------------------------------------------------------
   bool result(double bnorm, double rnorm, double eps)
   {
      return (rnorm < eps * bnorm || (IsDZero(bnorm) && IsDOne(1.0 + rnorm)));
   }

   // -------------------------------------------------------------------------
   // conjugate gradient solve of A * x = b
   // with A stored in symmetric mode
   // -------------------------------------------------------------------------
   double *cgiters(double *b, double tol, int itmax, int *iter, double *err)
   {
      // define local variables
      double alpha, bnorm, beta, rho, rhoold = 0.0;
      double *x, *r, *p, *q;
      int n = rows;

      // alloctae temporary vectors
      x = new double[n];
      r = new double[n];
      p = new double[n];
      q = new double[n];

      // reset math error
      setmatherrflag(false);

      // initialize x
      for (int i=0; i<n; i++) x[i] = b[i]/values[i][0];

      // compute initial residual r = b - A * x
      vecmuls(r, x, n);
      for (int i=0; i<n; i++) r[i] = b[i] - r[i];
      
      // setup for iteration
      *iter = 0;
      bnorm = l2norm(n, b);
      rho   = pow(l2norm(n, r), 2.0);      
 
      while (!result(bnorm, l2norm(n, r), tol) && (++(*iter) < itmax)) {
      	 // get last error
         if (!IsDZero(bnorm)) *err = l2norm(n, r) / bnorm;

          // check for first step
         if (*iter == 1) {
            // p = r
            for (int i=0; i<n; i++) p[i] = r[i];
         } 
         else {
            beta = rho / rhoold;
            // p = r + beta * p
            for (int i=0; i<n; i++) p[i] = r[i] + beta * p[i];
         }
         // q = A * p
         vecmuls(q, p, n);
         // alpha = rho / p * q
         alpha = 0.0;
         for (int i=0; i<n; i++) alpha += p[i] * q[i];
         alpha = rho / alpha;
         // x = x + alpha * p
         // r = r - alpha * q
         for (int i=0; i<n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
         }
         rhoold = rho;

         // get new rho als L2Norm ^ 2
         rho = pow(l2norm(n, r), 2.0);

         // calculate error and check stopping criterion.
         //*err = sqrt(rho) / bnorm;
         //if (ismatherr() || *err <= tol) break;
      }
      delete[] r;
      delete[] p;
      delete[] q;

      return x;
   }

   // -------------------------------------------------------------------------
   // bcgiter help function
   // -------------------------------------------------------------------------
   inline void asolve(int n, real *b, real *x, int itrnsp)
   {
      for(int i=0;i<n;i++) x[i]=(sa[i] != 0.0 ? b[i]/sa[i] : b[i]);
      // The matrix A~ is the diagonal part of A, stored in the first n 
      // elements of sa. Since the
      // transpose matrix has the same diagonal, the flag itrnsp is not used.
   }

   // -------------------------------------------------------------------------
   // bcgiter help function
   // input x
   // output r
   // itrnsp = 0, Matrix
   // itrnsp = 1, Transposed
   // -------------------------------------------------------------------------
   inline void atimes(int n, real *x, real *r, int itrnsp)
   {
      if (itrnsp) vectmul(r, x ,n);
      else        vecmul(r, x, n); 
   }

   // -------------------------------------------------------------------------
   // bcgiter help function
   // Compute one of two norms for a vector sx[1..n],
   // as signaled by itol. Used by bcgiter
   // -------------------------------------------------------------------------
   inline static double snrm(int n, double sx[], int itol)
   {
      int isamax;
      double ans;
      if (itol <= 3) {
         ans = 0.0;
         // Vector magnitude norm.
         //for (int i=0;i<n;i++) ans += sx[i]*sx[i];
         //return sqrt(ans);
         return snrm2(n, sx);
      }
      else {
         isamax=1;
         // Largest component norm.
         for (int i=0;i<n;i++) {
            if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
         }
         return fabs(sx[isamax]);
      }
   }

   // -------------------------------------------------------------------------
   // Performs euclidean norm computation,
   // avoiding floating point errors
   // 
   // @param af
   // @return
   // -------------------------------------------------------------------------
  inline static double snrm2(int n, double sx[])
  {
    if (n == 1)
      return abs(sx[0]);

    double scale = 0.0;
    double ssq   = 1.0;

    for (int i = 0; i < n; i++)
    {
      double absx = sx[i];

      if (absx != 0.0)
      {
        if (absx < 0.0)
          absx = -absx;

        if (scale < absx)
        {
          double tmp = scale / absx;
          ssq        = 1.0 + ssq * tmp * tmp;
          scale      = absx;
        }
        else
        {
          double tmp = absx / scale;
          ssq       += tmp * tmp;
        }
      }
    }
    return scale * sqrt(ssq);
  }


   // -------------------------------------------------------------------------
   // bijective gradient solve of A * x = b
   // -------------------------------------------------------------------------
   double *bcgiter(double *b, int itol, double tol, int itmax, int *iter, double *err)
   {
      double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
      double *p,*pp,*r,*rr,*z,*zz, *x;
      int n = rows;
      double EPS = (real)1e-8;

      p  = new double[n];
      pp = new double[n];
      r  = new double[n];
      rr = new double[n];
      z  = new double[n];
      zz = new double[n];
      sa = new double[n];
      x  = new double[n];

      // store diagonal
      for (int i=0; i<n; i++) sa[i] = get(i, i);

      // Calculate initial residual.
      *iter=0;
      atimes(n,x,r,0);
      for (int j=0;j<n;j++) {
         r[j]=b[j]-r[j];
         rr[j]=r[j];
      }

      // uncomment for minimum residual variant
      // atimes(n,r,rr,0);
      if (itol == 1) {
         bnrm=snrm(n,b,itol);
         asolve(n,r,z,0);
      }
      else if (itol == 2) {
         asolve(n,b,z,0);
         bnrm=snrm(n,z,itol);
         asolve(n,r,z,0);
      }
      // else if (itol == 3 || itol == 4)
      else {
         asolve(n,b,z,0);
         bnrm=snrm(n,z,itol);
         asolve(n,r,z,0);
         znrm=snrm(n,z,itol);
      } 

      // Main loop.
      while (*iter < itmax) { 
         ++(*iter);
         // Final 1 indicates use of transpose matrix AT
         asolve(n,rr,zz,1);

         bknum=0.0;
         for (int j=0;j<n;j++) bknum += z[j]*rr[j];
         // Calculate coefficient bk and direction vectors p and pp.
         if (*iter == 1) {
            for (int j=0;j<n;j++) {
               p[j]=z[j];
               pp[j]=zz[j];
            }
         }
         else {
            bk=bknum/bkden;
            for (int j=0;j<n;j++) {
               p[j]=bk*p[j]+z[j];
               pp[j]=bk*pp[j]+zz[j];
            }
         }

         // Calculate coecient ak, new iterate x, and new
         // residuals r and rr. 
         bkden=bknum;
         atimes(n,p,z,0);
         akden=0.0;
         for (int j=0;j<n;j++) akden += z[j]*pp[j];
         ak=bknum/akden;
         atimes(n,pp,zz,1);

         for (int j=0;j<n;j++) {
            x[j] += ak*p[j];
            r[j] -= ak*z[j];
            rr[j] -= ak*zz[j];
         }
         // Solve A~ * z = r and check stopping criterion.
         asolve(n,r,z,0);
         if (itol == 1)
            *err=snrm(n,r,itol)/bnrm;
         else if (itol == 2)
            *err=snrm(n,z,itol)/bnrm;
         else if (itol == 3 || itol == 4) {
            zm1nrm=znrm;
            znrm=snrm(n,z,itol);
            if (fabs(zm1nrm-znrm) > EPS*znrm) {
               dxnrm=fabs(ak)*snrm(n,p,itol);
               *err=znrm/fabs(zm1nrm-znrm)*dxnrm;
            } else {
               // Error may not be accurate, so loop again.
               *err=znrm/bnrm;
               continue;
            }
            xnrm=snrm(n,x,itol);
            if (*err <= 0.5*xnrm) *err /= xnrm;
            else {
               // Error may not be accurate, so loop again.
               *err=znrm/bnrm;
               continue;
            }
         }
         // printf("iter=%4d err=%12.6f\n",*iter,*err);
         if (*err <= tol) break;
      }

      delete[] p;
      delete[] pp;
      delete[] r;
      delete[] rr;
      delete[] z;
      delete[] zz;
      delete[] sa;

      return x;
   }

   // -------------------------------------------------------------------------
   // solve A * x = b
   // -------------------------------------------------------------------------
   double *sor(double *b, double omega, int itol, double tol, int itmax, int *iter, double *err)
   {
      double diff2, res, factor1, factor2;
      double *x, *rowvalues;
      int n, size;
      int *rowindices;

      // set vectors
      n  = rows;
      x  = new double[rows];
      sa = new double[n];

      // calculate initial guess
      for (int i=0; i<rows; i++) {
         sa[i] = get(i, i);
         x[i]  = b[i] / sa[i];
      }

      *iter    = 0;
      // main loop
      while (++(*iter) < itmax) {
         // reset error
         *err = 0.0;
         // apply gauss seidel
         for (int i=0; i<n; i++) {
           res     = ( 1.0 - omega ) * x[i];
           factor1 = omega / sa[i];
           factor2 = b[i];

           size       = sizes[i];
           rowvalues  = values[i];
           rowindices = indices[i];
           //
           for (int j=0; j<size; j++) {
              if (*rowindices!=i) {
                 factor2 -= *rowvalues * x[*rowindices];
              }
              rowvalues++;
              rowindices++;
           }
           res += factor1 * factor2;

           diff2 = x[i] - res;
           *err += diff2 * diff2;
           x[i]  = res;
         }
         // calculate error and check stop criteria
         *err = sqrt(*err);
         if (*err < tol) break;
      }
      // delete tmporary arrays
      delete[] sa;

      return x;
   }

   // -------------------------------------------------------------------------
   // print matrix
   // -------------------------------------------------------------------------
   void print()
   {
	   int i, j;

	   for (i=0; i<rows; i++){
         printf("%4d (%4d): ", i, sizes[i]);
		   real *valuep = values[i];
         int *indexp = indices[i];
         int size = sizes[i];
         int *indexend = &indexp[size];

         for (j=0; j<cols; j++) {
            if (!valuep || j<*indexp || indexp>=indexend) {
               //printf("\t0.00");
            }
            else {
               printf("\t%4d : %1.2lf", j, *valuep);
               valuep++;
               indexp++;
            }
         }
		   printf("\n");
	   }
      printf("\n");
   }

};


#endif