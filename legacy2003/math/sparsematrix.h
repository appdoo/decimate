// ============================================================================
//
// Class SMat : Sparse Matrix
//
// Definition and Implementation
//
// $Author: Andor $
// $Date: 10.12.03 0:55 $
// $Revision: 4 $
//
// ============================================================================
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

// C++ Includes
#include <stdio.h>
#include <windows.h>

// Local Includes
#include "math3D.h"

typedef double real;


typedef struct Sparse
{
public:
   // value of the sparse element
   real value;
   // position of the element
   int row;
   int col;
   // next row / col elements
   struct Sparse *next_row;
   struct Sparse *next_col;
} Sparse;


/******************************************************************************
 * sparse matrix implementation
 *****************************************************************************/
class SMat
{
public:
   // number of rows
   int m;
   // number of columns
   int n;
   // Pointer to row and column elements   
   Sparse ** rows;
   Sparse ** cols;

public:
   /***************************************************************************
    * construct a standard m by n sparse matrix
    *
    * m rows
    * n columns
    *
    **************************************************************************/
   SMat(int m, int n)
   {
      init(m, n);
   }

   /***************************************************************************
    * copy constructor implementation
    **************************************************************************/
   SMat(const SMat & mat)
   {
      // initialization
      init(mat.m, mat.n);
      // iterate over all rows
      copy(mat);
   }

   /***************************************************************************
    * initialization
    **************************************************************************/
   void init(int m, int n)
   {
      // set dimension
      this->m = m;
      this->n = n;
      // allocate base reference
      rows = new Sparse *[m];
      cols = new Sparse *[n];
      // clear 
      memset(rows, 0, sizeof(Sparse *) * m);
      memset(cols, 0, sizeof(Sparse *) * n);
   }

   /***************************************************************************
    * copy function
    **************************************************************************/
   void copy(const SMat & mat)
   {
      // iterate over all rows
      for (int i=0; i<m;i++) {
         Sparse *element = mat.rows[i];
         // iterate over all columns
         while (element) {
            insert(i, element->col, element->value);
            element = element->next_col;
         }
      }
   }

   /***************************************************************************
    * deletes the  sparse matrix
    **************************************************************************/
   ~SMat()
   {
      clear();
      //
      if (rows) delete[] rows;
      if (cols) delete[] cols;
   }

   /***************************************************************************
    * clear all elements from the matrix
    **************************************************************************/
   void clear()
   {
      // iterate over all rows
      for (int i=0; i<m;i++) {
         Sparse *element = rows[i];
         // iterate over all columns
         while (element) {
            Sparse *garbage = element;
            element = element->next_col;
            delete garbage;
         }
      }
      // clear row and column pointers
      memset(rows, 0, sizeof(Sparse *) * m);
      memset(cols, 0, sizeof(Sparse *) * n);   
   }

   /***************************************************************************
    * set elements from the matrix to identity
    **************************************************************************/
   void identity()
   {
      clear();
      // iterate over all rows
      for (int i=0; i<m;i++) {
         insert(i, i, 1.0f);
      }
   }

   /***************************************************************************
    * get element at specified position
    **************************************************************************/
   Sparse *get(int row, int column)
   {
      // find our element in rows and columns
      Sparse **left = &rows[row];
      Sparse **top  = &cols[column];
      //
      // row search
      while ((*left) && ((*left)->col < column))
         left = &((*left)->next_col);
      // check if element exists
      if ((*left) && ((*left)->col == column))
         return *left;
      // column search
      while ((*top) && ((*top)->row < row))
         top = &((*top)->next_row);
      // create a new element
      Sparse *element = new Sparse;
      // initialize it
      element->row = row;
      element->col = column;
      element->next_row = *top;
      element->next_col = *left;
      // update list
      *left = element;
      *top = element;
      // return element
      return element;
   }

   /***************************************************************************
    * get element at specified position
    **************************************************************************/
   real getvalue(int row, int column)
   {
      // find our element in rows and columns
      Sparse *left = rows[row];
      //
      // row search
      while (left && left->col < column)
         left = left->next_col;
      // check if element exists
      if (left && left->col==column)
         return left->value;
      // return element
      return 0.0f;
   }

   // -------------------------------------------------------------------------
   // multipliziert eine sparse reihe mit einem vektor
   // -------------------------------------------------------------------------
   real mul(int rowid, real *vec)
   {
      // variablen initialisieren
      real result = 0.0;
      Sparse *row = rows[rowid];

      while (row) {
         result += row->value * vec[row->col];
         row     = row->next_col;
      }
      return result;
   }

   // -------------------------------------------------------------------------
   // multipliziert die sparse matrix mit einem vektor
   // -------------------------------------------------------------------------
   void mul(real **dst, real *vec)
   {
      // variablen initialisieren
      real *tmp   = new real[m];

      for (int i=0; i<m; i++) {
         tmp[i] = mul(i, vec);
      }
      if (*dst) delete[] *dst;
      *dst = tmp;
   }

   /***************************************************************************
    * inserts an element
    **************************************************************************/
   void insert(int row, int column, real value)
   {
      // validate, !value is for 0 values -> null for int, real, and double
      if (row<0 || column<0 || row >= m || column>=n || value==0.0f)
         return;
      // get element at specified position
      get(row, column)->value = value;
   }


   /***************************************************************************
	 * multiply by row and column
	 **************************************************************************/
   real mulrowbycolumn(Sparse *row, Sparse *col) 
   {
	   real result = 0.0f;

	   while (row && col) {
         if (row->col == col->row){
			   result += row->value * col->value;
            row = row->next_col;
            col = col->next_row;
		   }
         else  {
            if (row->col < col->row) { 
               row = row->next_col;
            }
            else {
               col = col->next_row;
            }
         }
	   }
	   return result;
   }

   real mulrowbycolumn(Sparse *row, Sparse *col, int rowbreak) 
   {
	   real result = 0.0f;

      while (row && col && row->col<rowbreak && col->row<rowbreak) {
         if (row->col == col->row){
			   result += row->value * col->value;
            row = row->next_col;
            col = col->next_row;
		   }
         else  {
            if (row->col < col->row) { 
               row = row->next_col;
            }
            else {
               col = col->next_row;
            }
         }
	   }
	   return result;
   }

   /***************************************************************************
	 * multiply by row and column
	 **************************************************************************/
   static inline real mulrowbycolumn(real *row, Sparse *col) 
   {
	   real result = 0.0f;

	   while (col) {
         result += row[col->row] * col->value;
         col = col->next_row;
	   }
	   return result;
   }

   /***************************************************************************
	 * define matrix multiplication
	 **************************************************************************/
	SMat operator *(const SMat &m2)
	{
      SMat result(m, m2.n);
      real *frows = new real[n];

      for (int i=0; i<result.m; i++) {
         // copy row data for fast access
         memset(frows, 0, sizeof(real) * n);
         Sparse *loop = rows[i];
         while (loop) {
            frows[loop->col] = loop->value;
            loop = loop->next_col;
         }
         // get reference to result row
         Sparse **element = &result.rows[i];
		   for(int j=0; j<result.n; j++) {
            // claculate value
            real value = mulrowbycolumn(frows, m2.cols[j]);
            // check for valid value
            if (value != 0.0f) {
               // insert it
               Sparse *sp   = new Sparse;
               sp->col      = j;
               sp->next_col = NULL;
               sp->next_row = NULL;
               sp->row      = i;
               sp->value    = value;
               // set pointers
               *element = sp;
               element = &(sp->next_col);
            }
         }
		}
      delete[] frows;
      return result;
   }

   /***************************************************************************
	 * define matrix multiplication case Atransposed * A
	 **************************************************************************/
	SMat ata()
	{
      SMat result(n, n);

      for (int i=0; i<result.m; i++) {
         // get reference to result row
         Sparse **element = &result.rows[i];
		   for(int j=0; j<result.n; j++) {
            // claculate value
            real value = mulrowbycolumn(rows[j], cols[i]);
            // check for valid value
            if (value != 0.0f) {
               // insert it
               Sparse *sp   = new Sparse;
               sp->col      = j;
               sp->next_col = NULL;
               sp->next_row = NULL;
               sp->row      = i;
               sp->value    = value;
               // set pointers
               *element = sp;
               element = &(sp->next_col);
            }
         }
		}
      return result;
   }

   // -------------------------------------------------------------------------
   // transponieren der matrix
   // -------------------------------------------------------------------------
   SMat &transpose()
   {
      // next pointer und indices tauschen
      for (int i=0; i<m;i++) {
         Sparse *element = rows[i];
         // iterate over all columns
         while (element) {
            swap(element->col, element->row);
            swap(element->next_col, element->next_row);
            // next_col bereits getauscht
            element = element->next_row;
         }
      }
      // dimensionen und zeiger tauschen
      swap (n, m);
      swap (rows, cols);

      return *this;
   }
   
   /***************************************************************************
    * applies an unary - operator
    **************************************************************************/
   void neg()
   {
      // iterate over all rows
      for (int i=0; i<m;i++) {
         Sparse *element = rows[i];
         // iterate over all columns
         while (element) {
            element->value = -element->value;
            element = element->next_col;
         }
      }
   }

   /***************************************************************************
    * divide row
    **************************************************************************/
   void divide(Sparse *row, real divider)
   {
      while (row) {
         // divide
         row->value /= divider;
         // next element
         row = row->next_col;
      }
   }

   /***************************************************************************
    * divide row
    **************************************************************************/
   void divide(int row, real divider)
   {
      divide(rows[row], divider);
   }

   /***************************************************************************
    * subtract source from row
    **************************************************************************/
   static void subrow(Sparse **element, Sparse *src, real modifier)
   {
      // subtract src element from elements
      while (src) {
         // skip to element
         while (*element && (*element)->col<src->col) {
            element = &((*element)->next_col);
         }
         // check if weve found the right element
         if (!(*element) || (*element)->col>src->col) 
         {
            // insert a new element
            Sparse *newelem    = new Sparse;
            newelem->col       = src->col;
            newelem->next_col  = *element;
            newelem->value     = -modifier * src->value;
            *element           = newelem;
         }
         // check if we are equal
         else if (src->col == (*element)->col) 
         {
            (*element)->value -= modifier * src->value;
         }
         // check if we need to delete the element
         if (fabs((*element)->value) == 0.0f) {
            // delete element
            Sparse *garbage  = *element;
            *element = (*element)->next_col;
            delete garbage;
         }
         // next element
         src = src->next_col;
      }
   }

   /***************************************************************************
    * gauss elimination
    **************************************************************************/
   void eliminate(Sparse **element, Sparse *src)
   {
      // check column indices
      if ((*element)->col != src->col)
         return;
      // calculate modifier
      real modifier = (*element)->value / src->value;
      // discard first element
      Sparse *garbage = *element;
      *element = (*element)->next_col;
      delete garbage;
      // skip first src element
      src = src->next_col;
      // sub lines
      subrow(element, src, modifier);
   }

   /***************************************************************************
    * search for the maximum element in this column
    **************************************************************************/
   int pivot_search(Sparse **rows, int column, int start, int end)
   {
      real max = 0.0f;
      int pivot = start;
      // search pivot
      for (int i=start; i<end; i++) {
         // compare value to max value
         if (rows[i] && rows[i]->col==column && fabs(rows[i]->value)>max) {
            // set new pivot
            max = fabs(rows[i]->value);
            pivot = i;
            // TODO:
            // changed pivot search....
            return pivot;
         }
      }
      printf("No pivot found...\n");
      return pivot;
   }

   /***************************************************************************
    * calculate the inverse matrix
    **************************************************************************/
   SMat inverse(SMat *result)
   {
      // define a mtrix for inserting
      SMat inv(m, n * 2);
      // copy the elements
      inv.copy(*this);
      // add an identity matrix to the right half (gauss)
      for (int i=0; i<m; i++) {
         inv.insert(i, i + n, 1.0f);
      }
      // get a local reference to the rows
      Sparse **rows = inv.rows;
      // iterate over all rows
      for (int i=0; i<m-1; i++) {
         // search pivot
         int pivot = pivot_search(rows, i, i, m);
         swap(rows[i], rows[pivot]);
         // check if we have a value
         if (!rows[i] || rows[i]->col!=i || fabs(rows[i]->value)==0.0f) {
            printf("inversion error\n");
            return inv;
         }
         // eliminate all values from subsequent rows
         for (int j=i+1; j<m; j++) {
            // eliminate value from row
            if ((rows[j]) && (rows[j]->col == i) && (rows[j]->value)) {
               eliminate(&rows[j], rows[i]);
            }
         }
         printf("%d von %d\n", i, m);
      }
      // substitute back
      for (int i=m-1; i>=0; i--) {
         // select source - line
         Sparse *src = rows[i];
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
      return inv;
   }

   /***************************************************************************
    * solve A * x = b
    **************************************************************************/
   SMat solve_iterative(SMat *srcb)
   {
      // define a result vector
      SMat a(m, n);
      SMat b(*srcb);
      SMat x(b.m, 1);

      if (m != b.m || n!=b.m)
         printf ("matrix error A[%d, %d] b[%d, %d]\n", m, n, b.m, b.n); 

      printf("solve, transforming b\n");

      // transformation of b
      for (int i=0; i<b.m; i++) {
         real aii = getvalue(i, i);
         if (aii==0.0f)
            printf("trans b, divider error\n");
         if (b.rows[i]) {
            //printf("%lf\n", b.rows[i]->value);
            b.rows[i]->value /= aii;
         }
         else
         {
            b.insert(i, 0, 1.0f);
            b.rows[i]->value = 0.0f;
         }
         x.insert(i, 0, 1.0f);
         x.rows[i]->value = 0.0f;//b.rows[i]->value;
      }
      printf("solve, transforming a\n");

      // transformation of a
      for (int i=0; i<m; i++) {
         // get diagonal element
         real aii = getvalue(i, i);
         if (aii==0.0f)
            printf("trans a, divider error\n");
         for (int k=0; k<n; k++) {
            if (k!=i) {
               a.insert(i, k, getvalue(i, k) / aii);
            }
         }
      }
      printf("solving system\n");
      int maxiter = 10000;
      int iter = 0;
      real conv = 1.0E-4f;
      real diff = 0.0f;

      do {
         // set difference to 0
         diff = 0.0f;
         // apply gauss seidel
         for (int k=0; k<x.m; k++) {
            real newvalue = b.getvalue(k, 0) - mulrowbycolumn(a.rows[k], x.cols[0]);
            real diff2 = x.rows[k]->value - newvalue;
            diff = diff + diff2 * diff2;
            x.rows[k]->value = newvalue;
         }
         diff = sqrt(diff);
         printf("%lf\n", diff);
      } while (iter++<maxiter && diff>conv);
      // return result
      return x;
   }

   /***************************************************************************
    * L - R Triangulation
    **************************************************************************/
   void lrTriangulate(SMat &l, SMat &r)
   {
      // clear results
      l.identity();
      r.identity();

      // copy first row
      Sparse *loop = rows[0];
      while (loop) {
         r.insert(0, loop->col, loop->value);
         loop = loop->next_col;
      }
      // triangulate
      for (int i=0; i<m; i++) {
         printf ("triangulate %d von %d\n", i, m);
         // process l - elements
         for (int k=0; k<i; k++) {
            real value = getvalue(i, k);

            value -= mulrowbycolumn(l.rows[i], r.cols[k], k);
            /*
            loop = l.rows[i];
            while (loop && loop->col<k) {
               value -= loop->value * r.getvalue(loop->col, k);
               loop = loop->next_col;
            }
            */
            l.insert(i, k, value / r.getvalue(k, k));
         }

         // process k - elements
         for (int k=i; k<m; k++) {
            real value = getvalue(i, k);

            value -= mulrowbycolumn(l.rows[i], r.cols[k], i);
            /*
            loop = l.rows[i];
            while (loop && loop->col<i) {
               value -= loop->value * r.getvalue(loop->col, k);
               loop = loop->next_col;
            }
            */
            r.insert(i, k, value);
         }
      }
   }


   /***************************************************************************
    * print matrix
    **************************************************************************/
   void print()
   {
	   int i, j;
	   Sparse * p;

	   for (i=0; i<m; i++){
		   p=rows[i];
		   for (j=0; j<n; j++)
			   if ((!p)||(p->col>j))
				   printf("\t0.00");
			   else{
				   printf("\t%1.2lf", p->value);
               p=p->next_col;
			   }
		   printf("\n");
	   }
      printf("\n");
   }
};

#endif