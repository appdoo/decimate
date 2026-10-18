//=============================================================================
//
// CConformalMapper : Class Definition
//
// Abbildung von Charts in das TexturkoordinatenSystem
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 5 $
// 
//=============================================================================
#ifndef CCONFORMALMAPPER_H
#define CCONFORMALMAPPER_H

// local includes
// Local Includes
#include "..\math\math3D.h"
#include "..\math\halfedge.h"
#include "..\math\smcr.h"
#include "..\math\algorithm.h"
#include "..\utility.h"

#include "cchart.h"
#include "cchartset.h"



namespace chartification {

class CConformalMapper
{
   // Vec3 coordinate Set to Navigate on
   Vec3 *points;

public:
   // Constructs a CConformalMapper class
   CConformalMapper(Vec3 *points);

   // Deletes a CConformalMapper
   ~CConformalMapper();

   // Maps a chart
   void Parametrize(CChart *chart);

   // Maps a whole chartset
   void Parametrize(CChartSet *charts);

private:
   // calculates chart geometric diameter
   float Diameter(int *p1, int *p2, CChart *chart);

   // returns the triangle gradient
   void Gradient(Complex *w1, Complex *w2, Complex *w3, float *dT, Face *f, Vertex *vertices);

   // chooses the pinned coordinates
   int InitPinned(CChart *chart);

   // Matrizen a und B initialisieren
   void InitSMCR(SMCR &am, SMCR &bm, int numPinned, CChart *chart);

   // inserts a complex value into the matrix
   void Insert(SMCR *mat, Complex *c, int i, int j, int rows2, int columns2);
};

} // end namespace
#endif
