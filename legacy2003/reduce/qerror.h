//=============================================================================
//
// QError : QError Class Definition
//
// Definiert eine Klasse QError, welche eine Matrix zur Bestimmung
// der Quadrierten Abstände eines Punktes zu einer / mehrerer Ebenen
// aufnimmt
//
// $Author: Andor $
// $Date: 11.12.03 0:16 $
// $Revision: 4 $
// 
//=============================================================================
#ifndef QERROR_H
#define QERROR_H

#include "..\math\math3d.h"

namespace decimate {

class QError
{
private:
   // rechte, obere Hälfte der
   // symmetrischen Abstandsmatrix
   double aa, ab, ac, ad;
   double     bb, bc, bd;
   double         cc, cd;
   double             dd;

public:
   QError();
   QError(QError &q1, QError &q2);

   void   Init();
   void   Init(QError &q1, QError &q2);

   void   SetPlane(Vec3 &v1, Vec3 &v2, Vec3 &v3);
   void   Add(QError &qerror);
   double GetError(Vec3 &p);
   bool   Optimize(Vec3 &p);
};

} // end namespace
#endif
