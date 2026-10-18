//=============================================================================
//
// QError : QError Implementation
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

#include "qerror.h"

namespace decimate {

// ----------------------------------------------------------------------------
// Erstellt eine leere Fehlermatrix
// ----------------------------------------------------------------------------
QError::QError()
{
   Init();
}

// ----------------------------------------------------------------------------
// Erstellt eine Fehlermatrix aus der Addition der übergebenen Matrizen
// ----------------------------------------------------------------------------
QError::QError(QError &q1, QError &q2)
{
   Init(q1, q2);
}

// ----------------------------------------------------------------------------
// cleares the error matrix
// ----------------------------------------------------------------------------
void QError::Init()
{
   // Fehlermatrix mit 0 initialisieren
   aa = ab = ac = ad = 0.0f;
        bb = bc = bd = 0.0f;
             cc = cd = 0.0f;
                  dd = 0.0f;
}


// ----------------------------------------------------------------------------
// initializes to the sum of both matrices
// ----------------------------------------------------------------------------
void QError::Init(QError &q1, QError &q2)
{
  aa = q1.aa + q2.aa;
  ab = q1.ab + q2.ab;
  ac = q1.ac + q2.ac;
  ad = q1.ad + q2.ad;

  bb = q1.bb + q2.bb;
  bc = q1.bc + q2.bc;
  bd = q1.bd + q2.bd;

  cc = q1.cc + q2.cc;
  cd = q1.cd + q2.cd;

  dd = q1.dd + q2.dd;
}


// ----------------------------------------------------------------------------
// creates plane equation from the given coordinates
// ----------------------------------------------------------------------------
void QError::SetPlane(Vec3 &p1, Vec3 &p2, Vec3 &p3)
{
   // vektoren ermitteln 
   Vec3 a = p2 - p1;
   Vec3 b = p3 - p1;

   // normalenvektor auf die flaeche bestimmen
   Vec3 nrm = a ^ b;
   nrm.norm();

   // distanz - faktor ermitteln
   double d = -(nrm.x * p1.x + nrm.y * p1.y + nrm.z * p1.z);

   // matrix werte eintragen
   aa = nrm.x * nrm.x;
   ab = nrm.x * nrm.y;
   ac = nrm.x * nrm.z;
   ad = nrm.x * d;

   bb = nrm.y * nrm.y;
   bc = nrm.y * nrm.z;
   bd = nrm.y * d;

   cc = nrm.z * nrm.z;
   cd = nrm.z * d;

   dd = d * d;
}


// ----------------------------------------------------------------------------
// adds the error matrix
// ----------------------------------------------------------------------------
void QError::Add(QError &qerror)
{
  aa += qerror.aa;
  ab += qerror.ab;
  ac += qerror.ac;
  ad += qerror.ad;

  bb += qerror.bb;
  bc += qerror.bc;
  bd += qerror.bd;

  cc += qerror.cc;
  cd += qerror.cd;

  dd += qerror.dd;
}


// ----------------------------------------------------------------------------
// liefert den Fehler fuer einen Punkt p
//
//
// Der Fehler errechnet sich anhand, und ergibt die quadrierte Summe
// der Abstände von p zu allen Flächen wieder
//
//     | a2 ab ac ad |
//     | ab b2 bc bd |     (T)
// p * | ac bc c2 cd | * p
//     | ad bd cd d2 |
//
//
// ----------------------------------------------------------------------------
double QError::GetError(Vec3 &p)
{
   double px = p.x;
   double py = p.y;
   double pz = p.z;

   double x = aa * px + ab * py + ac * pz + ad;
   double y = ab * px + bb * py + bc * pz + bd;
   double z = ac * px + bc * py + cc * pz + cd;
   double w = ad * px + bd * py + cd * pz + dd;

   return x * px + y * py + z * pz + w;
}


// ----------------------------------------------------------------------------
// versucht den Fehler des angegebenen Punktes durch Positiosänderung zu 
// verbessern
//
//     | a2 ab ac ad | inverse      
//     | ab b2 bc bd |                       T
//     | ac bc c2 cd |           * [0 0 0 1]
//     | 0  0  0   1 |
//
// ----------------------------------------------------------------------------
bool QError::Optimize(Vec3 &p)
{
   double m[4][4];
   double inverse[4][4];

	// Werte in die matrix eintragen
	m[0][0] = aa;
	m[1][1] = bb;
	m[2][2] = cc;
	//
	m[1][0] = m[0][1]= ab;
	m[2][0] = m[0][2]= ac;
	m[0][3] = ad;
	//
	m[2][1] = m[1][2]= bc;
	m[1][3] = bd;
	//
	m[2][3] = cd;

	m[3][0] = 0.0;
   m[3][1] = 0.0;
   m[3][2] = 0.0;
	m[3][3] = 1.0;

   if (!invertMatrixd((double *)inverse, (double *)m))
      return false;

   // punkt ergibt sich durch die multiplikation mit [0 0 0 1]
   p.x = (float) inverse[0][3];
   p.y = (float) inverse[1][3];
   p.z = (float) inverse[2][3];

   // Inverse konnte gebildet werden
   return true;
}

} // end namespace


