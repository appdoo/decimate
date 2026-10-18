/**
 * Algorithm3D
 *
 * Hilfsfunktionen unter anderem für die Berechnung
 * von Normalenvektoren, Auffunden aller Kanten in
 * einem Dreiecksmesh etc...
 *
 **/
#ifndef ALGORITHM3D_H
#define ALGORITHM3D_H
#include <windows.h>
#include <stdio.h>
#include <math.h>
#include "math3d.h"



class Algorithm3D
{
public:

   /***************************************************************************
    *
    * Berechnet die Normalenvektoren für ein gegebenes Mesh
    *
    * Der Normalenvektor wird für jede der gegebenen Flächen berechnet
    * und zu allen drei eckpunkten addiert.
    * Im zweiten Schritt werden die Normalen dann normiert
    *
    ***************************************************************************/
   static void createNormals(
      Vec3 **norm3,
      Vec3 *coord3,
      Face *face3,
      int numCoord,
      int numFaces
      )
   {
      // Alle normalen auf [0, 0, 0] initialisieren
      *norm3 = new Vec3[numCoord];
      memset(*norm3, 0, sizeof(Vec3) * numCoord);

      for (int i=0; i<numFaces; i++) 
      {
         Vec3 a = coord3[face3[i].b] - coord3[face3[i].a];
         Vec3 b = coord3[face3[i].c] - coord3[face3[i].a];

         Vec3 nrm = a ^ b;

         (*norm3)[face3[i].a] += nrm;
         (*norm3)[face3[i].b] += nrm;
         (*norm3)[face3[i].c] += nrm;
      }
      // normalen auf länge 1 bringen
      for (i=0; i<numCoord; i++)
      {
        (*norm3)[i].norm();
      }
   }



};



#endif