/**
 * Klasse zum Laden von '.iv' Modelldateien
 *
 * Inventor - File Format
 */
#ifndef OBJLDR_H
#define OBJLDR_H

#include "lexer.h"
#include <vector>

class ObjLdr
{
   static char *lookup[];
   static int tokens[];
   int line;
public:
   enum
   {
      OBJ_G = Lexer::USER_TOKEN,
	   SHARP,
	   VERTICES,
	   OBJ_V,
	   SURFACE,
	   POLYGONS,
	   OBJ_F,
      OBJ_SLASH,
   };

   ////////////////////////////////////////////////////////////////////////////
   // static entry point for loading
   ////////////////////////////////////////////////////////////////////////////
   static bool load(char *fileName, float **coord3, int **faces, int *numCoord, int *numFaces);

   ////////////////////////////////////////////////////////////////////////////
   // reads a face index
   ////////////////////////////////////////////////////////////////////////////
   static void read_face3(char *line, std::vector<int> *faces);

   ////////////////////////////////////////////////////////////////////////////
   // reads a vertex
   ////////////////////////////////////////////////////////////////////////////
   static void read_vertex3(char *line, std::vector<float> *coords);
};


#endif

