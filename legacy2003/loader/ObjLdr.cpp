/**
 * Klasse zum Laden von '.obj' Modelldateien
 *
 * File Format
 */
#include "objldr.h"
#include <vector>
#include <stdio.h>
#include <fstream>

////////////////////////////////////////////////////////////////////////////
// reads a face index
////////////////////////////////////////////////////////////////////////////
void ObjLdr::read_face3(char *line, std::vector<int> *faces)
{
   // face indices
   int indices[10];
   // skip first token
   char *token = strtok(line, " \t");
   int count = 0;
   // read data
   while ((token = strtok(NULL, " \t")) && count < 10) {
      // get next index
      indices[count++] = (token!=NULL ? atoi(token) : 1) - 1;
   }
   // check if we've got enough indices
   if (count < 3) return;
   // build the first triangle
   int f1 = indices[0];
   int f2 = indices[1];
   int f3 = indices[2];
   //
   faces->push_back(f1);
   faces->push_back(f2);
   faces->push_back(f3);

   for (int i=3; i<count; i++) {
      f2 = f3;
      f3 = indices[i];
      //
      faces->push_back(f1);
      faces->push_back(f2);
      faces->push_back(f3);
   }
}

////////////////////////////////////////////////////////////////////////////
// reads a vertex
////////////////////////////////////////////////////////////////////////////
void ObjLdr::read_vertex3(char *line, std::vector<float> *coords)
{
   // skip first token
   char *token = strtok(line, " \t");
   // read data
   for (int i=0; i<3; i++) {
      // get next token
      token = strtok(NULL, " \t");
      float number = token!=NULL ? (float) atof(token) : 0.0f;
      coords->push_back(number);
   }
}

////////////////////////////////////////////////////////////////////////////
// static entry point for loading
////////////////////////////////////////////////////////////////////////////
bool ObjLdr::load(char *fileName, float **coord3, int **face3, int *numCoord, int *numFaces)
{
   std::vector<float> coords;
   std::vector<int> faces;
   char row[2048];
   int line = 0;
   //
   // create file input stream
   std::ifstream in(fileName);
   // check file open success
   if (!in)
      throw ParseException("Unable to open file", "InventorLdr::load");

   while (in.getline(row, 2048))
   {
      line++;
      switch (row[0]) {
         case 'v':
            read_vertex3(row, &coords);
            break;

         case 'f':
            read_face3(row, &faces);
            break;
      }
   }
   // copy coordinates
   *coord3 = new float[3 * coords.size()];
   memcpy(*coord3, &coords[0], sizeof(float) * coords.size());
   *numCoord = coords.size() / 3;

   // copy faces
   *face3 = new int[3 * faces.size()];
   memcpy(*face3, &faces[0], sizeof(int) * faces.size());
   *numFaces = faces.size() / 3;
   return false;
}




