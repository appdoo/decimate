/**
 * Klasse zum Laden von '.iv' Modelldateien
 *
 * Inventor - File Format
 */
#include "inventorldr.h"
#include <vector>


char *InventorLdr::lookup[] =
{
	"#Inventor",
	"SEPARATOR",
	"COORDINATE3",
	"POINT",
	"MATERIAL",
	"DIFFUSECOLOR",
	"MATERIALBINDING",
	"VALUE",
	"PER_VERTEX_INDEXED",
	"INDEXEDFACESET",
	"COORDINDEX",
	"MATERIALINDEX",
	"ASCII",
};

int InventorLdr::tokens[] =
{
	INVENTOR,
	SEPARATOR,
	COORDINATE3,
	POINT,
	MATERIAL,
	DIFFUSECOLOR,
	MATERIALBINDING,
	VALUE,
	PER_VERTEX_INDEXED,
	INDEXEDFACESET,
	COORDINDEX,
	MATERIALINDEX,
	ASCII,
};


////////////////////////////////////////////////////////////////////////////
// static entry point for loading
////////////////////////////////////////////////////////////////////////////
bool InventorLdr::load(char *fileName, float **coord3, int **faces, int *numCoord, int *numFaces)
{
   Lexer lexer(fileName);

   if (!lexer)
      throw ParseException("Unable to open file", "InventorLdr::load");

   int token = 0;

   lexer.setUserLexems(lookup, tokens, sizeof(tokens) / sizeof(int));

   while ((token = lexer.lex()) != EOF)
   {
     if (token == SEPARATOR) {
        readSeparator(lexer, coord3, faces, numCoord, numFaces);
     }
   }
   return false;
}


////////////////////////////////////////////////////////////////////////////
// read single mesh
////////////////////////////////////////////////////////////////////////////
void InventorLdr::readSeparator(Lexer &lexer, float **coord3, int **faces, int *numCoord, int *numFaces)
{
   int token;

   lexer.match(Lexer::BROPEN);

   while ((token = lexer.lex()) != EOF) {
      // read a part of the Inventor model
      switch (token) {
         case COORDINATE3:
            readCoordinate3(lexer, coord3, numCoord);
            break;

         case MATERIAL:
            readMaterial(lexer);
            break;

         case MATERIALBINDING:
            readMaterialBinding(lexer);
            break;

         case INDEXEDFACESET:
            readIndexedFaceSet(lexer, faces, numFaces);
            break;

         case Lexer::BRCLOSE:
            return;
      }
   }
   throw ParseException("Unexpected end of file", "InventorLdr::readSeparator");
}


////////////////////////////////////////////////////////////////////////////
// read COORDINATES
////////////////////////////////////////////////////////////////////////////
void InventorLdr::readCoordinate3(Lexer &lexer, float **coords, int *numCoord)
{
   int token;
   std::vector<float> coord3;

   lexer.match(Lexer::BROPEN);
   lexer.match(POINT);
   lexer.match(Lexer::SQUAREOPEN);

   while ((token = lexer.lex()) != EOF) {
      // read a part of the Inventor model
      switch (token) {
         case Lexer::FLOAT:
            coord3.push_back(lexer.getFloat());

            lexer.match(Lexer::FLOAT); 
            coord3.push_back(lexer.getFloat());
            
            lexer.match(Lexer::FLOAT); 
            coord3.push_back(lexer.getFloat());
            break;

         case Lexer::COMMA:
            break;

         case Lexer::SQUARECLOSE:
            lexer.match(Lexer::BRCLOSE);

            *coords = new float[3 * coord3.size()];
            memcpy(*coords, &coord3[0], sizeof(float) * coord3.size());
            *numCoord = coord3.size() / 3;
            return;
      }
   }
   throw ParseException("Unexpected end of file", "InventorLdr::readCoordinate3");
}


////////////////////////////////////////////////////////////////////////////
// read material
////////////////////////////////////////////////////////////////////////////
void InventorLdr::readMaterial(Lexer &lexer)
{
   int token;

   lexer.match(Lexer::BROPEN);

   while ((token = lexer.lex()) != EOF) {
      // read a part of the Inventor model
      switch (token) {

         case Lexer::BRCLOSE:
            return;
      }
   }
   throw ParseException("Unexpected end of file", "InventorLdr::readMaterial");
}


////////////////////////////////////////////////////////////////////////////
// read material bindings
////////////////////////////////////////////////////////////////////////////
void InventorLdr::readMaterialBinding(Lexer &lexer)
{
   int token;

   lexer.match(Lexer::BROPEN);

   while ((token = lexer.lex()) != EOF) {
      // read a part of the Inventor model
      switch (token) {

         case Lexer::BRCLOSE:
            return;
      }
   }
   throw ParseException("Unexpected end of file", "InventorLdr::readMaterialBinding");
}


////////////////////////////////////////////////////////////////////////////
// read face set
////////////////////////////////////////////////////////////////////////////
void InventorLdr::readIndexedFaceSet(Lexer &lexer, int **faces, int *numFaces)
{
   int token;
   std::vector<int> face3;

   lexer.match(Lexer::BROPEN);
   lexer.match(COORDINDEX);
   lexer.match(Lexer::SQUAREOPEN);

   while ((token = lexer.lex()) != EOF) {
      // read a part of the Inventor model
      switch (token) {
         case Lexer::INT:
            face3.push_back(lexer.getInt());
            lexer.match(Lexer::COMMA);
            lexer.match(Lexer::INT); face3.push_back(lexer.getInt());
            lexer.match(Lexer::COMMA);
            lexer.match(Lexer::INT); face3.push_back(lexer.getInt());
            lexer.match(Lexer::COMMA);
            lexer.match(Lexer::INT);
            lexer.match(Lexer::COMMA);
            break;

         case Lexer::COMMA:
            break;

         case Lexer::SQUARECLOSE:
            *faces = new int[3 * face3.size()];
            memcpy(*faces, &face3[0], sizeof(int) * face3.size());
            *numFaces = face3.size() / 3;
            return;
      }
   }
   throw ParseException("Unexpected end of file", "InventorLdr::readIndexedFaceSet");
}



