/**
 * lexer
 *
 * basic lexer object
 */
#include "lexer.h"
#include <windows.h>
#include <stdio.h>
#include <fstream>


///////////////////////////////////////////////////////////////////////////////
// construct a Lexer
///////////////////////////////////////////////////////////////////////////////
Lexer::Lexer(char *fileName)
{
   numTokens = 0;
   tokens    = NULL;
   lexems    = NULL;
   //
   line      = 1;

	in = new std::ifstream(fileName);
	//
	ok = (in!=NULL) && (*in);
}


///////////////////////////////////////////////////////////////////////////////
// destruct the Lexer
///////////////////////////////////////////////////////////////////////////////
Lexer::~Lexer()
{
	if (in) delete in;
}


///////////////////////////////////////////////////////////////////////////////
// return the next recognized token
///////////////////////////////////////////////////////////////////////////////
int Lexer::lex()
{
	int status = 0;
	char c;

   lexem.clear();

	while (in->get(c)) {
		// case insensitiv scanning
		c = toupper(c);
		lexem.append(c);

		switch (status) {
			// initial state, start scanning
			case 0:
				if (isalpha(c)) status = 3;
				else if (isdigit(c) || c=='-') status = 1;
				else switch(c) {
               case '<':
                  return ANGLEOPEN;
               case '>':
                  return ANGLECLOSE;
               case '(':
                  return PARENOPEN;
               case ')':
                  return PARENCLOSE;
               case '[':
                  return SQUAREOPEN;
               case ']':
                  return SQUARECLOSE;
					case '{':
						return BROPEN;
					case '}':
						return BRCLOSE;
               case ',':
                  return COMMA;
               case '/':
                  return SLASH;
					case '"':
						status = 4;
						lexem--;
						break;
					case '\n': line++;
					case ' ':
					case '\t':
						lexem--;
						break;
					default:
						return UNDEFINED;
				}
				break;

			// numeric stuff found
			case 1:
				if (c=='.') status = 2;
            else if (c=='E' || c=='e') status = 5;
				else if (!isdigit(c)) {
					in->putback(c);
               lexem--;
					return INT;
				}
				break;

			// float found
			case 2:
				if (!isdigit(c)) {
               if (c=='E' || c=='e') status = 5;
					else {
                  in->putback(c);
					   return FLOAT;
               }
				}
				break;

			// token found
			case 3:
				if (!isalnum(c)) {
					in->putback(c);
               lexem--;
					return checkLexem(*lexem);
				}
				break;

			// string found
			case 4:
				if (c!='"')    lexem.append(c);
				else           return STRING;
				break;

         // case 5, float with exponent found
         case 5:
            if (c=='+' || c=='-' || isdigit(c)) status = 6;
            else {
               in->putback(c);
               return FLOAT;
            }
            break;

         // case 6, read exponent
         case 6:
				if (!isdigit(c)) {
               in->putback(c);
   			   return FLOAT;
   			}
				break;
		}
	}
	return EOF;
}


///////////////////////////////////////////////////////////////////////////////
// check if Lexer was successfully created
///////////////////////////////////////////////////////////////////////////////
bool Lexer::operator !()
{
	return (in==NULL)  ||  (!(*in));
}


///////////////////////////////////////////////////////////////////////////////
// returns the current lexem
///////////////////////////////////////////////////////////////////////////////
const char *Lexer::getLexem()
{
	return *lexem;
}


///////////////////////////////////////////////////////////////////////////////
// checks if a string is a user defined lexem 
///////////////////////////////////////////////////////////////////////////////
int Lexer::checkLexem(const char *lexem)
{
   for (int i=0; i<numTokens; i++)
   {
     if (stricmp(lexems[i], lexem) == 0)
        return tokens[i];
   }
   return STRING;
}


///////////////////////////////////////////////////////////////////////////////
// returns the current lexem as an integer
///////////////////////////////////////////////////////////////////////////////
int Lexer::getInt()
{
   return atoi(*lexem);
}


////////////////////////////////////////////////////////////////////////////
// returns the current lexem as a float
///////////////////////////////////////////////////////////////////////////////
float Lexer::getFloat()
{
   return (float) atof(*lexem);
}


///////////////////////////////////////////////////////////////////////////////
// match a token
///////////////////////////////////////////////////////////////////////////////
void Lexer::match(int token)
{
  if (lex() != token)
  {
     sprintf(errMsg, "Fehler beim Einlesen der Datei! Zeile: %d, token '%s'", line, *lexem);
     throw ParseException(errMsg, "Lexer::match");
  }
}


///////////////////////////////////////////////////////////////////////////////
// set the user defined tokens
///////////////////////////////////////////////////////////////////////////////
void Lexer::setUserLexems(char **lexems, int *tokens, int numTokens)
{
   this->numTokens = numTokens;
   this->lexems    = lexems;
   this->tokens    = tokens;
}