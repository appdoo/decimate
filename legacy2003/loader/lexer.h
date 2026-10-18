/**
 * md5scanner
 *
 * class for scanning md5 files
 *
 *
 */

#ifndef LEXER_H
#define LEXER_H

#include <fstream>
#include "stringbuffer.h"

#define DEFAULT_PARSE_MSG "An unknown display exception has occurred."
#define DEFAULT_PARSE_METHOD "Unknown"

///////////////////////////////////////////////////////////////////////////////
// A parser exception
///////////////////////////////////////////////////////////////////////////////
class ParseException
{
	// set a default string representation for the display exception
	char *msg;

	// by default, it is not known, which method raised this exception
	char *method;

public:

	/**
	 * construct a standard DisplayException
	 * this shouldn't be used - cause it contains
	 * nearly no information
	 *
	 * @author: Andor Greissl
	 * @creation: 29.09.2002
	 */
   ParseException() : msg(DEFAULT_PARSE_MSG), method(DEFAULT_PARSE_METHOD) {}

	
	/**
	 * construct a DisplayException
	 * with the supplied message text as an description
	 *
	 * @param msg   a text describing what has happened
	 *
	 * @author: Andor Greissl
	 * @creation: 29.09.2002
	 */
	ParseException(char *msg) : msg(msg), method(DEFAULT_PARSE_METHOD) {}


	/**
	 * construct a DisplayException
	 * with the supplied message text as an description
	 *
	 * @param msg      a text describing what has happened
	 * @param method   a text describing which method raised this exception
	 *
	 * @author: Andor Greissl
	 * @creation: 29.09.2002
	 */
	ParseException(char *msg, char *method) : msg(msg), method(method) {}


	/**
	 * returns the message text
	 *
	 * @return msg        a text describing what has happened
	 *
	 * @author: Andor Greissl
	 * @creation: 29.09.2002
	 */
	char *getMsg()
	{
		return msg;
	}


	/**
	 * returns the method text
	 *
	 * @param method   a text describing which method raised this exception
	 *
	 * @author: Andor Greissl
	 * @creation: 29.09.2002
	 */
	char *getMehtod()
	{
		return method;
	}
};



///////////////////////////////////////////////////////////////////////////////
// A very simple and basic lexer class
///////////////////////////////////////////////////////////////////////////////
class Lexer
{
private:
	// the actual lexem string
	StringBuffer lexem;

   // user defined lexem values
   int  numTokens;
   char **lexems;
   int  *tokens;

   // error message
   char errMsg[2048];

	// pointer to the actual fileinput stream
	std::ifstream *in;

	// constructor state
	bool ok;

public:
   // aktuelle Zeile
   int line;

   enum 
	{ 
		INT,
		FLOAT,	
		STRING,	
		BROPEN,	
		BRCLOSE, 
      PARENOPEN,
      PARENCLOSE,
      ANGLEOPEN,
      ANGLECLOSE,
      SQUAREOPEN,
      SQUARECLOSE,
      COMMA,
      LEXEM,
		UNDEFINED,
      SLASH,
      USER_TOKEN
	};



	///////////////////////////////////////////////////////////////////////////////
	// construct a Lexer
	///////////////////////////////////////////////////////////////////////////////
	Lexer(char *fileName);


	///////////////////////////////////////////////////////////////////////////////
	// destruct the Lexer
	///////////////////////////////////////////////////////////////////////////////
	~Lexer();


	///////////////////////////////////////////////////////////////////////////////
	// registers a lexem token specification
	///////////////////////////////////////////////////////////////////////////////
	void reg(char *lexem, int token);


	///////////////////////////////////////////////////////////////////////////////
	// return the next recognized token
	///////////////////////////////////////////////////////////////////////////////
	int lex();


	///////////////////////////////////////////////////////////////////////////////
	// check if MD5Scanner was successfully created
	///////////////////////////////////////////////////////////////////////////////
	bool operator !();

	
	///////////////////////////////////////////////////////////////////////////////
	// returns a pointer to the current lexem
	///////////////////////////////////////////////////////////////////////////////
	const char *getLexem();


	///////////////////////////////////////////////////////////////////////////////
	// returns the current lexem as an integer
	///////////////////////////////////////////////////////////////////////////////
   int getInt();


	///////////////////////////////////////////////////////////////////////////////
	// returns the current lexem as a float
	///////////////////////////////////////////////////////////////////////////////
   float getFloat();


	///////////////////////////////////////////////////////////////////////////////
	// checks if a string is a user defined lexem 
	///////////////////////////////////////////////////////////////////////////////
   int checkLexem(const char *lexem);


	///////////////////////////////////////////////////////////////////////////////
	// match a token
	///////////////////////////////////////////////////////////////////////////////
   void match(int token);


	///////////////////////////////////////////////////////////////////////////////
	// set the user defined tokens
	///////////////////////////////////////////////////////////////////////////////
   void setUserLexems(char **lexems, int *tokens, int numtoken);
};




#endif

