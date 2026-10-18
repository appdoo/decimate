/**
 * Class DisplayException
 *
 * Class for casting and receiving exceptions that
 * are throen by the display class
 *
 *
 */

#ifndef DISPLAY_EXCEPTION_H
#define DISPLAY_EXCEPTION_H

#define DEFAULT_MSG "An unknown display exception has occurred."
#define DEFAULT_METHOD "Unknown"

class DisplayException
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
	DisplayException() : msg(DEFAULT_MSG), method(DEFAULT_METHOD) {}

	
	/**
	 * construct a DisplayException
	 * with the supplied message text as an description
	 *
	 * @param msg   a text describing what has happened
	 *
	 * @author: Andor Greissl
	 * @creation: 29.09.2002
	 */
	DisplayException(char *msg) : msg(msg), method(DEFAULT_METHOD) {}


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
	DisplayException(char *msg, char *method) : msg(msg), method(method) {}


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





#endif







