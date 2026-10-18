/**
 * class Display
 *
 * a structure to hold the rendering context - that is
 * for e.g. a window handle for a windows - rendering context
 *
 **/

#ifndef ENGINE_DISPLAY_H
#define ENGINE_DISPLAY_H

// include all the stuff we need
#include <windows.h>
#include <gl\glu.h>
#include "displayexception.h"



class Display
{
private:
   bool       fullscreen;     // flag if we are in fullscreen mode

   int        x;              // left position
   int        y;              // top position
   int        width;          // width of our window
   int        height;         // height of our window

   int        depth;          // rendering color - depth

   HGLRC      hRC;            // Permanent Rendering Context
   HDC        hDC;            // Private GDI Device Context
   HWND       hWnd;           // Holds Our Window Handle
   HINSTANCE  hInstance;      // Holds The Instance Of The Application

   WNDPROC    WndProc;        // the WNDPROC for this Display

   DWORD      dwExStyle;      // Window Extended Style
   DWORD      dwStyle;        // Window Style

   RECT       WindowRect;     // Grabs Rectangle Upper Left / Lower Right Values

   GLuint     PixelFormat;  // Holds The Results After Searching For A Match


public:
   ///////////////////////////////////////////////////////////////////////////////
	// construct a Standard Display
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	Display();


   ///////////////////////////////////////////////////////////////////////////////
	// construct a Display with given Dimensions
	//
	// @param width       requested display width
	// @param height      requested display height
	// @param depth       requested depth
	// @param fullscreen  flag, if Display should be created in fullscreen mode
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
   Display(WNDPROC WndProc, HINSTANCE hInstance, int x, int y, int width, int height, int depth, bool fullscreen);


   ///////////////////////////////////////////////////////////////////////////////
	// deletes the Display
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	~Display();


   /////////////////////////////////////////////////////////////////////////////// 
   // returns the screen width
   /////////////////////////////////////////////////////////////////////////////// 
   int getScreenWidth ();


   /////////////////////////////////////////////////////////////////////////////// 
   // returns the screen height
   /////////////////////////////////////////////////////////////////////////////// 
   int getScreenHeight ();


   ///////////////////////////////////////////////////////////////////////////////
	// initializesd win32 structures, registers
	// the window class, sets fullscreen mode
	// and shows the Display
	//
	// @exception DisplayException
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	void create();


   ///////////////////////////////////////////////////////////////////////////////
   // creates a child - window as child of a 
   // parent window
   ///////////////////////////////////////////////////////////////////////////////
   HWND createChild(HWND parent, long id);


   ///////////////////////////////////////////////////////////////////////////////
	// swap the buffers when using double buffered painting
	// @exception DisplayException
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	void swapBuffers();


private:
   ///////////////////////////////////////////////////////////////////////////////
	// registers a window - class for the display
	//
	// @exception  DisplayException
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	void registerWindowClasses();


   ///////////////////////////////////////////////////////////////////////////////
	// method sets a Pixelformat for our Display
	//
	// @exception DisplayException
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	void Display::setPixelFormat(HDC hDC);


   ///////////////////////////////////////////////////////////////////////////////
	// create our rendering context for the display
	//
	// @exception DisplayException
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	void createOpenGlContext(HDC hDC);


   ///////////////////////////////////////////////////////////////////////////////
	// method tries to change the display settings
	// to class attributes width, height, depth
	//
	// @exception DisplayException
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	void changeResolution();


   ///////////////////////////////////////////////////////////////////////////////
	// method tries to reset everything
	// closeing the window
	// unregistering the window class
	//
	// @return success
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	void killDisplay();


   ///////////////////////////////////////////////////////////////////////////////
   // unregisters a windowclass
   ///////////////////////////////////////////////////////////////////////////////
   void unregisterWindowClass(HINSTANCE hInstance, const char *className);


   ///////////////////////////////////////////////////////////////////////////////
   // a Dummy WNDPROC implementaion - used for standard constructor
	// 
	// This is where all the window messages are dealt with. 
	// When we registred the Window Class we told it to jump to this 
	// section of code to deal with window messages.
	//
	// @param hWnd            Handle For This Window
	// @param uMsg            Message For This Window
	// @param wParam          Additional Message Information
	// @param lParam          Additional Message Information
	//
	// @author: Andor Greissl
	// @creation: 29.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	static LRESULT CALLBACK DefWndProc(HWND	hWnd,	UINT uMsg, WPARAM	wParam, LPARAM	lParam);


public:
   ///////////////////////////////////////////////////////////////////////////////
	// shows up a message dialog
	//
	// @param msg  the message text that should be displayed
	//
	// @author Andor Greissl
	// @creation: 27.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	static void showMessageDialog(char *msg);


   ///////////////////////////////////////////////////////////////////////////////
	// shows up a yes / no - dialog
	//
	// @param msg  the message text that should be displayed
	// @return true / false
	//
	// @author Andor Greissl
	// @creation: 27.09.2002
   ///////////////////////////////////////////////////////////////////////////////
	static bool showYesNoDialog(char *msg);
};

#endif