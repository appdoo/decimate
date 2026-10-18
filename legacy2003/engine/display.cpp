/**
 * display.cpp
 *
 * functions that deal with the rendering context
 *
 **/

#include <windows.h>
//#include "..\Resource.h"
#include "display.h"
#include "displayexception.h"



// define class name / title for display
const char *className = "GlRenderingEngineWndClass";
const char *childClassName = "GlRenderingEngineChildWndClass";
const char *windowTitle = "OpenGl Testbed 0.1";

// our global title
const char *dialogTitle = windowTitle;



///////////////////////////////////////////////////////////////////////////////
// construct a Standard Display
//
// @author: Andor Greissl
// @creation: 29.09.2002
///////////////////////////////////////////////////////////////////////////////
Display::Display()
{
	// create a standard window
	Display(DefWndProc, NULL, 0, 0, 640, 480, 16, false);
}


/**
 * construct a Display with given Dimensions
 *
 * @param width       requested display width
 * @param height      requested display height
 * @param depth       requested depth
 * @param fullscreen  flag, if Display should be created in fullscreen mode
 *
 * @author: Andor Greissl
 * @creation: 29.09.2002
 */
Display::Display(WNDPROC WndProc, HINSTANCE hInstance, int x, int y, int width, int height, int depth, bool fullscreen)
{
	// copy instantiation parameters
	this->WndProc    = WndProc;
	this->width      = width;
	this->height     = height;
	this->depth      = depth;
	this->fullscreen = fullscreen;

	// copy window dimensions
	WindowRect.left   = (long) (x);
	WindowRect.right  = (long) (x + width);
	WindowRect.top    = (long) (y);
	WindowRect.bottom = (long) (y + height);

	// set flags for windows
	// flags are for standard window (no fullscreen mode)
	dwExStyle = WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;
	dwStyle   = WS_OVERLAPPEDWINDOW;

	// clear all Handles and References
	hRC = NULL;
	hDC = NULL;
	hWnd = NULL;
	this->hInstance = hInstance;

  	// grab an Instance for our window
   if (!hInstance)
   	if(!(hInstance = GetModuleHandle(NULL)))
	   	throw DisplayException("No hInstance defined for Display.", "Display::Constructor");

	// register a windowclass for the Display
	registerWindowClasses();

   // check if we have a WndProc
	if(!WndProc)
		throw DisplayException("No WndProc defined for Display.", "Display::Constructor");
}


/**
 * deletes the Display
 *
 * @author: Andor Greissl
 * @creation: 29.09.2002
 */
Display::~Display()
{
	killDisplay();
}


/////////////////////////////////////////////////////////////////////////////// 
// returns the screen width
/////////////////////////////////////////////////////////////////////////////// 
int Display::getScreenWidth ()
{
//	if (GetSystemMetrics(SM_CMONITORS) < 2) 
	{
		return GetSystemMetrics (SM_CXSCREEN);
	}
//	else
	{
//		return GetSystemMetrics (SM_CXVIRTUALSCREEN);
	}
}


/////////////////////////////////////////////////////////////////////////////// 
// returns the screen height
/////////////////////////////////////////////////////////////////////////////// 
int Display::getScreenHeight ()
{
//	if (GetSystemMetrics(SM_CMONITORS) < 2) 
   {
      return GetSystemMetrics (SM_CYSCREEN);
   }
//	else
   {
//		return GetSystemMetrics (SM_CYVIRTUALSCREEN);
   }
}


/**
 * initializes win32 structures, registers
 * the window class, sets fullscreen mode
 * and shows the Display
 *
 * @exception DisplayException
 *
 * @author: Andor Greissl
 * @creation: 29.09.2002
 */
void Display::create()
{
	// change display resolution
	changeResolution();

	// Adjust Window To True Requested Size
	AdjustWindowRectEx(&WindowRect, dwStyle, FALSE, dwExStyle);

	// create a win32 window
	hWnd = CreateWindowEx(
		dwExStyle,                          // Extended Style For The Window
		className,                          // Class Name
		windowTitle,                        // Window Title
		WS_CLIPSIBLINGS |                   // Required Window Style
		WS_CLIPCHILDREN |                   // Required Window Style
		dwStyle,                            // Selected Window Style
		WindowRect.left,                    // Window Position
      WindowRect.top,                     // Window Position
		WindowRect.right - WindowRect.left,	// Calculate Adjusted Window Width
		WindowRect.bottom - WindowRect.top,	// Calculate Adjusted Window Height
		NULL,                               // No Parent Window
		NULL,                               // No Menu
		hInstance,                          // Instance
		NULL
	);


	// check against window creation success
	if (!hWnd)
		throw DisplayException("Unable to create window.", "Display::show");

	// get the device context
	if (!(hDC = GetDC(hWnd)))
		throw DisplayException("gdi::GetDC  Get Device Context failed.", "Display::show");

	// set the pixel format
	setPixelFormat(hDC);

	// create our rendering context
	createOpenGlContext(hDC);

	// show the window
	ShowWindow(hWnd,SW_SHOW);
   UpdateWindow(hWnd);
   SetForegroundWindow(hWnd);
	SetFocus(hWnd);
}


///////////////////////////////////////////////////////////////////////////////
// creates a child - window as child of a 
// parent window
///////////////////////////////////////////////////////////////////////////////
HWND Display::createChild(HWND parent, long id)
{
   unregisterWindowClass(hInstance, childClassName);
     
   WNDCLASS     wndclass ;
     
   wndclass.style         = CS_OWNDC;
   wndclass.lpfnWndProc   = WndProc ;
   wndclass.cbClsExtra    = 0 ;
   wndclass.hInstance     = (HINSTANCE) GetWindowLong (parent, GWL_HINSTANCE) ;
   wndclass.hIcon         = NULL ;
   wndclass.hCursor       = NULL ;
   wndclass.hbrBackground = (HBRUSH) GetStockObject (BLACK_BRUSH) ;
   wndclass.lpszMenuName  = NULL ;
   wndclass.lpszClassName = childClassName ;
   wndclass.cbWndExtra    = sizeof (long) ;

   RegisterClass (&wndclass) ;

   hWnd = CreateWindow(
      childClassName, 
      NULL,
      WS_CHILDWINDOW | WS_VISIBLE,
      0,
      0,
      0,
      0,
      parent, 
      (HMENU) id,
      (HINSTANCE) GetWindowLong (parent, GWL_HINSTANCE),
      NULL
   );

  	// check against window creation success
	if (!hWnd)
		throw DisplayException("Unable to create window.", "Display::createChild");

	// get the device context
	if (!(hDC = GetDC(hWnd)))
		throw DisplayException("gdi::GetDC  Get Device Context failed.", "Display::show");

	// set the pixel format
	setPixelFormat(hDC);

	// create our rendering context
	createOpenGlContext(hDC);

   return hWnd;
}



/**
 * swap the buffers when using double buffered painting
 * @exception DisplayException
 *
 * @author: Andor Greissl
 * @creation: 29.09.2002
 */
void Display::swapBuffers()
{
	SwapBuffers(hDC);	
}



/**
 * method sets a Pixelformat for our Display
 *
 * @exception DisplayException
 *
 * @author: Andor Greissl
 * @creation: 29.09.2002
 */
void Display::setPixelFormat(HDC hDC)
{
	// win32 Pixel format descriptor -> what we wan't to have
	static	PIXELFORMATDESCRIPTOR pfd=
	{
		sizeof(PIXELFORMATDESCRIPTOR),    // Size Of This Pixel Format Descriptor
		1,                                // Version Number
		PFD_DRAW_TO_WINDOW |              // Format Must Support Window
		PFD_SUPPORT_OPENGL |              // Format Must Support OpenGL
		PFD_DOUBLEBUFFER,                 // Must Support Double Buffering
		PFD_TYPE_RGBA,                    // Request An RGBA Format
		(unsigned char) depth,            // Select Our Color Depth
		0, 0, 0, 0, 0, 0,                 // Color Bits Ignored
		8,                                // No Alpha Buffer
		0,                                // Shift Bit Ignored
		0,                                // No Accumulation Buffer
		0, 0, 0, 0,                       // Accumulation Bits Ignored
		24,                               // 24 Bit Z-Buffer (Depth Buffer)
		1,                                // No Stencil Buffer
		0,                                // No Auxiliary Buffer
		PFD_MAIN_PLANE,                   // Main Drawing Layer
		0,                                // Reserved
		0, 0, 0                           // Layer Masks Ignored
	};

	// find a pixel format that matches our request
	if (!(PixelFormat = ChoosePixelFormat(hDC, &pfd)))
		throw DisplayException("No matching PixelFormat found.", "Display::setPixelFormat");


	// set the pixelformat in our device context
	if(!SetPixelFormat(hDC, PixelFormat, &pfd))
		throw DisplayException("Set PixelFormat failed.","Display::setPixelFormat");
}


/**
 * create our rendering context for the display
 *
 * @exception DisplayException
 *
 * @author: Andor Greissl
 * @creation: 29.09.2002
 */
void Display::createOpenGlContext(HDC hDC)
{
	// create a rendering context
	if (!(hRC = wglCreateContext(hDC)))
		throw DisplayException("wglCreateContext failed.","Display::createOpenGlContext");

	// make the context the current one
	if(!wglMakeCurrent(hDC, hRC))
		throw DisplayException("wglMakeCurrent failed.","Display::createOpenGlContext");
}


/**
 * method tries to change the display settings
 * to class attributes width, height, depth
 *
 * @return success
 *
 * @author: Andor Greissl
 * @creation: 29.09.2002
 */
void Display::changeResolution()
{
	// do nothing if fullscreen is not used
	if (!fullscreen)
		return;

	// set the currentstate (no fullscreen)
	fullscreen = false;

	// create Device Mode and clear it
	DEVMODE dmScreenSettings;
	memset(&dmScreenSettings, 0, sizeof(dmScreenSettings));

	dmScreenSettings.dmSize = sizeof(dmScreenSettings);
	dmScreenSettings.dmPelsWidth	= width;
	dmScreenSettings.dmPelsHeight	= height;
	dmScreenSettings.dmBitsPerPel	= depth;
	dmScreenSettings.dmFields     = DM_BITSPERPEL|DM_PELSWIDTH|DM_PELSHEIGHT;

	// Try To Set Selected Mode And Get Results.  NOTE: CDS_FULLSCREEN Gets Rid Of Start Bar.
	if (ChangeDisplaySettings(&dmScreenSettings, CDS_FULLSCREEN) != DISP_CHANGE_SUCCESSFUL)
		throw DisplayException("ChangeDisplaySettings failed.", "Display::ChangeResolution");

	// make some settings for fullscreen mode
	dwExStyle = WS_EX_APPWINDOW;
	dwStyle = WS_POPUP;
	ShowCursor(FALSE);

	// set current state to fullscreen
	fullscreen = true;
}


///////////////////////////////////////////////////////////////////////////////
// unregisters a windowclass
///////////////////////////////////////////////////////////////////////////////
void Display::unregisterWindowClass(HINSTANCE hInstance, const char *className)
{
   WNDCLASSEX wc;

   // check if the window class is already registered
	if (GetClassInfoEx (hInstance, className, &wc))	{
		UnregisterClass (className, hInstance);
	}
}


///////////////////////////////////////////////////////////////////////////////
// registers a windowclass for the display
//
// @exception DisplayException
//
// @author: Andor Greissl
// @creation: 29.09.2002
//////////////////////////////////////////////////////////////////////////////
void Display::registerWindowClasses()
{
   // check if window classes are already registered
   unregisterWindowClass(hInstance, className);
   unregisterWindowClass(hInstance, childClassName);

   // Windows Class Structure
   WNDCLASSEX wc;
   //
	wc.cbSize         = sizeof(WNDCLASSEX); 
	wc.style          = 0;  // redraw on move, and own DC for Window
	wc.lpfnWndProc    = WndProc;                             // WndProc Handles Messages
	wc.cbClsExtra     = 0;
	wc.cbWndExtra     = 0;
	wc.hInstance      = hInstance;
	wc.hIcon          = LoadIcon(NULL, IDI_WINLOGO);
	wc.hCursor        = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground  = (HBRUSH) COLOR_GRAYTEXT;
   wc.lpszMenuName   = NULL;
	wc.lpszClassName  = className;
	wc.hIconSm        = NULL;

	// try to register the window class
	if(!RegisterClassEx(&wc))
	  throw DisplayException("Error registering window class.", "Display::registerWindowClasses");

   // register a child window class
   wc.lpfnWndProc    = WndProc;
   wc.cbWndExtra     = sizeof (long);
   wc.hIcon          = NULL;
   wc.lpszClassName  = childClassName;
   wc.lpszMenuName   = NULL;
	wc.hIconSm        = NULL;

	// try to register the child window class
	if(!RegisterClassEx(&wc))
	  throw DisplayException("Error registering child window class.", "Display::registerWindowClasses");
}


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
void Display::killDisplay()
{
	// do we have to reset from fullscreen mode ?
	if (fullscreen) {
		ChangeDisplaySettings(NULL, 0);
		ShowCursor(TRUE);
	}

	// make the rendering context not current 
	if (hRC  &&  !wglMakeCurrent(NULL, NULL))
		showMessageDialog("Reset of the Rendering Context failed.");

	// delete the rendering context 
	if (hRC  &&  !wglDeleteContext(hRC))
		showMessageDialog("Couldn't delete the Rendering Context.");

	// release the device context
	if (hDC && !ReleaseDC(hWnd, hDC))
		showMessageDialog("Release Device Context Failed.");

	// destroy the window
	if (hWnd && !DestroyWindow(hWnd))
		showMessageDialog("Could Not Release hWnd.");

	// unregister the window classes
	unregisterWindowClass(hInstance, className);
	unregisterWindowClass(hInstance, childClassName);

	// reset values
	fullscreen = false;
	hRC = NULL;
	hDC = NULL;
	hWnd = NULL;
	hInstance = NULL;
}



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
LRESULT CALLBACK Display::DefWndProc(HWND	hWnd, UINT uMsg, WPARAM	wParam, LPARAM lParam)
{
	// dummy implementation - just sends a quit message
	switch (uMsg)
	{
	   case WM_COMMAND:
      case WM_ACTIVATE:
      case WM_SYSCOMMAND:
      case WM_LBUTTONDOWN:
      case WM_RBUTTONDOWN:
      case WM_MBUTTONDOWN:
      case WM_MOUSEMOVE:
      case WM_LBUTTONUP:
      case WM_RBUTTONUP:
      case WM_MBUTTONUP:
		case WM_KEYDOWN:	
		case WM_KEYUP:
		case WM_SIZE:
         break;
		case WM_CLOSE:
			PostQuitMessage(0);
			return 0;
	}

	// Pass All Unhandled Messages To DefWindowProc
	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}


/**
 * shows up a message dialog
 *
 * @param msg  the message text that should be displayed
 *
 * @author Andor Greissl
 * @creation: 27.09.2002
 */

void Display::showMessageDialog(char *msg)
{
  MessageBox(NULL, msg, dialogTitle, MB_OK|MB_ICONEXCLAMATION);
}



/**
 * shows up a yes / no - dialog
 *
 * @param msg  the message text that should be displayed
 * @return true / false
 *
 * @author Andor Greissl
 * @creation: 27.09.2002
 */

bool Display::showYesNoDialog(char *msg)
{
	return MessageBox(NULL, msg, dialogTitle, MB_YESNO|MB_ICONEXCLAMATION)==IDYES;
}





