// ============================================================================
// Class Engine 
//
// Definition of a OpenGl Rendering Engine
//
// $Author: Andor $
// $Date: 17.12.03 2:05 $
// $Revision: 4 $
//
// ============================================================================

#ifndef ENGINE_H
#define ENGINE_H

// C++ Includes
#include <windows.h>
#include <gl\gl.h>
#include <gl\glu.h>
#include <gl\glaux.h>

// Local Includes
#include "../math/math3D.h"
#include "texture.h"
#include "normalizationcubemap.h"

class Engine
{
   // define an Input - Model
   int   inputModel;
   //
   // attributes for input - model
   int   mouseDownX;
   int   mouseDownY;
   float curquat[4];
   float lastquat[4];

   // current viewport
   GLint x;
   GLint y;
   GLint width;
   GLint height;
   //
   GLfloat angle;
   GLfloat frontClippingPlane;
   GLfloat backClippingPlane;

   // current states
   GLint shadeMode;
   GLint fillMode;
   GLint faceMode;
   GLint matrixMode;
   //
   bool  depthTest;
   GLint depthFunc;

   // Hintergrundfarbe und aktuelle
   // Zeichenfarbe
   rgbf background;
   rgbf color;

   // Cubemap
   NormalizationCubeMap *cubeMap;

   // Aktuelle verfügbare Daten
   Texture *bumpTex;
   Texture *diffTex;
   Texture *specTex;

   int      numPoints;

   Vec3    *points;
   Vec3    *normals;
   rgbf    *rgbvalues;
   Tex2    *textures;

   Vec3    *lightmap;

   Face    *faces;
   int      numFaces;

   Edge    *edges;
   int      numEdges;

   Vec3     lpos;

   GLuint pointLightFragmentShaderIDATI;

   int drawMode;

public:
   // some input - model - types
   enum {NO_INPUT, TRACKBALL_INPUT};
   enum {
      WIREFRAME,
      COLORED,
      TEXTURED,
      LINED,
      BUMP
   };
   
   // -------------------------------------------------------------------------
	// construct a standard engine
   // -------------------------------------------------------------------------
	Engine();

   // -------------------------------------------------------------------------
   // delete engine
   // -------------------------------------------------------------------------
   ~Engine();

   // -------------------------------------------------------------------------
	// applies the inputmodel to the transformation stack
   // -------------------------------------------------------------------------
   void applyInput();

   // -------------------------------------------------------------------------
	//
	// In the next section of code we do all of the setup for OpenGL. 
	// We set what color to clear the screen to, we turn on the depth buffer,
	// enable smooth shading, etc. This routine will not be called until the OpenGL 
	// Window has been created. 
	// This procedure returns a value but because our initialization isn't that complex 
	// we wont worry about the value for now.
	//
	// Select a Polygon Model   GL_POINT
	//                          GL_LINE
	//                          GL_FILL
	//
	//    and a Shade Model     GL_FLAT
	//                          GL_SMOOTH
   //
   // -------------------------------------------------------------------------
	int init();
   
   // -------------------------------------------------------------------------
   // create an ATI Fragment Shader for Dot3 Bump Mapping
   // -------------------------------------------------------------------------
   void initATIDot3Bump();

   // -------------------------------------------------------------------------
   // tell the engine, that a mousebutton was pressed
   // -------------------------------------------------------------------------
   void mousePressed(int x, int y);

   // -------------------------------------------------------------------------
   // tell the engine, that a mousebutton was released
   // -------------------------------------------------------------------------
   void mouseReleased(int x, int y);

   // -------------------------------------------------------------------------
   // tell the engine, that the mouse was moved
   // -------------------------------------------------------------------------
   void mouseMoved(int x, int y);

   // -------------------------------------------------------------------------
   // Setzt alle verfügbaren Daten des Meshs
   // -------------------------------------------------------------------------
   void set(int numFaces, Face *faces, int numEdges, Edge *edges, int numPoints, Vec3 *points, Vec3 *normals, rgbf *rgbvalues, Tex2 *textures);

   // -------------------------------------------------------------------------
   // Setzt alle Texturen
   // -------------------------------------------------------------------------
   void set(Texture *bumpTex, Texture *diffTex, Texture *specTex);

   // -------------------------------------------------------------------------
   // Zeichenmodus wählen
   // -------------------------------------------------------------------------
   void setDrawMode(int drawMode);

   // -------------------------------------------------------------------------
   // Zeichnen
   // -------------------------------------------------------------------------
   void draw();

   // -------------------------------------------------------------------------
	// Resize And Initialize The GL Window
	//
	// The job of the next section of code is to resize the OpenGL scene whenever the window 
	// (assuming you are using a Window rather than fullscreen mode) has been resized.
	// Even if you are not able to resize the window (for example, you're in fullscreen mode), 
	// this routine will still be called at least once when the program is first run to set up 
	// our perspective view. The OpenGL scene will be resized based on the width and height of the window
	// it's being displayed in.
	//
	//
	// o  glMatrixMode(GL_PROJECTION) indicates that the next 2 lines of code
	//    will affect the projection matrix. The perspective matrix is responsible
	//    for adding perspective to our scene. glLoadIdentity() is similar to a reset.
	//    It restores the selected matrix to it's original state. After glLoadIdentity()
	//    has been called we set up our perspective view for the scene.
	//
	// o  glMatrixMode(GL_MODELVIEW) indicates that any new transformations will affect
	//    the modelview matrix. The modelview matrix is where our object information is stored.
	//    Lastly we reset the modelview matrix. 
	//    
   // -------------------------------------------------------------------------
	GLvoid setBounds(GLsizei x, GLsizei y, GLsizei width, GLsizei height);

   // -------------------------------------------------------------------------
   // translate
   // -------------------------------------------------------------------------
   void translate(float x, float y, float z);

   // -------------------------------------------------------------------------
   // scale
   // -------------------------------------------------------------------------
   void scale(float x, float y, float z);

   // -------------------------------------------------------------------------
   // rotate
   // -------------------------------------------------------------------------
   void rotate(float angle, float x, float y, float z);

   void setLightPosition(float x, float y, float z);

   void setColor(float r, float g, float b);

   // -------------------------------------------------------------------------
   // set background color
   // -------------------------------------------------------------------------
   void setClearColor(float r, float g, float b);

   // -------------------------------------------------------------------------
   // set z-buffer depth functionality
   // -------------------------------------------------------------------------
   void setDepthFunc(GLint func);

   // -------------------------------------------------------------------------
	// enable z-buffer depth test
   // -------------------------------------------------------------------------
	void setEnableDepthTest(bool enable);

   // -------------------------------------------------------------------------
	// select an input model
   // -------------------------------------------------------------------------
   void setInputModel(int inputModel);

   // -------------------------------------------------------------------------
   // sets a fillmode for opengl
   // -------------------------------------------------------------------------
   void setPolygonMode(GLint face, GLint mode);

   // -------------------------------------------------------------------------
   // selects a Matrix Mode for opengl
   // -------------------------------------------------------------------------
   void setMatrixMode(GLint mode);

   // -------------------------------------------------------------------------
   // sets a shademode for opengl
   // -------------------------------------------------------------------------
   void setShadeMode(GLint mode);

   // -------------------------------------------------------------------------
	// print a text on position 
	// x, y  -  a standard bitmapped
	// font with 8 pixels width 
	// is used
   // -------------------------------------------------------------------------
	void print(int x, int y, char *str);

   // -------------------------------------------------------------------------
   // reset all the stuff
   // -------------------------------------------------------------------------
	void clear();

   // -------------------------------------------------------------------------
   // Here's Where We Do All The Drawing
   // -------------------------------------------------------------------------
	int render();

   // -------------------------------------------------------------------------
   // visualization of edges that lie on an object
   // -------------------------------------------------------------------------
   void drawLines();

   // -------------------------------------------------------------------------
   // visualization of texture coordinates
   // -------------------------------------------------------------------------
   void drawTexture();

   // -------------------------------------------------------------------------
   // Wireframe visualization on a filled object
   // -------------------------------------------------------------------------
   void drawWireframe();

   // -------------------------------------------------------------------------
   // draw textured object
   // -------------------------------------------------------------------------
   void drawTextured();

   // -------------------------------------------------------------------------
   // visualization of edges that lie on an object
   // -------------------------------------------------------------------------
   void drawHiddenLines(float *coord3, int numFaces, int *faces, int numEdges, int *edges);

   // -------------------------------------------------------------------------
   // draw hidden lines onto a shaded object
   // -------------------------------------------------------------------------
   void drawColored();

   int CreateVertexArrayObjectATI(void *src, size_t size, int count);
   void DeleteVertexArrayObjectATI(int g_uiObjectBuffer);
   void RenderObjectATI(int uiVertexObject, int uiNormalObject, int faceCount, int *faces);

   // -------------------------------------------------------------------------
   // ATI Dot3 Bump - Mapping
   // -------------------------------------------------------------------------
   void Engine::renderDot3BumpATI();
};

#endif