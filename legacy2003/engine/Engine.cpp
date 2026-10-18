// ============================================================================
//
// Engine implementation
//
// ============================================================================
#include "engine.h"
#include "..\math\trackball.h"
#include "ext_gl.h"

   GLfloat mat_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
   GLfloat mat_shininess[] = { 100.0f };
   GLfloat white_light[] = { 0.8f, 0.8f, 0.8f, 1.0f };
   GLfloat lmodel_ambient[] = { 0.8f, 0.8f, 0.8f, 1.0f };
   GLfloat lightposition[] = { 1.0f, 1.0f, 1.0f, 0.0f };


// ----------------------------------------------------------------------------
// construct a standard engine
// ----------------------------------------------------------------------------
Engine::Engine()
{
  init();
}


// ----------------------------------------------------------------------------
// delete engine
// ----------------------------------------------------------------------------
Engine::~Engine()
{
   if (cubeMap)  delete   cubeMap;
   if (lightmap) delete[] lightmap;

   if (pointLightFragmentShaderIDATI > 0) glDeleteFragmentShaderATI(pointLightFragmentShaderIDATI);
}

// ----------------------------------------------------------------------------
// applies the inputmodel to the transformation stack
// ----------------------------------------------------------------------------
void Engine::applyInput()
{
   if (inputModel == TRACKBALL_INPUT) {
      float rotate[4][4];
      build_rotmatrix(rotate, curquat);
      glMultMatrixf(&rotate[0][0]);
   }
}


// ----------------------------------------------------------------------------
// tell the engine, that a mousebutton was pressed
// ----------------------------------------------------------------------------
void Engine::mousePressed(int x, int y)
{
   if (inputModel == TRACKBALL_INPUT) {
      trackball(lastquat, 0, 0, 0, 0);
      mouseDownX = x;
      mouseDownY = y;
   }
}


// ----------------------------------------------------------------------------
// tell the engine, that a mousebutton was released
// ----------------------------------------------------------------------------
void Engine::mouseReleased(int x, int y)
{
}


// ----------------------------------------------------------------------------
// tell the engine, that the mouse was moved
// ----------------------------------------------------------------------------
void Engine::mouseMoved(int x, int y)
{
   if (inputModel == TRACKBALL_INPUT ) {
      trackball(
         lastquat,
         (2.0f * width/2.0f - width) / width,
         (height - 2.0f * height/2.0f) / height,
         (2.0f * (x+(-mouseDownX+width/2.0f)) - width) / width,
         (height - 2.0f * (y+(-mouseDownY+height/2.0f))) / height
      );
      mouseDownX = x;
      mouseDownY = y;
      add_quats(lastquat, curquat, curquat);
   }
}


// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------
GLvoid Engine::setBounds(GLsizei x, GLsizei y, GLsizei width, GLsizei height)
{
   // height has always to be > 0 (Aspect Ratio)
	if (height==0)  height=1;

   // copy reference values
   this->x = x;
   this->y = y;
   this->width  = width;
   this->height = height;
   
   // set The Current Viewport, (window constants)
   glViewport(x, y, width, height);

   setMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Calculate The Aspect Ratio Of The Window
   gluPerspective(angle ,(GLfloat)width/(GLfloat)height , frontClippingPlane, backClippingPlane);

   setMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


// ----------------------------------------------------------------------------
// set z-buffer depth functionality
// ----------------------------------------------------------------------------
void Engine::setDepthFunc(GLint func)
{
   if (depthFunc!=func) {
      depthFunc = func;
      glDepthFunc(func);
   }
}

   
// ----------------------------------------------------------------------------
// enable z-buffer depth test
// ----------------------------------------------------------------------------
void Engine::setEnableDepthTest(bool enable)
{
   // check current state
   if (depthTest==enable) 
      return;

   // switch state
   depthTest = enable;
   if (enable) glEnable(GL_DEPTH_TEST);
   else        glDisable(GL_DEPTH_TEST);
}

   
// ----------------------------------------------------------------------------
// select an input model
// ----------------------------------------------------------------------------
void Engine::setInputModel(int inputModel)
{
   this->inputModel = inputModel;

   if (inputModel == TRACKBALL_INPUT) {
      trackball(curquat, -0.5, 0.0, 0.0, 0.0);
   }
}


// ----------------------------------------------------------------------------
// sets a fillmode for opengl
// ----------------------------------------------------------------------------
void Engine::setPolygonMode(GLint face, GLint mode)
{
   if (mode!=fillMode || face!=faceMode) {
      fillMode = mode;
      faceMode = face;
      glPolygonMode(face, mode);
   }
}


// ----------------------------------------------------------------------------
// selects a Matrix Mode for opengl
// ----------------------------------------------------------------------------
void Engine::setMatrixMode(int mode)
{
   if (mode != matrixMode) {
      matrixMode = mode;
      glMatrixMode(mode);
   }
}


// ----------------------------------------------------------------------------
// sets a shademode for opengl
// ----------------------------------------------------------------------------
void Engine::setShadeMode(int mode)
{
   if (mode != shadeMode) {
      shadeMode = mode;
      glShadeModel(mode);
   }
}

// -------------------------------------------------------------------------
// Setzt alle verfügbaren Daten des Meshs
// -------------------------------------------------------------------------
void Engine::set(int numFaces, Face *faces, int numEdges, Edge *edges, int numPoints, Vec3 *points, Vec3 *normals, rgbf *rgbvalues, Tex2 *textures)
{
   this->numFaces  = numFaces;
   this->faces     = faces;
   this->numEdges  = numEdges;
   this->edges     = edges;
   this->numPoints = numPoints;
   this->points    = points;
   this->normals   = normals;
   this->rgbvalues = rgbvalues;
   this->textures  = textures;

   if (lightmap) delete[] lightmap;
   lightmap = new Vec3[numPoints];
}

// -------------------------------------------------------------------------
// Setzt alle Texturen
// -------------------------------------------------------------------------
void Engine::set(Texture *bumpTex, Texture *diffTex, Texture *specTex)
{
   this->bumpTex = bumpTex;
   this->diffTex = diffTex;
   this->specTex = specTex;
}


// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------
int Engine::init()
{
   // clear state variable
   x = 0;
   y = 0;
   width  = 640;
   height = 480;
   //
   angle              = 45.0f;
   frontClippingPlane = 1.0f;
   backClippingPlane  = 20.0f;

   // current states
   shadeMode   = -1;
   fillMode    = -1;
   faceMode    = -1;
   matrixMode  = -1;
   depthFunc   = -1;

   lightmap    = NULL;

   initExtensions();
   set(NULL, NULL, NULL);
   set(0, NULL, 0, NULL, 0, NULL, NULL, NULL, NULL);

   // define default input model
   this->inputModel = NO_INPUT;

   // initialize openGl extensions
   initExtensions();

   // set actual sizes
	setBounds(0, 0, 640, 480);

	// set perspective correction mode
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	// could be GL_FRONT_AND_BACK, GL_FRONT, GL_BACK  | GL_LINE, GL_FILL
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
   glShadeModel(GL_SMOOTH);

	// enable texturing
	//GLstate::enableTexture2D(true);

	// enable z-buffer depth test
   depthTest = true;
   glEnable(depthTest);

	// set z-buffer depth functionality
   setDepthFunc(GL_LEQUAL);

	// clear background
	setClearColor(0.0f, 0.0f, 0.0f);
   glClearDepth(1.0f);

	// leave initialization in a defined mode
	setMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   // create a cubemap
   cubeMap = new NormalizationCubeMap(256);

   // init fragment shader
   initATIDot3Bump();

   return TRUE;
}

// ----------------------------------------------------------------------------
// Sets light position
// ----------------------------------------------------------------------------
void Engine::setLightPosition(float x, float y, float z)
{
   lpos.x = x;
   lpos.y = y;
   lpos.z = z;

   glLightfv(GL_LIGHT0, GL_POSITION, (float *)&lpos);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
   glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
}

// ----------------------------------------------------------------------------
// Sets Color
// ----------------------------------------------------------------------------
void Engine::setColor(float r, float g, float b)
{
   color.r = r;
   color.g = g;
   color.b = b;

   glColor3f(r, g, b);
}

/******************************************************************************
 * create an ATI Fragment Shader for Dot3 Bump Mapping
 *****************************************************************************/
void Engine::initATIDot3Bump()
{
   // generate a fragment shader id
   pointLightFragmentShaderIDATI = glGenFragmentShadersATI(1);
   if (pointLightFragmentShaderIDATI == 0)
      return;

   // bind and define and fragment shader
   glBindFragmentShaderATI(pointLightFragmentShaderIDATI);
   glBeginFragmentShaderATI();

      // Texture 0, decal
      // Texture 1, specular
      // Texture 2, bump
      // Texture 3, cube (light) 

   
      // get the outputs from the vertex shader
      glSampleMapATI(GL_REG_0_ATI, GL_TEXTURE0_ARB, GL_SWIZZLE_STR_ATI);    // decal
      glSampleMapATI(GL_REG_1_ATI, GL_TEXTURE1_ARB, GL_SWIZZLE_STR_ATI);    // specular
      glSampleMapATI(GL_REG_2_ATI, GL_TEXTURE2_ARB, GL_SWIZZLE_STR_ATI);    // N from bump map
      glSampleMapATI(GL_REG_3_ATI, GL_TEXTURE3_ARB, GL_SWIZZLE_STR_ATI);    // L (sample cubemap normalizer)

      //glPassTexCoordATI(GL_REG_2_ATI, GL_TEXTURE2_ARB, GL_SWIZZLE_STR_ATI); // light to vertex vector in light space  
      //glPassTexCoordATI(GL_REG_3_ATI, GL_TEXTURE3_ARB, GL_SWIZZLE_STR_ATI); // H

      // reg3 = N.L
      glColorFragmentOp2ATI(GL_DOT3_ATI,
                            GL_REG_3_ATI, GL_NONE, GL_NONE,
                            GL_REG_2_ATI, GL_NONE, GL_2X_BIT_ATI|GL_BIAS_BIT_ATI,
                            GL_REG_3_ATI, GL_NONE, GL_2X_BIT_ATI|GL_BIAS_BIT_ATI);

      // reg0 = 
      glColorFragmentOp2ATI(GL_ADD_ATI,
                            GL_REG_3_ATI, GL_NONE, GL_NONE,
                            GL_REG_3_ATI, GL_NONE, GL_NONE,
                            GL_CON_1_ATI, GL_NONE, GL_NONE);

      // reg1 = decal * light
      glColorFragmentOp2ATI(GL_MUL_ATI,
                            GL_REG_0_ATI, GL_NONE, GL_NONE,
                            GL_REG_0_ATI, GL_NONE, GL_NONE,
                            GL_REG_3_ATI, GL_NONE, GL_NONE);

      // reg3 = light * light
      glColorFragmentOp2ATI(GL_MUL_ATI,
                            GL_REG_3_ATI, GL_NONE, GL_NONE,
                            GL_REG_3_ATI, GL_NONE, GL_NONE,
                            GL_REG_3_ATI, GL_NONE, GL_NONE);

      // reg2 = specular * light
      glColorFragmentOp2ATI(GL_MUL_ATI,
                            GL_REG_1_ATI, GL_NONE, GL_NONE,
                            GL_REG_1_ATI, GL_NONE, GL_NONE,
                            GL_REG_3_ATI, GL_NONE, GL_NONE);


      // reg0 = specular + diffuse
      glColorFragmentOp2ATI(GL_ADD_ATI,
                            GL_REG_0_ATI, GL_NONE, GL_NONE,
                            GL_REG_0_ATI, GL_NONE, GL_NONE,
                            GL_REG_1_ATI, GL_NONE, GL_NONE);


   glEndFragmentShaderATI();
}


///////////////////////////////////////////////////////////////////////////
// set background color
/////////////////////////////////////////////////////////////////////////
void Engine::setClearColor(float r, float g, float b)
{
   background.r = r;
   background.g = g;
   background.b = b;
   glClearColor(background.r, background.g, background.b, 0.0f);
}

///////////////////////////////////////////////////////////////////////////
// reset all the stuff
///////////////////////////////////////////////////////////////////////////
void Engine::clear()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


///////////////////////////////////////////////////////////////////////////
// Here's Where We Do All The Drawing
///////////////////////////////////////////////////////////////////////////
int Engine::render()
{
	return TRUE;
}

// ----------------------------------------------------------------------------
// Zeichenmodus setzen
// ----------------------------------------------------------------------------
void Engine::setDrawMode(int drawMode)
{
   this->drawMode = drawMode;
}

// -------------------------------------------------------------------------
// translate
// -------------------------------------------------------------------------
void Engine::translate(float x, float y, float z)
{
   glTranslatef(x, y, z);
}

// -------------------------------------------------------------------------
// scale
// -------------------------------------------------------------------------
void Engine::scale(float x, float y, float z)
{
   glScalef(x, y, z);
}

// -------------------------------------------------------------------------
// rotate
// -------------------------------------------------------------------------
void Engine::rotate(float angle, float x, float y, float z)
{
   glRotatef(angle, x, y, z);
}


// ----------------------------------------------------------------------------
// Zeichnen
// ----------------------------------------------------------------------------
void Engine::draw()
{
   switch (drawMode) {
      case WIREFRAME: drawWireframe();
                      break;

      case COLORED:   drawColored();
                      break;

      case TEXTURED:  drawTextured();
                      break;

      case LINED:     drawLines();
                      break;

      case BUMP:      renderDot3BumpATI();
                      break;

      default:        drawLines();
                      break;

   }
}

// ----------------------------------------------------------------------------
// visualization of edges that lie on an object
// ----------------------------------------------------------------------------
void Engine::drawLines()
{
   if (points) glEnableClientState(GL_VERTEX_ARRAY);
   else        glDisableClientState(GL_VERTEX_ARRAY);

   glDisableClientState(GL_COLOR_ARRAY);
   glDisableClientState(GL_NORMAL_ARRAY);

   glDisable(GL_CULL_FACE);
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_LIGHTING);
   glDisable(GL_LIGHT0);
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

   if (points) glVertexPointer(3, GL_FLOAT, 0, points);
   if (faces)  glDrawElements(GL_TRIANGLES, numFaces * 3, GL_UNSIGNED_INT, faces);
}

// ----------------------------------------------------------------------------
// visualization of texture coordinates
// ----------------------------------------------------------------------------
void Engine::drawTexture()
{
   glEnableClientState(GL_VERTEX_ARRAY);
   glDisableClientState(GL_COLOR_ARRAY);
   glDisableClientState(GL_NORMAL_ARRAY);

   glDisable(GL_CULL_FACE);
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_LIGHTING);
   glDisable(GL_LIGHT0);
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
   
   if (textures) glVertexPointer(2, GL_FLOAT, 0, textures);
   if (faces)    glDrawElements(GL_TRIANGLES, numFaces * 3, GL_UNSIGNED_INT, faces);
}

// -------------------------------------------------------------------------
// Wireframe visualization on a filled object
// -------------------------------------------------------------------------
void Engine::drawWireframe()
{
   // enable states
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_NORMALIZE);
   glDisable(GL_CULL_FACE);
   //
   // set client states
   if (points)  glEnableClientState(GL_VERTEX_ARRAY);
   else         glDisableClientState(GL_VERTEX_ARRAY);

   if (normals) glEnableClientState(GL_NORMAL_ARRAY);
   else         glDisableClientState(GL_NORMAL_ARRAY);

   // check color array
   glDisableClientState(GL_COLOR_ARRAY);
   // use color for materials
   glEnable(GL_COLOR_MATERIAL);
   glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);
   glColorMaterial(GL_DIFFUSE, GL_DIFFUSE);

   // put polygons a little bit behind
   glEnable(GL_POLYGON_OFFSET_FILL);
   glPolygonOffset(2.01f, 2.0f);

   // draw body
   glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
   glDepthMask(GL_TRUE);
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   glColor3f(color.r, color.g, color.b);

   if (points)  glVertexPointer(3, GL_FLOAT, 0, points);
   if (normals) glNormalPointer(GL_FLOAT, 0, normals);
   if (faces)   glDrawElements(GL_TRIANGLES, numFaces * 3, GL_UNSIGNED_INT, faces);

   // no Polygon offset
   glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

   glDisable(GL_POLYGON_OFFSET_FILL);
   glDisable(GL_LIGHTING);
   glDisable(GL_COLOR_MATERIAL);
   glDisableClientState(GL_COLOR_ARRAY);
   glDisableClientState(GL_NORMAL_ARRAY);

   glColor3f(0.0f, 0.0f, 0.0f);
   glEnable(GL_LINE_WIDTH);
   glLineWidth(2.0f);
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
 
   // draw edges
   if (points) glVertexPointer(3, GL_FLOAT, 0, points);
   if (faces)  glDrawElements(GL_TRIANGLES, numFaces * 3, GL_UNSIGNED_INT, faces);
}

// ----------------------------------------------------------------------------
// draw textured object
// ----------------------------------------------------------------------------
void Engine::drawTextured()
{
   // enable states
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_NORMALIZE);
   glDisable(GL_CULL_FACE);
   //
   // set client states
   if (points)  glEnableClientState(GL_VERTEX_ARRAY);
   else         glDisableClientState(GL_VERTEX_ARRAY);

   if (normals)  glEnableClientState(GL_NORMAL_ARRAY);
   else          glDisableClientState(GL_NORMAL_ARRAY);

   if (textures) glEnableClientState(GL_TEXTURE_COORD_ARRAY);
   else          glDisableClientState(GL_TEXTURE_COORD_ARRAY);

   if (diffTex)  glEnable(GL_TEXTURE_2D);
   else          glDisable(GL_TEXTURE_2D);

   // draw body
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

   if (points)   glVertexPointer(3, GL_FLOAT, 0, points);
   if (textures) glTexCoordPointer(2, GL_FLOAT, 0, textures);
   if (normals)  glNormalPointer(GL_FLOAT, 0, normals);

   if (diffTex)  diffTex->bindTexture();

   if (faces)    glDrawElements(GL_TRIANGLES, numFaces * 3, GL_UNSIGNED_INT, faces);
}

///////////////////////////////////////////////////////////////////////////
// draw hidden lines onto a shaded object
///////////////////////////////////////////////////////////////////////////
void Engine::drawColored()
{
   // enable states
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_NORMALIZE);
   glEnable(GL_CULL_FACE);
   //
   // set client states
   if (points)     glEnableClientState(GL_VERTEX_ARRAY);
   else            glDisableClientState(GL_VERTEX_ARRAY);

   if (normals)    glEnableClientState(GL_NORMAL_ARRAY);
   else            glDisableClientState(GL_NORMAL_ARRAY);

   if (rgbvalues)  glEnableClientState(GL_COLOR_ARRAY);
   else            glDisableClientState(GL_COLOR_ARRAY);

   if (!rgbvalues) glColor3f(1.0f, 1.0f, 1.0f);

   // use color for materials
   glEnable(GL_COLOR_MATERIAL);
   glColorMaterial(GL_FRONT, GL_AMBIENT);
   glColorMaterial(GL_DIFFUSE, GL_DIFFUSE);

   // put polygons a little bit behind
   glEnable(GL_POLYGON_OFFSET_FILL);
   glPolygonOffset(2.01f, 2.0f);

   // draw body
   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   
   if (points)    glVertexPointer(3, GL_FLOAT, 0, points);
   if (rgbvalues) glColorPointer(3, GL_FLOAT, 0, (float *)rgbvalues);
   if (normals)   glNormalPointer(GL_FLOAT, 0, normals);
   if (faces)     glDrawElements(GL_TRIANGLES, numFaces * 3, GL_UNSIGNED_INT, faces);

   // no Polygon offset
   glDisable(GL_POLYGON_OFFSET_FILL);
   glDisable(GL_LIGHTING);
   glDisable(GL_COLOR_MATERIAL);
   glDisableClientState(GL_COLOR_ARRAY);
   glDisableClientState(GL_NORMAL_ARRAY);

   glColor3f(0.0f, 0.0f, 0.0f);
   glEnable(GL_LINE_WIDTH);
   glLineWidth(2.0f);

   if (points)    glVertexPointer(3, GL_FLOAT, 0, points);
   if (edges)     glDrawElements(GL_LINES, numEdges * 2, GL_UNSIGNED_INT, edges);
}

/******************************************************************************
 * ATI Dot3 Bump - Mapping
 *****************************************************************************/
void Engine::renderDot3BumpATI()
{
   glEnable(GL_FRAGMENT_SHADER_ATI);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_NORMALIZE);
   glEnable(GL_CULL_FACE);
   glDisable(GL_LIGHTING);

   glBindFragmentShaderATI(pointLightFragmentShaderIDATI);
   
   GLfloat lmodel_ambient[] = { 0.1f, 0.1f, 0.1f, 1.0f };

   glSetFragmentShaderConstantATI (GL_CON_1_ATI, lmodel_ambient);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// switch state to no lighning, no texture, just painting
	glEnable(GL_CULL_FACE);
   glCullFace(GL_BACK);

	float glmat[16];
	float gl_n_mat[16];
	Vec3 olight;
	glGetFloatv(GL_MODELVIEW_MATRIX,  gl_n_mat);
	glGetFloatv(GL_MODELVIEW_MATRIX,  glmat);
	invertMatrixf(glmat, gl_n_mat);

	// transform light position into object view
	olight.x = glmat[0] * lpos.x + glmat[4] * lpos.y + glmat[8]  * lpos.z + glmat[12];
	olight.y = glmat[1] * lpos.x + glmat[5] * lpos.y + glmat[9]  * lpos.z + glmat[13];
	olight.z = glmat[2] * lpos.x + glmat[6] * lpos.y + glmat[10] * lpos.z + glmat[14];

	if (points) glEnableClientState(GL_VERTEX_ARRAY);
	else        glDisableClientState(GL_VERTEX_ARRAY);

	glDisableClientState(GL_NORMAL_ARRAY);

   // transform light vector into tangent space
   for (int i=0; i<numPoints; i++) {
		// calculate L-Vectors
      lightmap[i] = (olight - points[i]).norm();
	}
   // Texture 0, decal
   // Texture 1, specular
   // Texture 2, bump
   // Texture 3, cube (light) 

   // set decal map
   if (diffTex && textures) {
	   glActiveTextureARB(GL_TEXTURE0_ARB);
      glClientActiveTextureARB(GL_TEXTURE0_ARB);
      
      glEnable(GL_TEXTURE_2D);
      glEnableClientState(GL_TEXTURE_COORD_ARRAY);

      diffTex->bindTexture();
      glTexCoordPointer(2, GL_FLOAT, 0, textures);
   }
   // set specular map
   if (specTex && textures) {
      glActiveTextureARB(GL_TEXTURE1_ARB);
      glClientActiveTextureARB(GL_TEXTURE1_ARB);

      glEnable(GL_TEXTURE_2D);
      glEnableClientState(GL_TEXTURE_COORD_ARRAY);

      specTex->bindTexture();
      glTexCoordPointer(2, GL_FLOAT, 0, textures);
   }
	// set normal map
   if (bumpTex && textures) {
      glActiveTextureARB(GL_TEXTURE2_ARB);
      glClientActiveTextureARB(GL_TEXTURE2_ARB);

      glEnable(GL_TEXTURE_2D);
      glEnableClientState(GL_TEXTURE_COORD_ARRAY);

      bumpTex->bindTexture();
      glTexCoordPointer(2, GL_FLOAT, 0, textures);
   }
	// set normalization cube map
   if (cubeMap && lightmap) {
	   glActiveTextureARB(GL_TEXTURE3_ARB);
	   glEnable(GL_TEXTURE_CUBE_MAP_EXT);
      glBindTexture(GL_TEXTURE_CUBE_MAP_EXT, cubeMap->texturename);

      glClientActiveTextureARB(GL_TEXTURE3_ARB);
      glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      glTexCoordPointer(3, GL_FLOAT, 0, lightmap);
   }
   // draw body
   if (points) glVertexPointer(3, GL_FLOAT, 0, points);
   if (faces)  glDrawElements(GL_TRIANGLES, numFaces * 3, GL_UNSIGNED_INT, faces);

	// bind specular map to texture 1
	glActiveTextureARB(GL_TEXTURE0_ARB);
	glDisable(GL_TEXTURE_2D);
	glClientActiveTextureARB(GL_TEXTURE0_ARB);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);

	glActiveTextureARB(GL_TEXTURE1_ARB);
	glDisable(GL_TEXTURE_2D);
	glClientActiveTextureARB(GL_TEXTURE1_ARB);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);

	glActiveTextureARB(GL_TEXTURE2_ARB);
	glDisable(GL_TEXTURE_2D);
	glClientActiveTextureARB(GL_TEXTURE2_ARB);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);

   glActiveTextureARB(GL_TEXTURE3_ARB);
	glDisable(GL_TEXTURE_CUBE_MAP_EXT);
	glClientActiveTextureARB(GL_TEXTURE3_ARB);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);

   glDisable(GL_FRAGMENT_SHADER_ATI);
}





