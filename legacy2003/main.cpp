/**
 * Testtreiber für die Polygonreduzierung
 *
 *
 *
 *
 *
 * @author Andor Greissl
 *
 */
#include <windows.h>
#include <stdio.h>

// include header for rendering subsystem
#include "engine/display.h"
#include "engine/engine.h"
#include "engine/normalizationcubemap.h"
#include "engine/texture.h"
#include "math/math3d.h"
// stuff for model loading
#include "loader/InventorLdr.h"
#include "loader/ObjLdr.h"

// polygon reduzierung
#include "reduce/edgecollapse.h"
#include "reduce/cvertexpaircollection.h"
#include "reduce/cquadricerrorsystem.h"

// Textur -Atlas Generierung
#include "split/splitter.h"
#include "split/texturesample.h"
#include "chartification/chartification.h"

// prototypes
LRESULT CALLBACK WndProc(HWND	hWnd, UINT uMsg, WPARAM	wParam, LPARAM lParam);
void render();
void setupTransformation(float *coord3, int numCoord);
void setupWireframeRendering();
void setupSolidRendering();

Vec3 rectPoints[] = {
   {-0.5f,-0.5f, 0.0f},
   {-0.5f, 0.5f, 0.0f},
   { 0.5f, 0.5f, 0.0f},
   { 0.5f,-0.5f, 0.0f}
};

Tex2 rectTextures[] = {
   { 0.0f, 1.0f },
   { 1.0f, 1.0f },
   { 1.0f, 0.0f },
   { 0.0f, 0.0f }
};

Face rectFaces[] = {
   {0, 1, 2},
   {0, 2, 3},
};

// global handles for display and engine
Engine  *engine = NULL;
Display *display = NULL;

//
Vec3   *coord3, *norm3, *vertexCoord, *vertexNormals, *lightmap;
Face   *faces, *vertexFaces;
Edge   *edges;
rgbf   *rgbvalues;
Tex2   *texture;
int     numCoord, numFaces, numEdges, numVertices, numVertexFaces;
Vertex *vertices = NULL;
rgb    *texmem;
int     texsize = 512;

float  *highCoord3, *highNorm3;
int    *highFaces;
int     numHighCoord, numHighFaces, numDesiredCharts;

Texture *diffTex, *flatTex, *specTex, *bumpTex;

float scaleFactor = 1.0f;
float middleX = 0.0f;
float middleY = 0.0f;
float middleZ = 0.0f;


void createNormals(Vec3 **norm3, Vec3 *coord3, Face *face3, int numCoord, int numFaces)
{
   *norm3 = new Vec3[numCoord];
   memset(*norm3, 0, sizeof(Vec3) * numCoord);
   // loop over all faces
   for (int i=0; i<numFaces; i++) {
      // calculate face normal
      Vec3 a = coord3[face3[i].b] - coord3[face3[i].a];
      Vec3 b = coord3[face3[i].c] - coord3[face3[i].a];

      Vec3 nrm = (a ^ b).norm();
      // cumulate face normal
      (*norm3)[face3[i].a] += nrm;
      (*norm3)[face3[i].b] += nrm;
      (*norm3)[face3[i].c] += nrm;
   }
   // bring normals to unit length
   for (int i=0; i<numCoord; i++)
   {
      (*norm3)[i].norm();
   }
}



int main()
{
   MSG  msg;
   char name[256];
   char path[1024];
   int  npoints;
   int  option;
   memset(&msg, 0, sizeof(MSG));

   coord3 = NULL;
   faces  = NULL;

   printf("Bitte geben Sie eine Modelldatei an: ");
   scanf("%s", name);

   try
   {
      strcpy (path, "e:\\models\\");
      strcat (path, name);

      printf ("Lade Modeldatei '%s'\n", path);

      if (strstr(name, ".iv"))
         InventorLdr::load(path, (float **)&coord3, (int **)&faces, &numCoord, &numFaces);
      else if (strstr(name, ".obj"))
         ObjLdr::load(path, (float **) &coord3, (int **)&faces, &numCoord, &numFaces);
      else
      {
         printf("Modelldatei konnte nicht geöffnet werden.\n\n");
         return 0;
      }

      do {
         printf("Model [Punkte %d, Dreiecke %d]\n\n", numCoord, numFaces);

         // Display a menu
         printf("Menue:\n");
         printf("1 Polygon reduction\n");
         printf("2 Feature curves\n");
         printf("3 Distance to features\n");
         printf("4 Charts\n");
         printf("5 Seeds\n");
         printf("6 Textur Atlas\n");
         printf("7 Textur Atlas inklusive Normalen\n");
         printf("8 Texturiertes Modell\n");
         printf("9 Bump Mapping\n\n");
         //
         printf("Choose an option: ");
         scanf("%d", &option);

         switch (option) {
            case 1: 
               {
               // Klasse für Polygonreduzierung instanziieren
               printf("Vorbereitung fuer Polygonreduzierung\n\n");
               printf("Bitte geben Sie die gewuenschte Punktanzahl an: ");
               scanf("%d", &npoints);

               decimate::CQuadricErrorSystem errorSystem;
               decimate::CVertexPairCollection vertexPairCollection;
               decimate::CEdgeCollapse edgeCollapse(numCoord, (Vec3 *) coord3, numFaces, (Face *) faces, &errorSystem, &vertexPairCollection);
               delete[] coord3;
               delete[] faces;
               edgeCollapse.DecimateToPoints(npoints, &numCoord, (Vec3 **) &coord3, &numFaces, (Face **) &faces);

               printf("\nReduziertes Model [Punkte %d, Dreiecke %d]\n\n", numCoord, numFaces);
               }
               break;

            case 2:
               {
               // Feature curves
               printf("\nFeature curves\n");
               Splitter splitter((Vec3 *) coord3, numCoord, (Face *) faces, numFaces);
               splitter.featureVisualization(0.05f, 15, 5, &numEdges, &edges);
               }
               break;

            case 3:
               {
               // Distance to features
               printf("\nDistance to features\n");
               Splitter splitter((Vec3 *) coord3, numCoord, (Face *) faces, numFaces);
               splitter.distance_to_features_visualization(0.05f, 15, 5, &numEdges, &edges);
               }
               break;

            case 4:
               {
               // Splitter test
               printf("\nCharts\n");
               printf("Bitte geben Sie die gewuenschte Anzahl an Charts an: ");
               scanf("%d", &numDesiredCharts);               

               chartification::CChartification chartification(numCoord, (Vec3 *) coord3, numFaces, (Face *) faces);

               delete[] coord3;
               delete[] faces;

               chartification.Visualization(numDesiredCharts, &numCoord, (Vec3 **) &coord3, &rgbvalues, &numFaces, (Face **) &faces, &numEdges, &edges);               
               }
               break;

            case 5:
               {
               // Splitter test
               printf("\nCharts\n");
               Splitter splitter((Vec3 *) coord3, numCoord, (Face *) faces, numFaces);
               splitter.seedVisualization(0.05f, 15, 5, &numEdges, &edges, &rgbvalues);
               }
               break;

            case 6:
               {
               printf("\nCharts\n");
               printf("Bitte geben Sie die gewuenschte Anzahl an Charts an: ");
               scanf("%d", &numDesiredCharts);               

               chartification::CChartification chartification(numCoord, (Vec3 *) coord3, numFaces, (Face *) faces);

               delete[] coord3;
               delete[] faces;

               chartification.TextureVisualization(numDesiredCharts, &numCoord, (Vec3 **) &coord3, &rgbvalues, &numFaces, (Face **) &faces, &numEdges, &edges);               
               }
               break;
            case 7:
            case 8:
            case 9:
               {
               printf("\nCharts\n");
               printf("Bitte geben Sie die gewuenschte Anzahl an Charts an: ");
               scanf("%d", &numDesiredCharts);

               chartification::CChartification chartification(numCoord, (Vec3 *) coord3, numFaces, (Face *) faces);

               printf("Normalen werden berechnet\n");
               createNormals((Vec3 **)&norm3, (Vec3 *)coord3, (Face *)faces, numCoord, numFaces);

               chartification.CreateVertices(numDesiredCharts, &numVertices, &vertices, &numVertexFaces, &vertexFaces);               

               vertexCoord   = new Vec3[numVertices];
               vertexNormals = new Vec3[numVertices];
               texture       = new Tex2[numVertices];

               for (int i=0; i<numVertices; i++) {
                  vertexCoord[i]   = coord3[vertices[i].index];
                  vertexNormals[i] = norm3[vertices[i].index];
                  texture[i].u     = vertices[i].u;
                  texture[i].v     = vertices[i].v;
               }
               delete[] coord3;
               delete[] norm3;
               delete[] vertices;
               delete[] faces;

               coord3    = vertexCoord;
               norm3     = vertexNormals;
               faces     = vertexFaces;
               numFaces  = numVertexFaces;
               numCoord  = numVertices;

               // Splitter test
               printf("\nTexturgenerierung\n");
               texmem = new rgb[texsize * texsize];

               printf("Lade Originaldatei...\n");
               if (strstr(name, ".iv"))  InventorLdr::load(path, &highCoord3, &highFaces, &numHighCoord, &numHighFaces);
               else if (strstr(name, ".obj")) ObjLdr::load(path, &highCoord3, &highFaces, &numHighCoord, &numHighFaces);

               createNormals((Vec3 **)&highNorm3, (Vec3 *)highCoord3, (Face *)highFaces, numHighCoord, numHighFaces);

               normalmap::mapnormals(texsize, texmem, numCoord, (Vec3 *)coord3, (Vec3 *)norm3, texture, numFaces, (Face*) faces, (Vec3 *)highCoord3, (Vec3 *)highNorm3, numHighFaces, (Face *)highFaces);

               if (highCoord3) delete[] highCoord3;
               if (highNorm3)  delete[] highNorm3;
               if (highFaces)  delete[] highFaces;
               }
               break;
         }
      } while (option == 1);

      // normalen erzeugen
      if (!norm3) {
         printf("Normalen werden berechnet\n");
         createNormals((Vec3 **)&norm3, (Vec3 *)coord3, (Face *)faces, numCoord, numFaces);
      }
      setupTransformation((float *)coord3, numCoord);

      printf("beendnet...\n\n");
   }
   catch (ParseException e)
   {
      printf("Exception while loading:\n%s, %s\n\n", e.getMsg(), e.getMehtod());
   }
   catch (...)
   {
      Display::showMessageDialog("Unknown Exception while creating data...");
   }
   display = NULL;
   engine  = NULL;

   try 
   {
      display = new Display(WndProc, NULL, 100, 100, 800, 600, 32, false);
      display->create();

      engine = new Engine();
      engine->setInputModel(Engine::TRACKBALL_INPUT);
      engine->setBounds(0, 0, 800, 600);

      engine->setClearColor(1.0f, 1.0f, 1.0f);
      engine->setColor(0.8f, 0.8f, 0.8f);
      engine->setLightPosition(1.0f, 1.0f, 1.0f);

      switch (option) {
         case 1:
            engine->set(numFaces, faces, 0, NULL, numCoord, coord3, norm3, NULL, NULL); 
            engine->setDrawMode(Engine::WIREFRAME);
            break;

         case 2:
         case 3:
            engine->set(numFaces, faces, numEdges, edges, numCoord, coord3, norm3, rgbvalues, NULL);
            engine->setDrawMode(Engine::COLORED);
            break;

         case 4:
         case 6:
            engine->set(numFaces, faces, numEdges, edges, numCoord, coord3, norm3, rgbvalues, NULL);
            engine->setDrawMode(Engine::COLORED);
            break;

         case 5:
            engine->set(numFaces, faces, numEdges, edges, numCoord, coord3, norm3, rgbvalues, NULL);
            engine->setDrawMode(Engine::COLORED);
            break;


         case 7:
            // reset zoom / scale factors
            scaleFactor = 1.0f;
            middleX = 0.0f;
            middleY = 0.0f;
            middleZ = 0.0f;

            diffTex = new Texture(texsize, texsize, 24, texmem);
            texmem  = NULL;

            engine->set(NULL, diffTex, NULL);
            engine->set(2, rectFaces, 0, NULL, 4, rectPoints, NULL, NULL, rectTextures);
            engine->setDrawMode(Engine::TEXTURED);
            break;

         case 8:
            diffTex = new Texture(texsize, texsize, 24, texmem);
            texmem  = NULL;

            engine->setClearColor(1.0f, 1.0f, 1.0f);
            engine->set(NULL, diffTex, NULL);
            engine->set(numFaces, faces, 0, NULL, numCoord, coord3, norm3, NULL, texture);
            engine->setDrawMode(Engine::TEXTURED);
            break;

         case 9:
            specTex = Texture::readTga("e:\\models\\color\\bunny.tga");
            diffTex = Texture::readTga("e:\\models\\color\\bunny.tga");
            flatTex = Texture::readTga("e:\\models\\bump\\flat.tga");
            bumpTex = new Texture(texsize, texsize, 24, texmem);

            if (!diffTex) printf("error");

            engine->set(bumpTex, diffTex, specTex);
            engine->set(numFaces, faces, 0, NULL, numCoord, coord3, norm3, NULL, texture);
            engine->setDrawMode(Engine::BUMP);
            break;
      }
      render();

      while (GetMessage(&msg, NULL, NULL, 0))
      {
         TranslateMessage(&msg);
         DispatchMessage(&msg);
      }
   }
   catch (DisplayException e)
   {
      Display::showMessageDialog(e.getMsg());
   }
   catch (...)
   {
      Display::showMessageDialog("Display::An unknwon Exception has occured");
   }

   // cleanup
   if (diffTex)  delete diffTex;
   if (flatTex)  delete flatTex;
   if (specTex)  delete specTex;
   if (bumpTex)  delete bumpTex;

   if (engine)  delete engine;
   if (display) delete display;
   //
   if (coord3)  delete[] coord3;
   if (norm3)   delete[] norm3;
   if (faces)   delete[] faces;
   if (rgbvalues) delete[] rgbvalues;

   return 0;
}


// ----------------------------------------------------------------------------
// scale model
// ----------------------------------------------------------------------------
#define MIN(x, y)((x)<(y) ? (x):(y))
#define MAX(x, y)((x)>(y) ? (x):(y))

void setupTransformation(float *coord3, int numCoord)
{
float minx, maxx;
float miny, maxy;
float minz, maxz;

   minx = maxx = coord3[0];
   miny = maxy = coord3[1];
   minz = maxz = coord3[2];

   for (int i=1; i<numCoord; i++)
   {
      float x = coord3[i * 3];
      float y = coord3[i * 3 + 1];
      float z = coord3[i * 3 + 2];

      minx = MIN(minx, x);
      miny = MIN(miny, y);
      minz = MIN(minz, z);

      maxx = MAX(maxx, x);
      maxy = MAX(maxy, y);
      maxz = MAX(maxz, z);
   }
   scaleFactor = maxx - minx;
   scaleFactor = MAX(scaleFactor, maxy - miny);
   scaleFactor = MAX(scaleFactor, maxz - minz);

   if (scaleFactor > 0.00001f)
      scaleFactor = 2.0f / scaleFactor;

   middleX = (minx + maxx) / 2.0f;
   middleY = (miny + maxy) / 2.0f;
   middleZ = (minz + maxz) / 2.0f;

   for (int i=0; i<numCoord; i++)
   {
      coord3[i * 3]     -= middleX;
      coord3[i * 3 + 1] -= middleY;
      coord3[i * 3 + 2] -= middleZ;
/*
      coord3[i * 3]     *= scaleFactor;
      coord3[i * 3 + 1] *= scaleFactor;
      coord3[i * 3 + 2] *= scaleFactor;*/
   }
   middleX = 0.0f;
   middleY = 0.0f;
   middleZ = 0.0f;

   //scaleFactor = 1.0f;
}

// ----------------------------------------------------------------------------
// render 
//
// renders the current data
//
// ----------------------------------------------------------------------------
void render()
{
   /* test if we have a engine instance */
   if (!engine)
      return;

   engine->clear();
   engine->translate(0.0f, 0.0f, -5.0f);
   engine->scale(scaleFactor, scaleFactor, scaleFactor);

   // apply data of the the engine input - model
   engine->applyInput();

   engine->rotate(-10.0f, 1.0f, 0.0f, 0.0f);
   engine->rotate(-40.0f, 0.0f, 1.0f, 0.0f);
   engine->rotate(10.0f, 1.0f, 0.0f, 0.0f);
   engine->translate(-middleX, -middleY, -middleZ);

   engine->draw();

   display->swapBuffers();
}


// ----------------------------------------------------------------------------
// windows callback method for Display
//
// @author: Andor Greissl
// @creation: 29,02.2002
// ----------------------------------------------------------------------------
LRESULT CALLBACK WndProc(HWND	hWnd, UINT uMsg, WPARAM	wParam, LPARAM lParam)
{
static POINTS oldMousePos;
static POINTS mousePos;
static bool   shift;

	// process win32 messages
	switch (uMsg)	{
      // Watch For Window Activate Message
      case WM_ACTIVATE:
         return 0;
         break;

      case WM_KEYUP:
         //switch (MapVirtualKey(wParam, 2))
         switch (wParam)
         {
            case VK_SHIFT:
               shift = false;
               break;
            default:
              	return DefWindowProc(hWnd, uMsg, wParam, lParam);
         }
         render();
         break;

      case WM_KEYDOWN:
         //switch (MapVirtualKey(wParam, 2))
         switch (wParam)
         {
            case VK_ESCAPE:
               PostQuitMessage(0);
               break;
            case VK_F1:
               engine->setDrawMode(Engine::WIREFRAME);
               break;
            case VK_F2:
               engine->setDrawMode(Engine::COLORED);
               break;
            case VK_F3:
               engine->setDrawMode(Engine::TEXTURED);
               break;
            case VK_F4:
               engine->setDrawMode(Engine::BUMP);
               break;
            case VK_SHIFT:
               shift = true;
               break;
            default:
              	return DefWindowProc(hWnd, uMsg, wParam, lParam);
         }
         render();
         break;

      case WM_CHAR:
         //switch (MapVirtualKey(wParam, 2))
         switch (wParam)
         {
            case '+': scaleFactor += scaleFactor * .1f;
               break;
            case '-': scaleFactor -= scaleFactor * .1f;
               break;
            default:
              	return DefWindowProc(hWnd, uMsg, wParam, lParam);
         }
         render();
         break;

		//handle mouse drags
      case WM_LBUTTONDOWN:
      case WM_RBUTTONDOWN:
      case WM_MBUTTONDOWN:
         mousePos = MAKEPOINTS(lParam);
         if (engine) engine->mousePressed(mousePos.x, mousePos.y);
         break;

      case WM_MOUSEMOVE:
         if (!shift && wParam & MK_LBUTTON) {
            mousePos = MAKEPOINTS(lParam);
            if (engine) engine->mouseMoved(mousePos.x, mousePos.y);
            render();
         }
         else if (wParam & MK_LBUTTON) {
            oldMousePos = MAKEPOINTS(lParam);
            middleX += (mousePos.x - oldMousePos.x) * scaleFactor * 0.1f;
            middleY -= (mousePos.y - oldMousePos.y) * scaleFactor * 0.1f;

            mousePos = oldMousePos;
            render();
         }
         break;

      case WM_LBUTTONUP:
      case WM_RBUTTONUP:
      case WM_MBUTTONUP:
         if (!(wParam & (MK_LBUTTON|MK_MBUTTON|MK_RBUTTON))) {
            mousePos = MAKEPOINTS(lParam);
            if (engine)  engine->mouseReleased(mousePos.x, mousePos.y);
            render();
   		}
         break;

		case WM_SIZE:
         if (engine) {
            engine->setBounds(0, 0, LOWORD(lParam),HIWORD(lParam));
            render();
         }
	   	break;

      case WM_CLOSE:
         PostQuitMessage(0);
         break;

      default:
      	return DefWindowProc(hWnd, uMsg, wParam, lParam);
	}

	// Pass All Unhandled Messages To DefWindowProc
	return 0;
}

