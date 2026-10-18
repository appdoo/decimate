/******************************************************************************
 * 
 * Datei: ext_gl.h
 * 
 *
 * stellt einem Anwendungsprogramm die erweiterten OpenGl Funktionen
 * zur Verfügung
 *
 * insbesondere:
 *                Multitextur Fähigkeit
 *                NVidia Pixel- und Vertex-Shader
 *                Ati Fragment- und Vertex-Shader
 *
 *
 * @author: Andor Greißl
 * @date:   23.04.2003
 *
 ******************************************************************************/

#ifndef EXT_GL_H
#define EXT_GL_H

#include <windows.h>
#include <gl\gl.h>
#include <gl\glaux.h>
#include "glnvidia.h"
#include "glati.h"


// extern definitionen für OpenGl Erweiterungen *******************************


// ARB_multitexture Funktionen
extern PFNGLMULTITEXCOORD2IARBPROC         glMultiTexCoord2iARB;
extern PFNGLMULTITEXCOORD2FARBPROC         glMultiTexCoord2fARB;
extern PFNGLMULTITEXCOORD3FARBPROC         glMultiTexCoord3fARB;
extern PFNGLMULTITEXCOORD3FVARBPROC        glMultiTexCoord3fvARB;
extern PFNGLACTIVETEXTUREARBPROC           glActiveTextureARB;

extern PFNGLCLIENTACTIVETEXTUREARBPROC     glClientActiveTextureARB;
extern PFNGLMULTITEXCOORD1DARBPROC         glMultiTexCoord1dARB;
extern PFNGLMULTITEXCOORD1DVARBPROC        glMultiTexCoord1dvARB;
extern PFNGLMULTITEXCOORD1FARBPROC         glMultiTexCoord1fARB;
extern PFNGLMULTITEXCOORD1FVARBPROC        glMultiTexCoord1fvARB;
extern PFNGLMULTITEXCOORD1IARBPROC         glMultiTexCoord1iARB;
extern PFNGLMULTITEXCOORD1IVARBPROC        glMultiTexCoord1ivARB;
extern PFNGLMULTITEXCOORD1SARBPROC         glMultiTexCoord1sARB;
extern PFNGLMULTITEXCOORD1SVARBPROC        glMultiTexCoord1svARB;
extern PFNGLMULTITEXCOORD2DARBPROC         glMultiTexCoord2dARB;
extern PFNGLMULTITEXCOORD2DVARBPROC        glMultiTexCoord2dvARB;
extern PFNGLMULTITEXCOORD2FVARBPROC        glMultiTexCoord2fvARB;
extern PFNGLMULTITEXCOORD2IVARBPROC        glMultiTexCoord2ivARB;
extern PFNGLMULTITEXCOORD2SARBPROC         glMultiTexCoord2sARB;
extern PFNGLMULTITEXCOORD2SVARBPROC        glMultiTexCoord2svARB;
extern PFNGLMULTITEXCOORD3DARBPROC         glMultiTexCoord3dARB;
extern PFNGLMULTITEXCOORD3DVARBPROC        glMultiTexCoord3dvARB;
extern PFNGLMULTITEXCOORD3IARBPROC         glMultiTexCoord3iARB;
extern PFNGLMULTITEXCOORD3IVARBPROC        glMultiTexCoord3ivARB;
extern PFNGLMULTITEXCOORD3SARBPROC         glMultiTexCoord3sARB;
extern PFNGLMULTITEXCOORD3SVARBPROC        glMultiTexCoord3svARB;
extern PFNGLMULTITEXCOORD4DARBPROC         glMultiTexCoord4dARB;
extern PFNGLMULTITEXCOORD4DVARBPROC        glMultiTexCoord4dvARB;
extern PFNGLMULTITEXCOORD4FARBPROC         glMultiTexCoord4fARB;
extern PFNGLMULTITEXCOORD4FVARBPROC        glMultiTexCoord4fvARB;
extern PFNGLMULTITEXCOORD4IARBPROC         glMultiTexCoord4iARB;
extern PFNGLMULTITEXCOORD4IVARBPROC        glMultiTexCoord4ivARB;
extern PFNGLMULTITEXCOORD4SARBPROC         glMultiTexCoord4sARB;
extern PFNGLMULTITEXCOORD4SVARBPROC        glMultiTexCoord4svARB;

// Open GL_EXT_vertex_shader func
extern PFNGLBEGINVERTEXSHADEREXTPROC         glBeginVertexShaderEXT;
extern PFNGLENDVERTEXSHADEREXTPROC           glEndVertexShaderEXT;
extern PFNGLBINDVERTEXSHADEREXTPROC          glBindVertexShaderEXT;
extern PFNGLGENVERTEXSHADERSEXTPROC          glGenVertexShadersEXT;
extern PFNGLDELETEVERTEXSHADEREXTPROC        glDeleteVertexShaderEXT;
extern PFNGLSHADEROP1EXTPROC                 glShaderOp1EXT;
extern PFNGLSHADEROP2EXTPROC                 glShaderOp2EXT;
extern PFNGLSHADEROP3EXTPROC                 glShaderOp3EXT;
extern PFNGLSWIZZLEEXTPROC                   glSwizzleEXT;
extern PFNGLWRITEMASKEXTPROC                 glWriteMaskEXT;
extern PFNGLINSERTCOMPONENTEXTPROC           glInsertComponentEXT;
extern PFNGLEXTRACTCOMPONENTEXTPROC          glExtractComponentEXT;
extern PFNGLGENSYMBOLSEXTPROC                glGenSymbolsEXT;
extern PFNGLSETINVARIANTEXTPROC              glSetInvariantEXT;
extern PFNGLSETLOCALCONSTANTEXTPROC          glSetLocalConstantEXT;
extern PFNGLVARIANTBVEXTPROC                 glVariantbvEXT;
extern PFNGLVARIANTSVEXTPROC                 glVariantsvEXT;
extern PFNGLVARIANTIVEXTPROC                 glVariantivEXT;
extern PFNGLVARIANTFVEXTPROC                 glVariantfvEXT;
extern PFNGLVARIANTDVEXTPROC                 glVariantdvEXT;
extern PFNGLVARIANTUBVEXTPROC                glVariantubvEXT;
extern PFNGLVARIANTUSVEXTPROC                glVariantusvEXT;
extern PFNGLVARIANTUIVEXTPROC                glVariantuivEXT;
extern PFNGLVARIANTPOINTEREXTPROC            glVariantPointerEXT;
extern PFNGLENABLEVARIANTCLIENTSTATEEXTPROC  glEnableVariantClientStateEXT;
extern PFNGLDISABLEVARIANTCLIENTSTATEEXTPROC glDisableVariantClientStateEXT;
extern PFNGLBINDLIGHTPARAMETEREXTPROC        glBindLightParameterEXT;
extern PFNGLBINDMATERIALPARAMETEREXTPROC     glBindMaterialParameterEXT;
extern PFNGLBINDTEXGENPARAMETEREXTPROC       glBindTexGenParameterEXT;
extern PFNGLBINDTEXTUREUNITPARAMETEREXTPROC  glBindTextureUnitParameterEXT;
extern PFNGLBINDPARAMETEREXTPROC             glBindParameterEXT;
extern PFNGLISVARIANTENABLEDEXTPROC          glIsVariantEnabledEXT;
extern PFNGLGETVARIANTBOOLEANVEXTPROC        glGetVariantBooleanvEXT;
extern PFNGLGETVARIANTINTEGERVEXTPROC        glGetVariantIntegervEXT;
extern PFNGLGETVARIANTFLOATVEXTPROC          glGetVariantFloatvEXT;
extern PFNGLGETVARIANTPOINTERVEXTPROC        glGetVariantPointervEXT;
extern PFNGLGETINVARIANTBOOLEANVEXTPROC      glGetInvariantBooleanvEXT;
extern PFNGLGETINVARIANTINTEGERVEXTPROC      glGetInvariantIntegervEXT;
extern PFNGLGETINVARIANTFLOATVEXTPROC        glGetInvariantFloatvEXT;
extern PFNGLGETLOCALCONSTANTBOOLEANVEXTPROC  glGetLocalConstantBooleanvEXT;
extern PFNGLGETLOCALCONSTANTINTEGERVEXTPROC  glGetLocalConstantIntegervEXT;
extern PFNGLGETLOCALCONSTANTFLOATVEXTPROC    glGetLocalConstantFloatvEXT;


// NVidia Vertex shader
extern PFNGLGENPROGRAMSNVPROC              glGenProgramsNV;
extern PFNGLDELETEPROGRAMSNVPROC           glDeleteProgramsNV;
extern PFNGLBINDPROGRAMNVPROC              glBindProgramNV;
extern PFNGLLOADPROGRAMNVPROC              glLoadProgramNV;
extern PFNGLPROGRAMPARAMETER4FNVPROC       glProgramParameter4fNV;
extern PFNGLTRACKMATRIXNVPROC              glTrackMatrixNV;

// NV_register_combiners 
extern PFNGLCOMBINERPARAMETERFVNVPROC               glCombinerParameter4fvNV;
extern PFNGLCOMBINERSTAGEPARAMETERFVNVPROC          glCombinerStageParameterfvNV;
extern PFNGLCOMBINERPARAMETERFVNVPROC               glCombinerParameterfvNV;
extern PFNGLCOMBINERPARAMETERIVNVPROC               glCombinerParameterivNV;
extern PFNGLCOMBINERPARAMETERFNVPROC                glCombinerParameterfNV;
extern PFNGLCOMBINERPARAMETERINVPROC                glCombinerParameteriNV;
extern PFNGLCOMBINERINPUTNVPROC                     glCombinerInputNV;
extern PFNGLCOMBINEROUTPUTNVPROC                    glCombinerOutputNV;
extern PFNGLFINALCOMBINERINPUTNVPROC                glFinalCombinerInputNV;
extern PFNGLGETCOMBINERINPUTPARAMETERFVNVPROC       glGetCombinerInputParameterfvNV;
extern PFNGLGETCOMBINERINPUTPARAMETERIVNVPROC       glGetCombinerInputParameterivNV;
extern PFNGLGETCOMBINEROUTPUTPARAMETERFVNVPROC      glGetCombinerOutputParameterfvNV;
extern PFNGLGETCOMBINEROUTPUTPARAMETERIVNVPROC      glGetCombinerOutputParameterivNV;
extern PFNGLGETFINALCOMBINERINPUTPARAMETERFVNVPROC  glGetFinalCombinerInputfvNV;
extern PFNGLGETFINALCOMBINERINPUTPARAMETERIVNVPROC  glGetFinalCombinerInputivNV;


// ATI fragment shader
extern PFNGLGENFRAGMENTSHADERSATIPROC        glGenFragmentShadersATI;
extern PFNGLBINDFRAGMENTSHADERATIPROC        glBindFragmentShaderATI;
extern PFNGLDELETEFRAGMENTSHADERATIPROC      glDeleteFragmentShaderATI;
extern PFNGLBEGINFRAGMENTSHADERATIPROC       glBeginFragmentShaderATI;
extern PFNGLENDFRAGMENTSHADERATIPROC         glEndFragmentShaderATI;
extern PFNGLPASSTEXCOORDATIPROC              glPassTexCoordATI;
extern PFNGLSAMPLEMAPATIPROC                 glSampleMapATI;
extern PFNGLCOLORFRAGMENTOP1ATIPROC          glColorFragmentOp1ATI;
extern PFNGLCOLORFRAGMENTOP2ATIPROC          glColorFragmentOp2ATI;
extern PFNGLCOLORFRAGMENTOP3ATIPROC          glColorFragmentOp3ATI;
extern PFNGLALPHAFRAGMENTOP1ATIPROC          glAlphaFragmentOp1ATI;
extern PFNGLALPHAFRAGMENTOP2ATIPROC          glAlphaFragmentOp2ATI;
extern PFNGLALPHAFRAGMENTOP3ATIPROC          glAlphaFragmentOp3ATI;
extern PFNGLSETFRAGMENTSHADERCONSTANTATIPROC glSetFragmentShaderConstantATI;

// ATI object buffers
extern PFNGLNEWOBJECTBUFFERATIPROC         glNewObjectBufferATI;
extern PFNGLISOBJECTBUFFERATIPROC          glIsObjectBufferATI;
extern PFNGLUPDATEOBJECTBUFFERATIPROC      glUpdateObjectBufferATI;
extern PFNGLGETOBJECTBUFFERFVATIPROC       glGetObjectBufferfvATI;
extern PFNGLGETOBJECTBUFFERIVATIPROC       glGetObjectBufferivATI;
extern PFNGLFREEOBJECTBUFFERATIPROC        glFreeObjectBufferATI;
extern PFNGLARRAYOBJECTATIPROC             glArrayObjectATI;
extern PFNGLGETARRAYOBJECTFVATIPROC        glGetArrayObjectfvATI;
extern PFNGLGETARRAYOBJECTIVATIPROC        glGetArrayObjectivATI;
extern PFNGLVARIANTARRAYOBJECTATIPROC      glVariantArrayObjectATI;
extern PFNGLGETVARIANTARRAYOBJECTFVATIPROC glGetVariantArrayObjectfvATI;
extern PFNGLGETVARIANTARRAYOBJECTIVATIPROC glGetVariantArrayObjectivATI;





///////////////////////////////////////////////////////////////////////////////
// Multitexture - Funktionen initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initMultiTextureARB();


///////////////////////////////////////////////////////////////////////////////
// NVidia Register Combiner Funktionen initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initNVRegisterCombiners();


///////////////////////////////////////////////////////////////////////////////
// NVidia Vertex Shader initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initNVVertexShader();


///////////////////////////////////////////////////////////////////////////////
// GL_EXT_vertex_program initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initVertexShader();


///////////////////////////////////////////////////////////////////////////////
// Ati Fragment Shader Funktionen initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initAtiFragmentShaders();


///////////////////////////////////////////////////////////////////////////////
// Ati Vertex Array object initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initAtiVertexArrayObject();


///////////////////////////////////////////////////////////////////////////////
// prüfen ob Texture-Cube-Mapping möglich ist
///////////////////////////////////////////////////////////////////////////////
bool isTextureCubeMapSupported();


/******************************************************************************
 * initialize opengl extensions
 *****************************************************************************/
bool initExtensions();


#endif