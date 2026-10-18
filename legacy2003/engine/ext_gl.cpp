/******************************************************************************
 * 
 * Datei: ext_gl.cpp
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

#include <stdio.h>
#include "ext_gl.h"


// globale variablen für OpenGl Erweiterungen *********************************

int maxMultiTextureUnits;
int maxGeneralCombiners;



// ARB_multitexture Funktionen
PFNGLMULTITEXCOORD2IARBPROC         glMultiTexCoord2iARB=NULL;
PFNGLMULTITEXCOORD2FARBPROC         glMultiTexCoord2fARB=NULL;
PFNGLMULTITEXCOORD3FARBPROC         glMultiTexCoord3fARB=NULL;
PFNGLMULTITEXCOORD3FVARBPROC        glMultiTexCoord3fvARB=NULL;
PFNGLACTIVETEXTUREARBPROC           glActiveTextureARB=NULL;

PFNGLCLIENTACTIVETEXTUREARBPROC     glClientActiveTextureARB=NULL;
PFNGLMULTITEXCOORD1DARBPROC         glMultiTexCoord1dARB=NULL;
PFNGLMULTITEXCOORD1DVARBPROC        glMultiTexCoord1dvARB=NULL;
PFNGLMULTITEXCOORD1FARBPROC         glMultiTexCoord1fARB=NULL;
PFNGLMULTITEXCOORD1FVARBPROC        glMultiTexCoord1fvARB=NULL;
PFNGLMULTITEXCOORD1IARBPROC         glMultiTexCoord1iARB=NULL;
PFNGLMULTITEXCOORD1IVARBPROC        glMultiTexCoord1ivARB=NULL;
PFNGLMULTITEXCOORD1SARBPROC         glMultiTexCoord1sARB=NULL;
PFNGLMULTITEXCOORD1SVARBPROC        glMultiTexCoord1svARB=NULL;
PFNGLMULTITEXCOORD2DARBPROC         glMultiTexCoord2dARB=NULL;
PFNGLMULTITEXCOORD2DVARBPROC        glMultiTexCoord2dvARB=NULL;
PFNGLMULTITEXCOORD2FVARBPROC        glMultiTexCoord2fvARB=NULL;
PFNGLMULTITEXCOORD2IVARBPROC        glMultiTexCoord2ivARB=NULL;
PFNGLMULTITEXCOORD2SARBPROC         glMultiTexCoord2sARB=NULL;
PFNGLMULTITEXCOORD2SVARBPROC        glMultiTexCoord2svARB=NULL;
PFNGLMULTITEXCOORD3DARBPROC         glMultiTexCoord3dARB=NULL;
PFNGLMULTITEXCOORD3DVARBPROC        glMultiTexCoord3dvARB=NULL;
PFNGLMULTITEXCOORD3IARBPROC         glMultiTexCoord3iARB=NULL;
PFNGLMULTITEXCOORD3IVARBPROC        glMultiTexCoord3ivARB=NULL;
PFNGLMULTITEXCOORD3SARBPROC         glMultiTexCoord3sARB=NULL;
PFNGLMULTITEXCOORD3SVARBPROC        glMultiTexCoord3svARB=NULL;
PFNGLMULTITEXCOORD4DARBPROC         glMultiTexCoord4dARB=NULL;
PFNGLMULTITEXCOORD4DVARBPROC        glMultiTexCoord4dvARB=NULL;
PFNGLMULTITEXCOORD4FARBPROC         glMultiTexCoord4fARB=NULL;
PFNGLMULTITEXCOORD4FVARBPROC        glMultiTexCoord4fvARB=NULL;
PFNGLMULTITEXCOORD4IARBPROC         glMultiTexCoord4iARB=NULL;
PFNGLMULTITEXCOORD4IVARBPROC        glMultiTexCoord4ivARB=NULL;
PFNGLMULTITEXCOORD4SARBPROC         glMultiTexCoord4sARB=NULL;
PFNGLMULTITEXCOORD4SVARBPROC        glMultiTexCoord4svARB=NULL;

// Open GL_EXT_vertex_shader func
PFNGLBEGINVERTEXSHADEREXTPROC         glBeginVertexShaderEXT=NULL;
PFNGLENDVERTEXSHADEREXTPROC           glEndVertexShaderEXT=NULL;
PFNGLBINDVERTEXSHADEREXTPROC          glBindVertexShaderEXT=NULL;
PFNGLGENVERTEXSHADERSEXTPROC          glGenVertexShadersEXT=NULL;
PFNGLDELETEVERTEXSHADEREXTPROC        glDeleteVertexShaderEXT=NULL;
PFNGLSHADEROP1EXTPROC                 glShaderOp1EXT=NULL;
PFNGLSHADEROP2EXTPROC                 glShaderOp2EXT=NULL;
PFNGLSHADEROP3EXTPROC                 glShaderOp3EXT=NULL;
PFNGLSWIZZLEEXTPROC                   glSwizzleEXT=NULL;
PFNGLWRITEMASKEXTPROC                 glWriteMaskEXT=NULL;
PFNGLINSERTCOMPONENTEXTPROC           glInsertComponentEXT=NULL;
PFNGLEXTRACTCOMPONENTEXTPROC          glExtractComponentEXT=NULL;
PFNGLGENSYMBOLSEXTPROC                glGenSymbolsEXT=NULL;
PFNGLSETINVARIANTEXTPROC              glSetInvariantEXT=NULL;
PFNGLSETLOCALCONSTANTEXTPROC          glSetLocalConstantEXT=NULL;
PFNGLVARIANTBVEXTPROC                 glVariantbvEXT=NULL;
PFNGLVARIANTSVEXTPROC                 glVariantsvEXT=NULL;
PFNGLVARIANTIVEXTPROC                 glVariantivEXT=NULL;
PFNGLVARIANTFVEXTPROC                 glVariantfvEXT=NULL;
PFNGLVARIANTDVEXTPROC                 glVariantdvEXT=NULL;
PFNGLVARIANTUBVEXTPROC                glVariantubvEXT=NULL;
PFNGLVARIANTUSVEXTPROC                glVariantusvEXT=NULL;
PFNGLVARIANTUIVEXTPROC                glVariantuivEXT=NULL;
PFNGLVARIANTPOINTEREXTPROC            glVariantPointerEXT=NULL;
PFNGLENABLEVARIANTCLIENTSTATEEXTPROC  glEnableVariantClientStateEXT=NULL;
PFNGLDISABLEVARIANTCLIENTSTATEEXTPROC glDisableVariantClientStateEXT=NULL;
PFNGLBINDLIGHTPARAMETEREXTPROC        glBindLightParameterEXT=NULL;
PFNGLBINDMATERIALPARAMETEREXTPROC     glBindMaterialParameterEXT=NULL;
PFNGLBINDTEXGENPARAMETEREXTPROC       glBindTexGenParameterEXT=NULL;
PFNGLBINDTEXTUREUNITPARAMETEREXTPROC  glBindTextureUnitParameterEXT=NULL;
PFNGLBINDPARAMETEREXTPROC             glBindParameterEXT=NULL;
PFNGLISVARIANTENABLEDEXTPROC          glIsVariantEnabledEXT=NULL;
PFNGLGETVARIANTBOOLEANVEXTPROC        glGetVariantBooleanvEXT=NULL;
PFNGLGETVARIANTINTEGERVEXTPROC        glGetVariantIntegervEXT=NULL;
PFNGLGETVARIANTFLOATVEXTPROC          glGetVariantFloatvEXT=NULL;
PFNGLGETVARIANTPOINTERVEXTPROC        glGetVariantPointervEXT=NULL;
PFNGLGETINVARIANTBOOLEANVEXTPROC      glGetInvariantBooleanvEXT=NULL;
PFNGLGETINVARIANTINTEGERVEXTPROC      glGetInvariantIntegervEXT=NULL;
PFNGLGETINVARIANTFLOATVEXTPROC        glGetInvariantFloatvEXT=NULL;
PFNGLGETLOCALCONSTANTBOOLEANVEXTPROC  glGetLocalConstantBooleanvEXT=NULL;
PFNGLGETLOCALCONSTANTINTEGERVEXTPROC  glGetLocalConstantIntegervEXT=NULL;
PFNGLGETLOCALCONSTANTFLOATVEXTPROC    glGetLocalConstantFloatvEXT=NULL;


// NVidia Vertex shader
PFNGLGENPROGRAMSNVPROC              glGenProgramsNV=NULL;
PFNGLDELETEPROGRAMSNVPROC           glDeleteProgramsNV=NULL;
PFNGLBINDPROGRAMNVPROC              glBindProgramNV=NULL;
PFNGLLOADPROGRAMNVPROC              glLoadProgramNV=NULL;
PFNGLPROGRAMPARAMETER4FNVPROC       glProgramParameter4fNV=NULL;
PFNGLTRACKMATRIXNVPROC              glTrackMatrixNV=NULL;

// NV_register_combiners 
PFNGLCOMBINERPARAMETERFVNVPROC               glCombinerParameter4fvNV=NULL;
PFNGLCOMBINERSTAGEPARAMETERFVNVPROC          glCombinerStageParameterfvNV=NULL;
PFNGLCOMBINERPARAMETERFVNVPROC               glCombinerParameterfvNV=NULL;
PFNGLCOMBINERPARAMETERIVNVPROC               glCombinerParameterivNV=NULL;
PFNGLCOMBINERPARAMETERFNVPROC                glCombinerParameterfNV=NULL;
PFNGLCOMBINERPARAMETERINVPROC                glCombinerParameteriNV=NULL;
PFNGLCOMBINERINPUTNVPROC                     glCombinerInputNV=NULL;
PFNGLCOMBINEROUTPUTNVPROC                    glCombinerOutputNV=NULL;
PFNGLFINALCOMBINERINPUTNVPROC                glFinalCombinerInputNV=NULL;
PFNGLGETCOMBINERINPUTPARAMETERFVNVPROC       glGetCombinerInputParameterfvNV=NULL;
PFNGLGETCOMBINERINPUTPARAMETERIVNVPROC       glGetCombinerInputParameterivNV=NULL;
PFNGLGETCOMBINEROUTPUTPARAMETERFVNVPROC      glGetCombinerOutputParameterfvNV=NULL;
PFNGLGETCOMBINEROUTPUTPARAMETERIVNVPROC      glGetCombinerOutputParameterivNV=NULL;
PFNGLGETFINALCOMBINERINPUTPARAMETERFVNVPROC  glGetFinalCombinerInputfvNV=NULL;
PFNGLGETFINALCOMBINERINPUTPARAMETERIVNVPROC  glGetFinalCombinerInputivNV=NULL;


// ATI fragment shader
PFNGLGENFRAGMENTSHADERSATIPROC        glGenFragmentShadersATI=NULL;
PFNGLBINDFRAGMENTSHADERATIPROC        glBindFragmentShaderATI=NULL;
PFNGLDELETEFRAGMENTSHADERATIPROC      glDeleteFragmentShaderATI=NULL;
PFNGLBEGINFRAGMENTSHADERATIPROC       glBeginFragmentShaderATI=NULL;
PFNGLENDFRAGMENTSHADERATIPROC         glEndFragmentShaderATI=NULL;
PFNGLPASSTEXCOORDATIPROC              glPassTexCoordATI=NULL;
PFNGLSAMPLEMAPATIPROC                 glSampleMapATI=NULL;
PFNGLCOLORFRAGMENTOP1ATIPROC          glColorFragmentOp1ATI=NULL;
PFNGLCOLORFRAGMENTOP2ATIPROC          glColorFragmentOp2ATI=NULL;
PFNGLCOLORFRAGMENTOP3ATIPROC          glColorFragmentOp3ATI=NULL;
PFNGLALPHAFRAGMENTOP1ATIPROC          glAlphaFragmentOp1ATI=NULL;
PFNGLALPHAFRAGMENTOP2ATIPROC          glAlphaFragmentOp2ATI=NULL;
PFNGLALPHAFRAGMENTOP3ATIPROC          glAlphaFragmentOp3ATI=NULL;
PFNGLSETFRAGMENTSHADERCONSTANTATIPROC glSetFragmentShaderConstantATI=NULL;

// ATI object buffers
PFNGLNEWOBJECTBUFFERATIPROC         glNewObjectBufferATI = NULL;
PFNGLISOBJECTBUFFERATIPROC          glIsObjectBufferATI = NULL;
PFNGLUPDATEOBJECTBUFFERATIPROC      glUpdateObjectBufferATI = NULL;
PFNGLGETOBJECTBUFFERFVATIPROC       glGetObjectBufferfvATI = NULL;
PFNGLGETOBJECTBUFFERIVATIPROC       glGetObjectBufferivATI = NULL;
PFNGLFREEOBJECTBUFFERATIPROC        glFreeObjectBufferATI = NULL;
PFNGLARRAYOBJECTATIPROC             glArrayObjectATI = NULL;
PFNGLGETARRAYOBJECTFVATIPROC        glGetArrayObjectfvATI = NULL;
PFNGLGETARRAYOBJECTIVATIPROC        glGetArrayObjectivATI = NULL;
PFNGLVARIANTARRAYOBJECTATIPROC      glVariantArrayObjectATI = NULL;
PFNGLGETVARIANTARRAYOBJECTFVATIPROC glGetVariantArrayObjectfvATI = NULL;
PFNGLGETVARIANTARRAYOBJECTIVATIPROC glGetVariantArrayObjectivATI = NULL;








///////////////////////////////////////////////////////////////////////////////
// prüft, ob die angegebene OpenGL Erweiterung verfügbar ist
//
// hierzu wird über glGetString(GL_EXTENSIONS) eine Liste
// aller verfügbarer Erweiterungen abgefragt und nach dem
// gewünschten String durchsucht
//
///////////////////////////////////////////////////////////////////////////////
bool extensionSupported(char *name)
{
static char *extensions = NULL;

	// first check if name is valid - no spaces inside
	char *pos = strchr(name, ' ');
	int len = strlen(name);

	if (pos || *name=='\0')	
		return false;


	// load OpenGL extension string
   if (!extensions)
	   extensions = (char *) glGetString(GL_EXTENSIONS);

   pos = extensions;

   // check for GL_GETSTRING SUCCEEDED
   if (!extensions)
      return false;

	// search the extensions
	// take care of correct start and stop condition
	for (;;) {
		if(!(pos = strstr(pos, name)))
			break;

		char *start = pos - 1;
		char *stop  = pos + len;

		if (pos==extensions  ||  *start==' ') {
			if (*stop==' ' || *stop=='\n')
				return true;
		}
		pos = stop;
	}
	return false;
}


///////////////////////////////////////////////////////////////////////////////
// Multitexture - Funktionen initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initMultiTextureARB()
{
   if (extensionSupported("GL_ARB_multitexture"))
	{
	   glMultiTexCoord2iARB    = (PFNGLMULTITEXCOORD2IARBPROC) wglGetProcAddress("glMultiTexCoord2iARB");
	   glMultiTexCoord2fARB    = (PFNGLMULTITEXCOORD2FARBPROC) wglGetProcAddress("glMultiTexCoord2fARB");
	   glMultiTexCoord3fARB    = (PFNGLMULTITEXCOORD3FARBPROC) wglGetProcAddress("glMultiTexCoord3fARB");
	   glMultiTexCoord3fvARB   = (PFNGLMULTITEXCOORD3FVARBPROC) wglGetProcAddress("glMultiTexCoord3fvARB");
	   glActiveTextureARB      = (PFNGLACTIVETEXTUREARBPROC) wglGetProcAddress("glActiveTextureARB");

		glClientActiveTextureARB= (PFNGLCLIENTACTIVETEXTUREARBPROC) wglGetProcAddress("glClientActiveTextureARB");
		glMultiTexCoord1dARB    = (PFNGLMULTITEXCOORD1DARBPROC) wglGetProcAddress("glMultiTexCoord1dARB");
		glMultiTexCoord1dvARB   = (PFNGLMULTITEXCOORD1DVARBPROC) wglGetProcAddress("glMultiTexCoord1dvARB");
		glMultiTexCoord1fARB    = (PFNGLMULTITEXCOORD1FARBPROC) wglGetProcAddress("glMultiTexCoord1fARB");
		glMultiTexCoord1fvARB   = (PFNGLMULTITEXCOORD1FVARBPROC) wglGetProcAddress("glMultiTexCoord1fvARB");
		glMultiTexCoord1iARB    = (PFNGLMULTITEXCOORD1IARBPROC) wglGetProcAddress("glMultiTexCoord1iARB");
		glMultiTexCoord1ivARB   = (PFNGLMULTITEXCOORD1IVARBPROC) wglGetProcAddress("glMultiTexCoord1ivARB");
		glMultiTexCoord1sARB    = (PFNGLMULTITEXCOORD1SARBPROC) wglGetProcAddress("glMultiTexCoord1sARB");
		glMultiTexCoord1svARB   = (PFNGLMULTITEXCOORD1SVARBPROC) wglGetProcAddress("glMultiTexCoord1svARB");
		glMultiTexCoord2dARB    = (PFNGLMULTITEXCOORD2DARBPROC) wglGetProcAddress("glMultiTexCoord2dARB");
		glMultiTexCoord2dvARB   = (PFNGLMULTITEXCOORD2DVARBPROC) wglGetProcAddress("glMultiTexCoord2dvARB");
		glMultiTexCoord2fARB    = (PFNGLMULTITEXCOORD2FARBPROC) wglGetProcAddress("glMultiTexCoord2fARB");
		glMultiTexCoord2fvARB   = (PFNGLMULTITEXCOORD2FVARBPROC) wglGetProcAddress("glMultiTexCoord2fvARB");
		glMultiTexCoord2iARB    = (PFNGLMULTITEXCOORD2IARBPROC) wglGetProcAddress("glMultiTexCoord2iARB");
		glMultiTexCoord2ivARB   = (PFNGLMULTITEXCOORD2IVARBPROC) wglGetProcAddress("glMultiTexCoord2ivARB");
		glMultiTexCoord2sARB    = (PFNGLMULTITEXCOORD2SARBPROC) wglGetProcAddress("glMultiTexCoord2sARB");
		glMultiTexCoord2svARB   = (PFNGLMULTITEXCOORD2SVARBPROC) wglGetProcAddress("glMultiTexCoord2svARB");
		glMultiTexCoord3dARB    = (PFNGLMULTITEXCOORD3DARBPROC) wglGetProcAddress("glMultiTexCoord3dARB");
		glMultiTexCoord3dvARB   = (PFNGLMULTITEXCOORD3DVARBPROC) wglGetProcAddress("glMultiTexCoord3dvARB");
		glMultiTexCoord3fARB    = (PFNGLMULTITEXCOORD3FARBPROC) wglGetProcAddress("glMultiTexCoord3fARB");
		glMultiTexCoord3fvARB   = (PFNGLMULTITEXCOORD3FVARBPROC) wglGetProcAddress("glMultiTexCoord3fvARB");
		glMultiTexCoord3iARB    = (PFNGLMULTITEXCOORD3IARBPROC) wglGetProcAddress("glMultiTexCoord3iARB");
		glMultiTexCoord3ivARB   = (PFNGLMULTITEXCOORD3IVARBPROC) wglGetProcAddress("glMultiTexCoord3ivARB");
		glMultiTexCoord3sARB    = (PFNGLMULTITEXCOORD3SARBPROC) wglGetProcAddress("glMultiTexCoord3sARB");
		glMultiTexCoord3svARB   = (PFNGLMULTITEXCOORD3SVARBPROC) wglGetProcAddress("glMultiTexCoord3svARB");
		glMultiTexCoord4dARB    = (PFNGLMULTITEXCOORD4DARBPROC) wglGetProcAddress("glMultiTexCoord4dARB");
		glMultiTexCoord4dvARB   = (PFNGLMULTITEXCOORD4DVARBPROC) wglGetProcAddress("glMultiTexCoord4dvARB");
		glMultiTexCoord4fARB    = (PFNGLMULTITEXCOORD4FARBPROC) wglGetProcAddress("glMultiTexCoord4fARB");
		glMultiTexCoord4fvARB   = (PFNGLMULTITEXCOORD4FVARBPROC) wglGetProcAddress("glMultiTexCoord4fvARB");
		glMultiTexCoord4iARB    = (PFNGLMULTITEXCOORD4IARBPROC) wglGetProcAddress("glMultiTexCoord4iARB");
		glMultiTexCoord4ivARB   = (PFNGLMULTITEXCOORD4IVARBPROC) wglGetProcAddress("glMultiTexCoord4ivARB");
		glMultiTexCoord4sARB    = (PFNGLMULTITEXCOORD4SARBPROC) wglGetProcAddress("glMultiTexCoord4sARB");
		glMultiTexCoord4svARB   = (PFNGLMULTITEXCOORD4SVARBPROC) wglGetProcAddress("glMultiTexCoord4svARB");

		glGetIntegerv(GL_MAX_TEXTURE_UNITS_ARB, &maxMultiTextureUnits);

		return true;
	}
	return false;
}


///////////////////////////////////////////////////////////////////////////////
// NVidia Register Combiner Funktionen initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initNVRegisterCombiners()
{
   if (extensionSupported("GL_NV_register_combiners"))
   {
		glCombinerParameter4fvNV         = (PFNGLCOMBINERPARAMETERFVNVPROC) wglGetProcAddress("glCombinerParameter4fvNV");
		glCombinerStageParameterfvNV     = (PFNGLCOMBINERSTAGEPARAMETERFVNVPROC) wglGetProcAddress("glCombinerStageParameterfvNV");
		glCombinerParameterfvNV          = (PFNGLCOMBINERPARAMETERFVNVPROC) wglGetProcAddress("glCombinerParameterfvNV");
		glCombinerParameterivNV          = (PFNGLCOMBINERPARAMETERIVNVPROC) wglGetProcAddress("glCombinerParameterivNV");
		glCombinerParameterfNV           = (PFNGLCOMBINERPARAMETERFNVPROC) wglGetProcAddress("glCombinerParameterfNV");
		glCombinerParameteriNV           = (PFNGLCOMBINERPARAMETERINVPROC) wglGetProcAddress("glCombinerParameteriNV");
		glCombinerInputNV                = (PFNGLCOMBINERINPUTNVPROC) wglGetProcAddress("glCombinerInputNV");
		glCombinerOutputNV               = (PFNGLCOMBINEROUTPUTNVPROC) wglGetProcAddress("glCombinerOutputNV");
		glFinalCombinerInputNV           = (PFNGLFINALCOMBINERINPUTNVPROC) wglGetProcAddress("glFinalCombinerInputNV");
		glGetCombinerInputParameterfvNV  = (PFNGLGETCOMBINERINPUTPARAMETERFVNVPROC) wglGetProcAddress("glGetCombinerInputParameterfvNV");
		glGetCombinerInputParameterivNV  = (PFNGLGETCOMBINERINPUTPARAMETERIVNVPROC) wglGetProcAddress("glGetCombinerInputParameterivNV");
		glGetCombinerOutputParameterfvNV = (PFNGLGETCOMBINEROUTPUTPARAMETERFVNVPROC) wglGetProcAddress("glGetCombinerOutputParameterfvNV");
		glGetCombinerOutputParameterivNV = (PFNGLGETCOMBINEROUTPUTPARAMETERIVNVPROC) wglGetProcAddress("glGetCombinerOutputParameterivNV");
		glGetFinalCombinerInputfvNV      = (PFNGLGETFINALCOMBINERINPUTPARAMETERFVNVPROC) wglGetProcAddress("glGetFinalCombinerInputfvNV");
		glGetFinalCombinerInputivNV      = (PFNGLGETFINALCOMBINERINPUTPARAMETERIVNVPROC) wglGetProcAddress("glGetFinalCombinerInputivNV");

		glGetIntegerv(GL_MAX_GENERAL_COMBINERS_NV, &maxGeneralCombiners);

		return true;
	}
	return false;
}


///////////////////////////////////////////////////////////////////////////////
// NVidia Vertex Shader initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initNVVertexShader()
{
   if(extensionSupported("GL_NV_vertex_program")) 
   {
		glGenProgramsNV         = (PFNGLGENPROGRAMSNVPROC)wglGetProcAddress("glGenProgramsNV");
		glDeleteProgramsNV      = (PFNGLDELETEPROGRAMSNVPROC)wglGetProcAddress("glDeleteProgramsNV");
		glBindProgramNV         = (PFNGLBINDPROGRAMNVPROC)wglGetProcAddress("glBindProgramNV");
		glLoadProgramNV         = (PFNGLLOADPROGRAMNVPROC)wglGetProcAddress("glLoadProgramNV");
		glProgramParameter4fNV  = (PFNGLPROGRAMPARAMETER4FNVPROC)wglGetProcAddress("glProgramParameter4fNV");
		glTrackMatrixNV         = (PFNGLTRACKMATRIXNVPROC)wglGetProcAddress("glTrackMatrixNV");

		return true;
	}
	return false;
}




///////////////////////////////////////////////////////////////////////////////
// GL_EXT_vertex_program initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initVertexShader()
{
   if (extensionSupported("GL_EXT_vertex_shader"))
   {
      glBeginVertexShaderEXT         = (PFNGLBEGINVERTEXSHADEREXTPROC) wglGetProcAddress("glBeginVertexShaderEXT");
      glEndVertexShaderEXT           = (PFNGLENDVERTEXSHADEREXTPROC) wglGetProcAddress("glEndVertexShaderEXT");
      glBindVertexShaderEXT          = (PFNGLBINDVERTEXSHADEREXTPROC) wglGetProcAddress("glBindVertexShaderEXT");
      glGenVertexShadersEXT          = (PFNGLGENVERTEXSHADERSEXTPROC) wglGetProcAddress("glGenVertexShadersEXT");
      glDeleteVertexShaderEXT        = (PFNGLDELETEVERTEXSHADEREXTPROC) wglGetProcAddress("glDeleteVertexShaderEXT");
      glShaderOp1EXT                 = (PFNGLSHADEROP1EXTPROC) wglGetProcAddress("glShaderOp1EXT");
      glShaderOp2EXT                 = (PFNGLSHADEROP2EXTPROC) wglGetProcAddress("glShaderOp2EXT");
      glShaderOp3EXT                 = (PFNGLSHADEROP3EXTPROC) wglGetProcAddress("glShaderOp3EXT");
      glSwizzleEXT                   = (PFNGLSWIZZLEEXTPROC) wglGetProcAddress("glSwizzleEXT");
      glWriteMaskEXT                 = (PFNGLWRITEMASKEXTPROC) wglGetProcAddress("glWriteMaskEXT");
      glInsertComponentEXT           = (PFNGLINSERTCOMPONENTEXTPROC) wglGetProcAddress("glInsertComponentEXT");
      glExtractComponentEXT          = (PFNGLEXTRACTCOMPONENTEXTPROC) wglGetProcAddress("glExtractComponentEXT");
      glGenSymbolsEXT                = (PFNGLGENSYMBOLSEXTPROC) wglGetProcAddress("glGenSymbolsEXT");
      glSetInvariantEXT              = (PFNGLSETINVARIANTEXTPROC) wglGetProcAddress("glSetInvariantEXT");
      glSetLocalConstantEXT          = (PFNGLSETLOCALCONSTANTEXTPROC) wglGetProcAddress("glSetLocalConstantEXT");
      glVariantbvEXT                 = (PFNGLVARIANTBVEXTPROC) wglGetProcAddress("glVariantbvEXT");
      glVariantsvEXT                 = (PFNGLVARIANTSVEXTPROC) wglGetProcAddress("glVariantsvEXT");
      glVariantivEXT                 = (PFNGLVARIANTIVEXTPROC) wglGetProcAddress("glVariantivEXT");
      glVariantfvEXT                 = (PFNGLVARIANTFVEXTPROC) wglGetProcAddress("glVariantfvEXT");
      glVariantdvEXT                 = (PFNGLVARIANTDVEXTPROC) wglGetProcAddress("glVariantdvEXT");
      glVariantubvEXT                = (PFNGLVARIANTUBVEXTPROC) wglGetProcAddress("glVariantubvEXT");
      glVariantusvEXT                = (PFNGLVARIANTUSVEXTPROC) wglGetProcAddress("glVariantusvEXT");
      glVariantuivEXT                = (PFNGLVARIANTUIVEXTPROC) wglGetProcAddress("glVariantuivEXT");
      glVariantPointerEXT            = (PFNGLVARIANTPOINTEREXTPROC) wglGetProcAddress("glVariantPointerEXT");
      glEnableVariantClientStateEXT  = (PFNGLENABLEVARIANTCLIENTSTATEEXTPROC) wglGetProcAddress("glEnableVariantClientStateEXT");
      glDisableVariantClientStateEXT = (PFNGLDISABLEVARIANTCLIENTSTATEEXTPROC) wglGetProcAddress("glDisableVariantClientStateEXT");
      glBindLightParameterEXT        = (PFNGLBINDLIGHTPARAMETEREXTPROC) wglGetProcAddress("glBindLightParameterEXT");
      glBindMaterialParameterEXT     = (PFNGLBINDMATERIALPARAMETEREXTPROC) wglGetProcAddress("glBindMaterialParameterEXT");
      glBindTexGenParameterEXT       = (PFNGLBINDTEXGENPARAMETEREXTPROC) wglGetProcAddress("glBindTexGenParameterEXT");
      glBindTextureUnitParameterEXT  = (PFNGLBINDTEXTUREUNITPARAMETEREXTPROC) wglGetProcAddress("glBindTextureUnitParameterEXT");
      glBindParameterEXT             = (PFNGLBINDPARAMETEREXTPROC) wglGetProcAddress("glBindParameterEXT");
      glIsVariantEnabledEXT          = (PFNGLISVARIANTENABLEDEXTPROC) wglGetProcAddress("glIsVariantEnabledEXT");
      glGetVariantBooleanvEXT        = (PFNGLGETVARIANTBOOLEANVEXTPROC) wglGetProcAddress("glGetVariantBooleanvEXT");
      glGetVariantIntegervEXT        = (PFNGLGETVARIANTINTEGERVEXTPROC) wglGetProcAddress("glGetVariantIntegervEXT");
      glGetVariantFloatvEXT          = (PFNGLGETVARIANTFLOATVEXTPROC) wglGetProcAddress("glGetVariantFloatvEXT");
      glGetVariantPointervEXT        = (PFNGLGETVARIANTPOINTERVEXTPROC) wglGetProcAddress("glGetVariantPointervEXT");
      glGetInvariantBooleanvEXT      = (PFNGLGETINVARIANTBOOLEANVEXTPROC) wglGetProcAddress("glGetInvariantBooleanvEXT");
      glGetInvariantIntegervEXT      = (PFNGLGETINVARIANTINTEGERVEXTPROC) wglGetProcAddress("glGetInvariantIntegervEXT");
      glGetInvariantFloatvEXT        = (PFNGLGETINVARIANTFLOATVEXTPROC) wglGetProcAddress("glGetInvariantFloatvEXT");
      glGetLocalConstantBooleanvEXT  = (PFNGLGETLOCALCONSTANTBOOLEANVEXTPROC) wglGetProcAddress("glGetLocalConstantBooleanvEXT");
      glGetLocalConstantIntegervEXT  = (PFNGLGETLOCALCONSTANTINTEGERVEXTPROC) wglGetProcAddress("glGetLocalConstantIntegervEXT");
      glGetLocalConstantFloatvEXT    = (PFNGLGETLOCALCONSTANTFLOATVEXTPROC) wglGetProcAddress("glGetLocalConstantFloatvEXT");

      return true;
	}
	return false;
}


///////////////////////////////////////////////////////////////////////////////
// Ati Fragment Shader Funktionen initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initAtiFragmentShaders()
{
   if (extensionSupported("GL_ATI_fragment_shader"))
   {
      glGenFragmentShadersATI        = (PFNGLGENFRAGMENTSHADERSATIPROC) wglGetProcAddress("glGenFragmentShadersATI");
      glBindFragmentShaderATI        = (PFNGLBINDFRAGMENTSHADERATIPROC) wglGetProcAddress("glBindFragmentShaderATI");
      glDeleteFragmentShaderATI      = (PFNGLDELETEFRAGMENTSHADERATIPROC) wglGetProcAddress("glDeleteFragmentShaderATI");
      glBeginFragmentShaderATI       = (PFNGLBEGINFRAGMENTSHADERATIPROC) wglGetProcAddress("glBeginFragmentShaderATI");
      glEndFragmentShaderATI         = (PFNGLENDFRAGMENTSHADERATIPROC) wglGetProcAddress("glEndFragmentShaderATI");
      glPassTexCoordATI              = (PFNGLPASSTEXCOORDATIPROC) wglGetProcAddress("glPassTexCoordATI");
      glSampleMapATI                 = (PFNGLSAMPLEMAPATIPROC) wglGetProcAddress("glSampleMapATI");
      glColorFragmentOp1ATI          = (PFNGLCOLORFRAGMENTOP1ATIPROC) wglGetProcAddress("glColorFragmentOp1ATI");
      glColorFragmentOp2ATI          = (PFNGLCOLORFRAGMENTOP2ATIPROC) wglGetProcAddress("glColorFragmentOp2ATI");
      glColorFragmentOp3ATI          = (PFNGLCOLORFRAGMENTOP3ATIPROC) wglGetProcAddress("glColorFragmentOp3ATI");
      glAlphaFragmentOp1ATI          = (PFNGLALPHAFRAGMENTOP1ATIPROC) wglGetProcAddress("glAlphaFragmentOp1ATI");
      glAlphaFragmentOp2ATI          = (PFNGLALPHAFRAGMENTOP2ATIPROC) wglGetProcAddress("glAlphaFragmentOp2ATI");
      glAlphaFragmentOp3ATI          = (PFNGLALPHAFRAGMENTOP3ATIPROC) wglGetProcAddress("glAlphaFragmentOp3ATI");
      glSetFragmentShaderConstantATI = (PFNGLSETFRAGMENTSHADERCONSTANTATIPROC) wglGetProcAddress("glSetFragmentShaderConstantATI");

      return true;
   }
	return false;
}


///////////////////////////////////////////////////////////////////////////////
// Ati Vertex Array object initialisieren
///////////////////////////////////////////////////////////////////////////////
bool initAtiVertexArrayObject()
{
   if (extensionSupported("GL_ATI_vertex_array_object"))
   {
      glNewObjectBufferATI          = (PFNGLNEWOBJECTBUFFERATIPROC) wglGetProcAddress("glNewObjectBufferATI");
      glIsObjectBufferATI           = (PFNGLISOBJECTBUFFERATIPROC) wglGetProcAddress("glIsObjectBufferATI");
      glUpdateObjectBufferATI       = (PFNGLUPDATEOBJECTBUFFERATIPROC) wglGetProcAddress("glUpdateObjectBufferATI");
      glGetObjectBufferfvATI        = (PFNGLGETOBJECTBUFFERFVATIPROC) wglGetProcAddress("glGetObjectBufferfvATI");
      glGetObjectBufferivATI        = (PFNGLGETOBJECTBUFFERIVATIPROC) wglGetProcAddress("glGetObjectBufferivATI");
      glFreeObjectBufferATI         = (PFNGLFREEOBJECTBUFFERATIPROC) wglGetProcAddress("glFreeObjectBufferATI");
      glArrayObjectATI              = (PFNGLARRAYOBJECTATIPROC) wglGetProcAddress("glArrayObjectATI");
      glGetArrayObjectfvATI         = (PFNGLGETARRAYOBJECTFVATIPROC) wglGetProcAddress("glGetArrayObjectfvATI");
      glGetArrayObjectivATI         = (PFNGLGETARRAYOBJECTIVATIPROC) wglGetProcAddress("glGetArrayObjectivATI");
      glVariantArrayObjectATI       = (PFNGLVARIANTARRAYOBJECTATIPROC) wglGetProcAddress("glVariantArrayObjectATI");
      glGetVariantArrayObjectfvATI  = (PFNGLGETVARIANTARRAYOBJECTFVATIPROC) wglGetProcAddress("glGetVariantArrayObjectfvATI");
      glGetVariantArrayObjectivATI  = (PFNGLGETVARIANTARRAYOBJECTIVATIPROC) wglGetProcAddress("glGetVariantArrayObjectivATI");

      return true;
   }
   return false;
}


///////////////////////////////////////////////////////////////////////////////
// prüfen ob Texture-Cube-Mapping möglich ist
///////////////////////////////////////////////////////////////////////////////
bool isTextureCubeMapSupported()
{
  return extensionSupported("GL_EXT_texture_cube_map") != 0;
}



/******************************************************************************
 * initialize opengl extensions
 *****************************************************************************/
bool initExtensions()
{
	if (!initMultiTextureARB()) {
		printf("OpenGL GLMultiTexture Extension not available.\n");
	}

	if (!initNVRegisterCombiners()) {
		printf("OpenGL RegisterCombinersNV Extension not available.\n");
	}

	if (!initVertexShader()) {
		printf("OpenGL Vertex Shader Extension not available.\n");
	}

	if (!initNVVertexShader()) {
		printf("OpenGL NVidia Vertex Shader Extension not available.\n");
	}

   if (!isTextureCubeMapSupported()) {
		printf("OpenGL TextureCubeMap Extension not available.\n");
	}

	if (!initAtiVertexArrayObject()) {
		printf("OpenGL ATI vertex object Extension not available.\n");
	}

	if (!initAtiFragmentShaders()) {
		printf("OpenGL ATI fragment Extension not available.\n");
	}
   return true;
}




