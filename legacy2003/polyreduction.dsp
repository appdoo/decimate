# Microsoft Developer Studio Project File - Name="polyreduction" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=polyreduction - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "polyreduction.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "polyreduction.mak" CFG="polyreduction - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "polyreduction - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "polyreduction - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "polyreduction - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib opengl32.lib glu32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "polyreduction - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /G6 /Zp1 /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib opengl32.lib glut32.lib glu32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "polyreduction - Win32 Release"
# Name "polyreduction - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "engine"

# PROP Default_Filter ".cpp"
# Begin Source File

SOURCE=.\engine\display.cpp
# End Source File
# Begin Source File

SOURCE=.\engine\Engine.cpp
# End Source File
# Begin Source File

SOURCE=.\engine\ext_gl.cpp
# End Source File
# End Group
# Begin Group "math"

# PROP Default_Filter ".cpp"
# Begin Source File

SOURCE=.\math\trackball.cpp
# End Source File
# End Group
# Begin Group "loader"

# PROP Default_Filter ".cpp"
# Begin Source File

SOURCE=.\loader\InventorLdr.cpp
# End Source File
# Begin Source File

SOURCE=.\loader\lexer.cpp
# End Source File
# Begin Source File

SOURCE=.\loader\ObjLdr.cpp
# End Source File
# Begin Source File

SOURCE=.\loader\stringbuffer.cpp
# End Source File
# End Group
# Begin Group "reduce"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\reduce\minheap.cpp
# End Source File
# Begin Source File

SOURCE=.\reduce\qerror.cpp
# End Source File
# Begin Source File

SOURCE=.\reduce\reduction.cpp

!IF  "$(CFG)" == "polyreduction - Win32 Release"

!ELSEIF  "$(CFG)" == "polyreduction - Win32 Debug"

# ADD CPP /G6 /W3 /Gm- /Gi- /vd1 /GR- /Ob2

!ENDIF 

# End Source File
# End Group
# Begin Source File

SOURCE=.\main.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Group "engine_h"

# PROP Default_Filter ".h"
# Begin Source File

SOURCE=.\engine\display.h
# End Source File
# Begin Source File

SOURCE=.\engine\displayexception.h
# End Source File
# Begin Source File

SOURCE=.\engine\engine.h
# End Source File
# Begin Source File

SOURCE=.\engine\ext_gl.h
# End Source File
# Begin Source File

SOURCE=.\engine\glati.h
# End Source File
# Begin Source File

SOURCE=.\engine\glnvidia.h
# End Source File
# End Group
# Begin Group "math_h"

# PROP Default_Filter ".h"
# Begin Source File

SOURCE=.\math\algorithm3D.h
# End Source File
# Begin Source File

SOURCE=.\math\math3D.h
# End Source File
# Begin Source File

SOURCE=.\math\trackball.h
# End Source File
# End Group
# Begin Group "loader_h"

# PROP Default_Filter ".h"
# Begin Source File

SOURCE=.\loader\InventorLdr.h
# End Source File
# Begin Source File

SOURCE=.\loader\lexer.h
# End Source File
# Begin Source File

SOURCE=.\loader\ObjLdr.h
# End Source File
# Begin Source File

SOURCE=.\loader\stringbuffer.h
# End Source File
# End Group
# Begin Group "reduce_h"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\reduce\minheap.h
# End Source File
# Begin Source File

SOURCE=.\reduce\qerror.h
# End Source File
# Begin Source File

SOURCE=.\reduce\reduction.h
# End Source File
# End Group
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
