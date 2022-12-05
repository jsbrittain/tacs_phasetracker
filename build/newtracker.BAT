@echo off
Rem ###
Rem ### Batch file to build AIonePoint.EXE (Borland C++ 32-bit)
Rem ###

Rem #
Rem # Set Target
Rem #
Set TARGET=newtracker2

Rem ####################################################
Rem ### You may want to change this for your system  ###
Rem ### This is where your Borland C++ Binaries are. ###
Rem ####################################################
set BorlandDir=C:\Progra~1\Borland\CBuilder3
if _%1_==_/?_ goto showUsage

Rem #
Rem # Set Environment
Rem #
set BorlandBinDir=%BorlandDir%\Bin\
set BorlandIncDir=%BorlandDir%\Include\
set BorlandLibDir=%BorlandDir%\Lib\
Set BorlandCC32=Bcc32.exe
Set BorlandTLink32=TLink32.exe
Set CCcfgFile=BccW32.cfg

if not exist %BorlandBinDir%%BorlandCC32% goto CompilerMissing

Rem ###
Rem ### NOTE:
Rem ### the NI-DAQ include files are in the ..\..\..\Include directory.
Rem ### the NI-DAQ libraries are in the ..\..\..\Lib directory.
Rem ###

Rem #
Rem # Set Options for Command Line Tools
Rem #
set CC32Opts=+%CCcfgFile% -IInclude;%BorlandIncDir%
set LINK32Opts=-v -Tpe -ap -c -LLib;%BorlandLibDir% -x
set LINK32Objs=c0x32.obj %TARGET%.OBJ
set LINK32Libs=nidex32b.LIB nidaq32b.LIB import32.lib cw32i.lib
set LINK32Target=%TARGET%.EXE


Rem ###
Rem ### Create Borland CC config file
Rem ###
echo -w > .\%CCcfgFile%
echo -R >>.\%CCcfgFile%
echo -v >>.\%CCcfgFile%
echo -vi>>.\%CCcfgFile%
echo -H >>.\%CCcfgFile%
echo -H=%TARGET%.csm >>.\%CCcfgFile%
echo -WC>>.\%CCcfgFile%
echo -P->>.\%CCcfgFile%
echo -c>>.\%CCcfgFile%
echo -D_RTLDLL>>.\%CCcfgFile%
Rem ###
Rem ### Clean up old stuff
Rem ###

if exist %TARGET%.OBJ del %TARGET%.OBJ
if exist %TARGET%.EXE del %TARGET%.EXE

Rem #
Rem # Now the real work!
Rem #

:Compile
echo ###
echo ### Compiling %TARGET%.c...
echo ###
%BorlandBinDir%%BorlandCC32% %CC32Opts% -o%TARGET%.OBJ %TARGET%.C
if NOT Exist %TARGET%.OBJ goto FailCompile

:Link
echo ###
echo ### Linking %TARGET%.EXE...
echo ###
%BorlandBinDir%%BorlandTLink32% %LINK32Opts% %LINK32Objs%, %LINK32Target%, , %LINK32Libs%,,
if not Exist %TARGET%.EXE goto FailLink
echo ###
echo ### %TARGET%.EXE built successfully!
echo ###
echo ### NOTE: To run this batch file again, make sure you run it
echo ###       by its FULL name '%TARGET%.BAT'. 
echo ###
GoTo OutaHere

:FailCompile
echo ###
echo ### ERROR: The source code '%TARGET%.C' failed to compile!
echo ###        Check for its existence or directory names in %TARGET%.BAT.
echo ###
GoTo OutaHere

:FailLink
echo ###
echo ### ERROR: The program '%TARGET%.EXE' failed to link!
echo ###        Check directory names in this %TARGET%.BAT file.
echo ###
GoTo OutaHere

:CompilerMissing
echo ###
echo ### ERROR: %BorlandCC32% not found!
echo ###        Please set 'BorlandDir' variable in %TARGET%.BAT to match your
echo ###        Borland C directory.
echo ###
Goto OutaHere

:showUsage
echo 
echo ###
echo ### Usage:   %TARGET%.bat [/?]
echo ###
echo ### Options: (no flags) - create %TARGET%.exe
echo ###          /?         - show this help
echo ###
echo ### Note:    Make sure the 'BorlandDir' environment variable in this
echo ###          batch file is set to your Borland C++ directory.
echo ###
echo ### -- Tips for using the Borland C++ IDE to compile this executable --
echo ###
echo ### (1) The project makes a Win32 Console Application, which can be set
echo ###     when you create a new project in the IDE via Project--New Project.
echo ### (2) The .c file is compiled with the -c flag (not linked by compiler)
echo ### (3) The compiler INCLUDE path needs ..\..\..\include appended to
echo ###     it, which can be added in the IDE via Options--Project--Directories.
echo ### (4) The linker library path needs ..\..\..\lib appended to it, which
echo ###     can be added in the IDE via Options--Project--Directories.
echo ### (5) The two import libraries that are linked are NIDEX32B.LIB
echo ###     and NIDAQ32B.LIB, which are both in the ..\..\..\lib directory.
echo ###     Note that the directory names are relative to the current path.
echo ###
goto OutaHere

Rem #
Rem # Clean up all system vars and Borland config file
Rem #
:OutaHere
if exist %CCcfgFile%  del %CCcfgFile%
if exist %TARGET%.CSM del %TARGET%.CSM
set BorlandDir=
set BorlandBinDir=
set BorlandIncDir=
set BorlandLibDir=
set BorlandCC32=
set BorlandTLink32=
set CC32Opts=
set CCcfgFile=
set LINK32Opts=
set LINK32Objs=
set LINK32Libs=
set LINK32Target=
set TARGET=

Rem ###
Rem ### End of Batch file: AIonePoint.BAT
Rem ###
