@echo off
rem 
rem build and install the myInterpolator software on Windows
rem 

if "%1"=="" GOTO HELP
if "%2"=="" GOTO HELP

setlocal
set BATDIR="%~dp0"
set BUILDDIR="%1"
set INSTALLDIR="%2"
set CCBYSA=%3
set CMAKEMODULESDIR=%~dp0CMakeModules

rem test if install directory exists
IF NOT EXIST %INSTALLDIR% (
    echo Installation direction %INSTALLDIR% must exist already
    goto END 
)

rem temporarily enhance CMAKE_MODULE_PATH
IF "%CMAKE_MODULE_PATH%"=="" (
    set CMAKE_MODULE_PATH=%CMAKEMODULESDIR%
) 

echo Using CMAKE_MODULE_PATH set to %CMAKE_MODULE_PATH% for the build

rem build all binaries
:BUILD

	rem 
	rem myInterpolator
	rem 
	
	set CURPROJECT=myInterpolator

:NORMALBUILD
    set CURPROJECT_BUILDDIR=%BUILDDIR%\%CURPROJECT%
    echo Doing NORMAL UNRESTRICTED build in %CURPROJECT_BUILDDIR%
    cmake -G%CMAKE_GENERATOR% -DCMAKE_INSTALL_PREFIX=%INSTALLDIR% -B %CURPROJECT_BUILDDIR% -S %BATDIR%
    GOTO COMPILE

:COMPILE
	mkdir %CURPROJECT_BUILDDIR%
	cmake --build %CURPROJECT_BUILDDIR% --config Debug --target INSTALL --parallel
	cmake --build %CURPROJECT_BUILDDIR% --config Release --target INSTALL --parallel
	
goto END

:HELP
echo.
echo usage: "%~nx0" ^<buildDir^> ^<installDir^>
goto END


:END
endlocal
@echo on
