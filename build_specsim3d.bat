@echo off
rem ============================================================================
rem  build_specsim3d.bat - builds SPECSIM3D (Fortran 90, OpenMP) on Windows. Source: src\specsim3d.f90, output: specsim3d.exe
rem  Autore: Stefano Bondua
rem
rem  Compiler search order:
rem    1) ifx / ifort already in PATH (e.g. "Intel oneAPI command prompt for Visual Studio")
rem    2) Intel oneAPI setvars.bat (loads compiler and Visual Studio linker)
rem    3) gfortran (MSYS2 / MinGW-w64)
rem
rem  Usage: build_specsim3d.bat         (optimized for THIS cpu, OpenMP: use it for the benchmarks)
rem         build_specsim3d.bat portable  (runs on any x64 CPU with SSE4.2 and copies the Intel runtime DLLs
rem                                        next to the exe, so the whole folder can be copied to another PC
rem                                        that has no compiler; slightly slower than the native build)
rem               build_specsim3d.bat debug   (no optimization, runtime checks)
rem
rem  Note: the Intel executable needs the oneAPI runtime DLLs (libiomp5md.dll,
rem  libifcoremd.dll ...). If a DLL is missing at start-up, run the program from a
rem  oneAPI prompt or run setvars.bat first.
rem ============================================================================
setlocal EnableExtensions
cd /d "%~dp0"
set "SRC=src\specsim3d.f90"
set "EXE=specsim3d.exe"
set "FC="
set "MODE=%~1"

if not exist "%SRC%" (
  echo ERROR: %SRC% not found in %CD%
  exit /b 1
)

where ifx >nul 2>&1 && set "FC=ifx"
if not defined FC where ifort >nul 2>&1 && set "FC=ifort"
if defined FC goto :build_intel

set "SETVARS="
if exist "%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat" set "SETVARS=%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat"
if exist "%ProgramFiles%\Intel\oneAPI\setvars.bat" set "SETVARS=%ProgramFiles%\Intel\oneAPI\setvars.bat"
if exist "C:\Intel\oneAPI\setvars.bat" set "SETVARS=C:\Intel\oneAPI\setvars.bat"
if not defined SETVARS goto :try_gfortran

echo Loading Intel oneAPI environment: %SETVARS%
call "%SETVARS%" intel64 >nul 2>&1
where ifx >nul 2>&1 && set "FC=ifx"
if not defined FC where ifort >nul 2>&1 && set "FC=ifort"
if defined FC goto :build_intel

:try_gfortran
where gfortran >nul 2>&1 && goto :build_gfortran
echo.
echo ERROR: no Fortran compiler found (ifx, ifort, gfortran).
echo  - Intel: install "Intel oneAPI HPC Toolkit" (includes ifx and the integration with
echo    Visual Studio; Visual Studio or Build Tools are also needed for the linker),
echo  - or gfortran (MSYS2/MinGW-w64) in PATH.
exit /b 1

:build_intel
echo Compiler: %FC%
if /i "%MODE%"=="debug" (
  %FC% /nologo /Od /Zi /traceback /check:all /Qopenmp /heap-arrays %SRC% /Fe:%EXE%
) else if /i "%MODE%"=="portable" (
  %FC% /nologo /O3 /arch:SSE4.2 /Qopenmp /heap-arrays %SRC% /Fe:%EXE%
) else (
  %FC% /nologo /O3 /QxHost /Qopenmp /heap-arrays %SRC% /Fe:%EXE%
)
if errorlevel 1 goto :fail
if /i "%MODE%"=="portable" call :copy_dlls
call :write_info
goto :ok

:build_gfortran
echo Compilatore: gfortran
if /i "%MODE%"=="debug" (
  gfortran -O0 -g -fcheck=all -fopenmp -o %EXE% %SRC%
) else (
  gfortran -O3 -march=native -fopenmp -o %EXE% %SRC%
)
if errorlevel 1 goto :fail
goto :ok

:fail
echo.
echo COMPILAZIONE FALLITA
exit /b 1

:ok
del /q *.obj *.mod *.pdb *.ilk >nul 2>&1
echo.
echo OK: creato %EXE%
echo Uso: %EXE% file.par [--vario variofile.par]
exit /b 0

:copy_dlls
rem  copies the Intel runtime DLLs (OpenMP, Fortran, math) needed by the portable exe
rem  NOTE: no paths inside parenthesised blocks ("Program Files (x86)" would close the block)
set "ROOT1=%ProgramFiles(x86)%\Intel\oneAPI"
if not exist "%ROOT1%\compiler" set "ROOT1=%ProgramFiles%\Intel\oneAPI"
set "RED="
for /d %%D in ("%ROOT1%\compiler\*") do if exist "%%D\windows\redist\intel64_win\compiler\libiomp5md.dll" set "RED=%%D\windows\redist\intel64_win\compiler"
if defined RED goto :dll_copy
echo WARNING: Intel runtime DLLs not found under oneAPI: copy libiomp5md.dll, libifcoremd.dll,
echo          libifportmd.dll, libmmd.dll, svml_dispmd.dll next to specsim3d.exe by hand.
exit /b 0
:dll_copy
for %%F in (libiomp5md.dll libifcoremd.dll libifportmd.dll libmmd.dll svml_dispmd.dll libmmdd.dll) do if exist "%RED%\%%F" copy /y "%RED%\%%F" . >nul
echo Intel runtime DLLs copied from:
echo %RED%
exit /b 0

:write_info
rem  records compiler and flags next to the exe (run_full_validation.py and make_figures.py print them in the reports)
set "BI=%~dp0build_info.txt"
%FC% /QV > "%TEMP%\specsim_qv.txt" 2>&1
set "QV="
for /f "usebackq delims=" %%L in ("%TEMP%\specsim_qv.txt") do if not defined QV set "QV=%%L"
set "FL=/O3 /QxHost /Qopenmp"
if /i "%MODE%"=="portable" set "FL=/O3 /arch:SSE4.2 /Qopenmp"
if /i "%MODE%"=="debug" set "FL=/Od /Qopenmp"
> "%BI%" echo %QV% [%FL%, %MODE% build]
exit /b 0
