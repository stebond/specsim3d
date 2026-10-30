@echo off
rem Run all SPECSIM3D regression tests, one after the other.
rem test_scaling.py is Linux-only; on Windows the benchmark is docs\make_figures.py --full
cd /d "%~dp0"
if not exist specsim3d.exe (
  echo specsim3d.exe not found - run build_specsim3d.bat first.
  pause
  exit /b 1
)
set PY=python
where py >nul 2>&1 && set PY=py -3
if exist specsim_python.txt set /p PY=<specsim_python.txt
%PY% -c "import numpy, scipy" 2>nul
if errorlevel 1 (
  echo Python with numpy and scipy not found. Install: py -m pip install numpy scipy matplotlib
  pause
  exit /b 1
)
call :run test_conditional
call :run test_rotation
call :run test_bias
call :run test_vario_fortran
%PY% -c "import gstools" 2>nul
if errorlevel 1 (
  echo.
  echo test_gauss and test_cases SKIPPED: gstools not installed ^(py -m pip install gstools^)
) else (
  call :run test_gauss
  call :run test_cases
)
echo.
echo All tests done.
pause
exit /b 0

:run
echo.
echo ===== %1 =====
%PY% tests\%1.py
if errorlevel 1 echo *** %1 FAILED ***
exit /b 0
