@echo off
rem Runs the Yao (1998) example: 64x64x16 grid, spherical 20/10/5, major axis N45E, 20 realizations.
rem Writes output\example_yao_spherical_padded_grid.bin (raw float32) and
rem output\example_yao_spherical_padded_grid.vti (VTK XML file with all the realizations, opens directly in ParaView).
setlocal EnableExtensions
cd /d "%~dp0"
if not exist specsim3d.exe call build_specsim3d.bat
if not exist specsim3d.exe exit /b 1
if not exist output mkdir output
cd output
..\specsim3d.exe ..\par\example_yao_spherical_padded_grid.par --vario ..\par\variogram_directions_yao.par
echo.
echo Done. Open output\example_yao_spherical_padded_grid.vti in ParaView (arrays sim_1 ... sim_20).
