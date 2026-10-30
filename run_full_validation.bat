@echo off
rem Complete validation: tests, GSTools (released + modified), figures/benchmark, machine info.
rem Run it from the Intel oneAPI command prompt so that ifx is found (its version is recorded).
rem Options: --quick  --skip-tests  --skip-gstools  --skip-figures
cd /d "%~dp0"
set PY=python
where py >nul 2>&1 && set PY=py -3
if exist specsim_python.txt set /p PY=<specsim_python.txt
%PY% run_full_validation.py %*
pause
