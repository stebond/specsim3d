@echo off
rem Regenerates all manual figures (needs Python with numpy, scipy, matplotlib; gstools optional)
cd /d "%~dp0"
python make_figures.py %*
