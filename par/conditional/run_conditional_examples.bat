@echo off
rem Run from this folder; the executable is two folders up.
for %%n in (example_cond1_gaussian_3d example_cond2_lognormal_nscore_3d example_cond3_2d_nugget) do (
  ..\..\specsim3d.exe %%n.par
)
pause
