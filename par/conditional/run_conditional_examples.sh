#!/bin/sh
# Run from this folder; the executable is two folders up.
EXE=../../specsim3d
for n in example_cond1_gaussian_3d example_cond2_lognormal_nscore_3d example_cond3_2d_nugget; do
  $EXE $n.par
done
