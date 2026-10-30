SPECSIM3D - conditional simulation examples
===========================================
Run (from this folder):   run_conditional_examples.bat   (Windows)   or   sh run_conditional_examples.sh

Usage:   specsim3d model.par
  model.par : the SPECSIM3D parameter file. After the output codes it contains the flag icond
              (0 = unconditional, 1 = conditional). With icond = 1 the file continues with the
              conditioning block (9 lines, GSLIB-style, annotated in the examples), so that the
              conditioning data and settings used for a run are always recorded in its parameter file.
  Older parameter files (without the output and icond lines) still work with the command line options
  --bin, --vtk and --cond condfile.par.

Method: kriging of the residual. Each unconditional spectral realization z_s is corrected with
    z_c(x) = z*_SK(x; data) + [ z_s(x) - z*_SK(x; z_s at the data nodes) ]
Neighbours are searched with a super block search inside the search ellipsoid (at most ndmax data).
Data are assigned to the nearest grid node (data in the same node are averaged; data outside the grid are
ignored), so data are reproduced exactly at their nodes. Nodes without neighbours remain unconditional.

Important
  * The variogram model of model.par must describe the (normal score) values of the data: for itrans = 0 the
    total sill (nugget + sills) is the variance of the data around skmean, normally 1 with Gaussian data.
  * itrans = 1 applies the GSLIB normal score transform (nscore) to the data, conditions in Gaussian space,
    and back transforms with the GSLIB backtr rules (linear interpolation in the table, tail options in
    cond.par). The model must then be the variogram of the normal scores with sill 1.
  * Choose the search radii of the order of the variogram ranges and ndmax between 12 and 32: larger values
    cost much more time and are rarely needed.
Extra output files:  <base>_nscore.trn (transform table, itrans = 1),  <base>_conddata.dat (data used).

Files
  example_cond1_gaussian_3d      : Gaussian data, exponential anisotropic model (64 x 64 x 16, 150 data)
  example_cond2_lognormal_nscore : lognormal data, normal score transform and back transform
  example_cond3_2d_nugget        : 2D grid (nz = 1), nugget effect, 80 data at off-grid coordinates
