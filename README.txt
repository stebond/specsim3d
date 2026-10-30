SPECSIM3D - spectral simulation of anisotropic 3D Gaussian random fields (unconditional and conditional)
Author: Stefano Bonduà (DICAM, University of Bologna), stefano.bondua@unibo.it
Manuscript: SPECSIM3D: a Fortran program for spectral simulation of anisotropic three-dimensional Gaussian random fields
            (submitted to Computers & Geosciences)
Licence: MIT (file LICENSE)

FOLDERS
  src/    Fortran source (specsim3d.f90)
  par/    parameter files: examples and variogram direction files
  setup_python_env.bat  installs Python + packages on a clean Windows machine
  run_all_tests.bat / run_full_validation.bat  quick tests / complete validation (tests, GSTools, figures, machine info)
  tests/  Python validation scripts (compare_gstools, test_gauss, test_cases, test_bias, test_rotation, test_conditional, ...)
          generated files go to tests/work, results to tests/results
  docs/   manual (SPECSIM3D_Manual.docx), make_figures.py (regenerates all figures), figures, tables

BUILD
  Windows (Visual Studio + Intel oneAPI, or gfortran):  build_specsim3d.bat     -> specsim3d.exe
  Linux/macOS:  gfortran -O3 -fopenmp -o specsim3d src/specsim3d.f90
  (-std=f2008 cannot be used: the binary VTK output uses CONVERT='big_endian')

USE
  specsim3d file.par [--vario variofile.par]
  Run from the folder where the output files must be written (the output base name is relative).
  Example: run_example.bat  (writes into output\)
  Everything about one run is in the parameter file (see the annotated files in par/):
    header              optional: free text (title, notes) up to a line 'START OF PARAMETERS:'; the parameters
                        follow it. Write your own comments BEFORE that line (or on the same line as a parameter).
                        Files without that line are read from the first line.
    ... model, grid, iseed nreal ipad, output base name
    nout                number of outputs (1 to 9), followed by nout lines with one code each:
                          1 = GSLIB grid <base>.out
                          2 = binary VTK <base>_<r>.vtk for every realization (ParaView)
                          3 = raw float32 <base>.bin (realizations consecutive, x fastest)
                          4 = <base>.vti, ONE VTK XML file with all the realizations (arrays sim_1 ... sim_n);
                              ParaView opens it directly, grid and spacing are in the file
    icond               0 = unconditional; 1 = conditional, and the 9 lines of the conditioning block
                        follow (data file, columns, transform, search ellipsoid, ndmax)
  --vario experimental variograms computed in Fortran -> <base>_vario.csv (see manual, chapter 5)
  Conditional simulation: kriging of the residual, super block search, optional normal score transform
  and back transform; see par/conditional/README_conditional.txt and manual, chapter 8.
  Older parameter files that end after the output base name are still read; for them the options
  --bin (.bin + binary VTK of realization 1), --vtk and --cond condfile.par select the outputs and the
  conditioning, and the default is <base>.out (GSLIB) + ASCII VTK.

PARAMETER FILES (par/)
  par/conditional/                        conditional simulation examples (data files, parameter files with conditioning block, run scripts)
  example_yao_spherical_padded_grid.par   Yao (1998) case, 64x64x16, spherical 20/10/5, N45E, ipad=1
  example_yao_spherical_legacy_grid.par   same, legacy SPECSIM2 grid rule (ipad=0)
  example_2d_exponential_nugget.par       2D case (nz=1), exponential, nugget 0.3, as the original SPECSIM2 test
  variogram_directions_yao.par            --vario file: N45E, N45W, vertical
  variogram_directions_benchmark.par      --vario file used for the speed benchmark

NOTES
- Practical ranges as in GSLIB (exponential exp(-3r), gaussian exp(-3r^2)); SPECSIM2 used exp(-h/a), exp(-h^2/a^2).
- ipad=1 (recommended): periodic grid large enough for an exact embedding. ipad=0 reproduces the legacy rule
  (np = power of 2 >= n): the field is periodic, sample mean/variance are almost deterministic.
- nz=1 gives the 2D case (validated against SPECSIM2.EXE output).
- Positive dip = major axis towards +z (verified numerically; compare with your own convention).
- GSTools 1.7.0 note: wave numbers are drawn by inversion only for Gaussian/Exponential in 1D and 2D; in 3D all
  models use MCMC (emcee). The Spherical 3D variogram is not reproduced because the numerical Hankel spectrum is
  negative at large k and spectral_rad_pdf takes abs() of it (non-integrable ~1/k tail, chain drifts to huge k).
- The Windows build (build_specsim3d.bat) was tested by the author on Windows and works.
