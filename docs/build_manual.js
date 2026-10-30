// Builds SPECSIM3D_Manual.docx from the figures/tables produced by make_figures.py
// Author: Stefano Bonduà
const fs = require('fs'), path = require('path');
const D = require('docx');
const { Document, Packer, Paragraph, TextRun, HeadingLevel, ImageRun, Table, TableRow, TableCell, WidthType, ShadingType,
  AlignmentType, LevelFormat, Footer, PageNumber, TableOfContents, BorderStyle, PageBreak } = D;
const FIG = path.join(__dirname, 'figures'), TAB = path.join(FIG, 'tables');
const sys = JSON.parse(fs.readFileSync(path.join(FIG, 'system_info.json')));
sys.compiler = 'Intel ifx 2023.0.0, /O3 /arch:SSE4.2 /Qopenmp (portable build)';   // the validation script could not detect the compiler on the benchmark machine
const scal = JSON.parse(fs.readFileSync(path.join(FIG, 'scaling_cached.json')));
const csv = n => fs.readFileSync(path.join(TAB, n + '.csv'), 'utf8').trim().split('\n').map(l => l.match(/("[^"]*"|[^,]+)/g).map(s => s.replace(/^"|"$/g, '')));
const W = 9360 + 400; // A4 content width in DXA with 2.0 cm margins ~ 9900
const CW = 9900;

const P = (t, o = {}) => new Paragraph({ spacing: { after: 120 }, ...o, children: runs(t) });
function runs(t) { // **bold**, `code`, *italic*
  return String(t).split(/(\*\*[^*]+\*\*|`[^`]+`|\*[^*]+\*)/).filter(Boolean).map(s =>
    s.startsWith('**') ? new TextRun({ text: s.slice(2, -2), bold: true }) :
    s.startsWith('`') ? new TextRun({ text: s.slice(1, -1), font: 'Consolas', size: 19 }) :
    s.startsWith('*') ? new TextRun({ text: s.slice(1, -1), italics: true }) : new TextRun(s)); }
const H = (t, l) => new Paragraph({ heading: l, children: [new TextRun(t)], keepNext: true });
const B = t => new Paragraph({ numbering: { reference: 'b', level: 0 }, spacing: { after: 60 }, children: runs(t) });
const code = lines => lines.split('\n').map((l, i, a) => new Paragraph({ spacing: { after: i === a.length - 1 ? 140 : 0 },
  shading: { type: ShadingType.CLEAR, fill: 'F2F2F2' }, children: [new TextRun({ text: l || ' ', font: 'Consolas', size: 18 })] }));
function pngSize(f) { const b = fs.readFileSync(f); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
function fig(name, caption, wCm = 16.5) {
  const f = path.join(FIG, name + '.png'); const [w, h] = pngSize(f); const wp = Math.round(wCm / 2.54 * 96), hp = Math.round(wp * h / w);
  return [new Paragraph({ alignment: AlignmentType.CENTER, keepNext: true, spacing: { before: 120, after: 60 },
      children: [new ImageRun({ type: 'png', data: fs.readFileSync(f), transformation: { width: wp, height: hp }, altText: { title: name, description: caption, name } })] }),
    new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 }, children: [new TextRun({ text: caption, italics: true, size: 18 })] })];
}
const bd = { style: BorderStyle.SINGLE, size: 4, color: 'BBBBBB' }, borders = { top: bd, bottom: bd, left: bd, right: bd };
function table(rows, widths, caption) {
  const tot = widths.reduce((a, b) => a + b, 0);
  const mk = (r, hd) => new TableRow({ tableHeader: hd, children: r.map((c, i) => new TableCell({ borders, width: { size: widths[i], type: WidthType.DXA },
    margins: { top: 40, bottom: 40, left: 80, right: 80 }, shading: hd ? { type: ShadingType.CLEAR, fill: 'DCE6F1' } : undefined,
    children: [new Paragraph({ children: [new TextRun({ text: String(c), bold: hd, size: 17 })] })] })) });
  const out = [];
  if (caption) out.push(new Paragraph({ keepNext: true, spacing: { before: 120, after: 60 }, children: [new TextRun({ text: caption, bold: true, size: 19 })] }));
  out.push(new Table({ width: { size: tot, type: WidthType.DXA }, columnWidths: widths, rows: rows.map((r, i) => mk(r, i === 0)) }));
  out.push(new Paragraph({ spacing: { after: 160 }, children: [] }));
  return out;
}
const f4 = x => (+x).toFixed(4).replace(/0+$/, '').replace(/\.$/, '.0');

// ---------------- content ----------------
const c = [];
c.push(new Paragraph({ spacing: { before: 1800, after: 200 }, alignment: AlignmentType.CENTER, children: [new TextRun({ text: 'SPECSIM3D', bold: true, size: 64 })] }));
c.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 }, children: [new TextRun({ text: 'Three-dimensional spectral simulation of anisotropic Gaussian random fields', size: 30 })] }));
c.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 100 }, children: [new TextRun({ text: 'User and validation manual', size: 22 })] }));
c.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 1200, after: 60 }, children: [new TextRun({ text: 'Stefano Bonduà', bold: true, size: 26 })] }));
c.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 60 }, children: [new TextRun('DICAM – Department of Civil, Chemical, Environmental and Materials Engineering, University of Bologna')] }));
c.push(new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun(`Manual generated ${sys.date.slice(0, 10)}`)] }));
c.push(new Paragraph({ children: [new PageBreak()] }));
c.push(H('Contents', HeadingLevel.HEADING_1));
c.push(new TableOfContents('Contents', { hyperlink: true, headingStyleRange: '1-2' }));
c.push(new Paragraph({ children: [new PageBreak()] }));

c.push(H('1. Introduction', HeadingLevel.HEADING_1));
c.push(P('SPECSIM3D generates unconditional and conditional realizations of stationary Gaussian random fields in one, two or three dimensions with a prescribed covariance model (spherical, exponential or Gaussian, nested structures and nugget) and arbitrary 3D geometric anisotropy. It is a 3D extension of the Fortran 77 program SPECSIM2 by Pardo-Igúzquiza and Dowd (2003), which was limited to 2D, and it keeps the same spectral (Fourier integral) principle: the covariance is sampled on a grid, its Fourier transform gives the spectral density, and each realization is obtained with one inverse FFT after assigning random phases.'));
c.push(P('The program is written in modern free-form Fortran (single source file, no external libraries), uses its own radix-2 FFT parallelized with OpenMP, and provides outputs that can be opened directly in ParaView (binary legacy VTK) or in GSLIB-compatible tools. A built-in module (`--vario`) computes experimental variograms of the realizations in Fortran, which is several times faster than equivalent NumPy code.'));
c.push(P('This manual describes the method, the build and use of the program, the file formats, and the validation against analytic models, against the original SPECSIM2 executable and against the Python library GSTools. All figures and tables are produced by the script `docs/make_figures.py`.'));
c.push(H('1.1 Main features', HeadingLevel.HEADING_2));
['3D (and 2D with nz = 1) spectral simulation, any number of realizations per run, own FFT with OpenMP.',
 'Spherical, exponential and Gaussian models with GSLIB conventions (practical ranges), nested structures and nugget.',
 'Full 3D anisotropy through GSLIB angles (azimuth, dip, plunge).',
 'Exact embedding (option `ipad = 1`): the periodic grid is enlarged so that realizations have the correct covariance up to the largest lag of the output grid; the compact power-of-two grid is available as an option (`ipad = 0`).',
 'Outputs: GSLIB ASCII, raw float32 binary, ASCII or binary legacy VTK for ParaView.',
 'Experimental variogram module in Fortran (any direction expressed as an integer step in grid cells).',
 'Conditional simulation (icond = 1 in the parameter file): kriging of the residual with a GSLIB-style super block search, optional normal score transform of the data and back transform of the realizations (GSLIB nscore/backtr rules), OpenMP parallel.'].forEach(t => c.push(B(t)));

c.push(H('2. Method', HeadingLevel.HEADING_1));
c.push(P('Let C(h) be the covariance of the field on a regular grid with spacing (Δx, Δy, Δz). SPECSIM3D evaluates C on a periodic grid of N = Nx·Ny·Nz points, using for every axis the minimum-image lag (the distance to the nearest periodic copy of the origin). Because C(h) = C(−h), the discrete Fourier transform S = F[C] is real. Its non-negativity is checked and the most negative value is reported as a diagnostic (small negative values are clipped to zero). For a covariance that is exactly periodic on the grid S ≥ 0 and its sum equals C(0), the variance of the field.'));
c.push(P('A realization is obtained as z = F⁻¹[ √S · exp(iφ) ], with φ uniformly distributed in [0, 2π) and Hermitian symmetry φ(−f) = −φ(f) so that z is real. The self-conjugate frequencies (f = −f modulo the grid) receive a random sign instead of a random phase. The resulting field is Gaussian with zero mean and covariance equal to the sampled one; it is the exact discrete counterpart of the spectral representation used by Pardo-Igúzquiza and Chica-Olmo (1993) and by the circulant-embedding method (Dietrich and Newsam, 1993). Only the first Nx·Ny·Nz nodes of the periodic grid are written to the output.'));
c.push(...fig('fig02_algorithm', 'Figure 1 – SPECSIM3D workflow.', 15.5));
c.push(H('2.1 Periodic grid and the ipad option', HeadingLevel.HEADING_2));
c.push(P('A field simulated on a periodic grid is periodic: values at opposite faces are correlated even when they are far apart, and, when the periodic grid is not much larger than the output grid, the sample mean and variance of a single realization are nearly deterministic. With **ipad = 0** (compact grid) the periodic size is the smallest power of two not smaller than the grid size, which is the rule of SPECSIM2 and keeps the memory small. With **ipad = 1** (recommended) each periodic size is chosen as a power of two not smaller than max(n + reach, 2·reach), where *reach* is the extent of the covariance support (the lag at which the correlation falls below a cutoff: 1.0 of the range for the spherical model, 2.3 for the exponential and 1.6 for the Gaussian model) along that axis. The embedding is then (practically) exact: the simulated field has the correct covariance for any pair of grid nodes, and the realization-to-realization variability of mean and variance is the natural one. The price is memory: padded grids can be up to eight times larger than the output grid (Section 7).'));
c.push(...fig('fig07_legacy_vs_padded', 'Figure 2 – Left: mean variogram of a Gaussian-model field on a 64×64×16 grid with the compact (ipad = 0) and padded (ipad = 1) periodic grids; the compact grid shows the periodic wrap-around at large lags. Centre and right: sample variance and sample mean of 300 spherical-model realizations; the field on the compact grid has an almost deterministic variance.', 16.5));
c.push(H('2.2 Covariance models and ranges', HeadingLevel.HEADING_2));
c.push(P('Models follow the GSLIB convention with **practical ranges** a: with r = h/a, the spherical model is C = σ²(1 − 1.5r + 0.5r³) for r < 1; the exponential model is C = σ² exp(−3r); the Gaussian model is C = σ² exp(−3r²). SPECSIM2 used exp(−h/a) and exp(−h²/a²) instead: to reproduce an old parameter file multiply the exponential range by 3 and the Gaussian range by √3. Here r = |R·h|/a_hmax, where R is the rotation-and-scaling matrix built from the three ranges and the angles. A nugget effect is added only at zero lag.'));
c.push(H('2.3 Anisotropy and angles', HeadingLevel.HEADING_2));
c.push(P('The three ranges are *ahmax* (major horizontal direction), *ahmin* (perpendicular horizontal direction) and *avert* (third axis). Angles follow GSLIB (setrot): the azimuth is measured clockwise from North (y axis), with x = East, y = North, z = up; a positive dip tilts the major axis toward +z; the plunge rotates the minor axes around the major axis. With all angles equal to zero the major axis is the y axis. The conventions were verified numerically with an independent test (variograms along integer directions, test_rotation.py). Users should check the sign of the dip against their own convention.'));
c.push(...fig('fig01_geometry', 'Figure 3 – Anisotropy ellipsoid (ranges 20/10/5) and angle convention for two combinations of azimuth, dip and plunge.', 14.5));

c.push(H('3. Installation', HeadingLevel.HEADING_1));
c.push(H('3.1 Windows with Visual Studio and Intel oneAPI', HeadingLevel.HEADING_2));
c.push(P('Install Visual Studio (Desktop development with C++, which provides the linker) and the Intel oneAPI HPC Toolkit (ifx or ifort). Then run `build_specsim3d.bat`: it initializes the Visual Studio and oneAPI environments, compiles `specsim3d.f90` with `/O3 /QxHost /Qopenmp /heap-arrays` and produces `specsim3d.exe`. If the Intel compiler is not found the script falls back to gfortran. `run_example.bat` runs the Yao example. The Windows build was tested by the author and works.'));
c.push(H('3.2 gfortran (Linux, macOS, MSYS2)', HeadingLevel.HEADING_2));
c.push(...code('gfortran -O3 -fopenmp -o specsim3d specsim3d.f90'));
c.push(P('The option `-std=f2008` cannot be combined with the `CONVERT=` extension used for the big-endian VTK files. The number of OpenMP threads is set with the environment variable `OMP_NUM_THREADS`.'));

c.push(H('3.3 Python environment for tests, figures and benchmarks', HeadingLevel.HEADING_2));
c.push(P('The validation scripts, the figure scripts and the benchmark need Python 3.9 or later with **numpy, scipy, matplotlib, gstools and psutil** (GSTools is optional: without it the comparisons with GSTools are skipped). On a machine without Python run, from a normal Command Prompt (PowerShell is not needed to start it; the script calls it internally):'));
c.push(...code('setup_python_env.bat'));
c.push(P('The script (i) looks for a working Python 3.9+ (the Microsoft Store alias is recognised as not usable), (ii) if none is found, downloads the latest stable Python from python.org and installs it silently for the current user (no administrator rights; `/force` always installs a fresh one), (iii) installs the packages with pip using binary wheels only and, if a package has no wheel yet for the newest Python, automatically tries the previous minor versions, (iv) checks the imports and writes the interpreter path in `specsim_python.txt` (read by `run_all_tests.bat` and `run_full_validation.bat`) and the package versions in `python_env.txt`. It needs internet access to python.org and pypi.org. It does not install the Fortran compiler (Section 3.1).'));
c.push(...code('py -m pip install numpy scipy matplotlib gstools psutil     # manual alternative'));

c.push(H('4. Running the program', HeadingLevel.HEADING_1));
c.push(...code('specsim3d file.par [--vario variofile.par]'));
c.push(P('Everything that defines a run, including the outputs and the conditioning, is in the parameter file (Section 4.1), so that the file always documents how the results were produced. The only command-line option is --vario; the options of earlier versions are still accepted for parameter files without the output lines.'));
table([['Option', 'Effect'], ['--vario f', 'computes the experimental variograms defined in f and writes <base>_vario.csv'],
  ['--bin, --vtk, --cond f', 'only for parameter files that end after the base name (legacy): --bin writes <base>.bin and a binary VTK of realization 1; --vtk writes a binary VTK for every realization; --cond f conditions the realizations as described in the conditioning file f. Without options the legacy default is <base>.out and an ASCII VTK of realization 1. They are ignored (with a message) when the parameter file contains the output lines']], [2200, 7700], 'Table 1 – Command-line options.').forEach(x => c.push(x));
c.push(H('4.1 Parameter file', HeadingLevel.HEADING_2));
c.push(P('The file starts with an optional header in the GSLIB style: any lines of text (title, notes, comments) up to and including a line containing START OF PARAMETERS:. The parameters follow that line. All free text of the user (description of the case, date, author) must therefore be written before START OF PARAMETERS:; after it, only the comment on the same line of a parameter is allowed. Files without the START OF PARAMETERS: line are read from the first line.'));
c.push(...code(`Parameters for SPECSIM3D
************************

START OF PARAMETERS:
xmin ymin zmin            origin of the grid (cell centres)
dx dy dz                  grid spacing
nx ny nz                  number of nodes (nz = 1 gives the 2D case)
nugget nstruct            nugget variance, number of nested structures
itype sill ahmax ahmin avert az dip plg     one line per structure
                          itype: 1 spherical, 2 exponential, 3 gaussian
iseed nreal ipad          random seed, number of realizations, 0 compact / 1 padded grid
outbase                   base name of the output files (no comment on this line)
nout                      number of outputs
code                      nout lines, one output code each:
                          1 = GSLIB grid <base>.out, 2 = binary VTK <base>_<r>.vtk
                          for every realization, 3 = raw float32 <base>.bin,
                          4 = <base>.vti, one VTK XML file with all the realizations
icond                     0 = unconditional, 1 = conditional (conditioning block follows)
[conditioning block]      only if icond = 1: the 9 lines of Section 8.1`));
c.push(P('Every line except the base name may be followed by a free comment, which is ignored. The example parameter files in the folder `par` are annotated in this way.'));
c.push(P('Example (Yao, 1998 test case, 64×64×16 grid, spherical model, ranges 20/10/5, major axis N45E):'));
c.push(...code(`Parameters for SPECSIM3D
************************

START OF PARAMETERS:
0 0 0
1 1 1
64 64 16
0.0 1
1 1.0 20 10 5 45 0 0
31791191 20 1
yao
2
3
2
0`));
c.push(P('This file writes yao.bin and one binary VTK file per realization, unconditionally. With the code 4 instead of 2 the 20 realizations are written to a single file yao.vti.'));
c.push(P('The random number generator is the Lehmer/Park–Miller generator of SPECSIM2 (multiplier 16807), implemented with 64-bit integers; the same seed therefore always gives the same field for the same grid and parameters, independently of the number of threads.'));
c.push(H('4.2 Output files', HeadingLevel.HEADING_2));
[ '**GSLIB (.out, code 1)**: header with the title and the number of variables, followed by the values, x fastest, then y, then z, one column per realization.',
  '**Raw binary (.bin, code 3)**: float32, little-endian, realization after realization; in NumPy: `np.fromfile(f, "float32").reshape(nreal, nz, ny, nx)`.',
  '**VTK XML (.vti, code 4)**: one VTK XML ImageData file with all the realizations as point-data arrays sim_1 … sim_n, stored as raw little-endian float32. The XML header contains the dimensions, the origin and the spacing, so ParaView opens the file directly (File > Open) and the realization is chosen in the colouring menu; in Python it is read with `pyvista.read("base.vti")`. It is the recommended format for visualization. The file size is 4·nx·ny·nz·nreal bytes, so for large grids use a smaller nreal or the code 2.',
  '**VTK (code 2)**: legacy STRUCTURED_POINTS with the scalar field "sim". The binary version is written big-endian as required by the legacy format; ParaView opens it directly and loads it much faster than the ASCII one. Spacing and origin are taken from the parameter file.'].forEach(t => c.push(B(t)));
c.push(H('4.3 Screen output', HeadingLevel.HEADING_2));
c.push(P('The program prints the periodic grid size, the memory requirement, the smallest spectral value and the sum of the spectrum (which must equal the total variance), the times for covariance + spectrum, realizations and variograms, and the sample mean and variance of each realization.'));

c.push(H('5. Variogram module (--vario)', HeadingLevel.HEADING_1));
c.push(P('The module computes the experimental semivariogram γ(h) = ½·mean[(z(x+h) − z(x))²] using **all** the pairs of the grid for each lag vector. Directions are defined by an integer step (sx, sy, sz) in grid cells, and lag k corresponds to the vector k·(sx, sy, sz). This is exact for any direction that is a multiple of a grid vector, requires no tolerances and is the same estimator used for the Python validation. The computation is parallelized with OpenMP over (direction, lag) tasks.'));
c.push(...code(`ndir iperreal            number of directions; 1 = also write the variogram of every realization
sx sy sz nlag            one line per direction (e.g. 1 1 0 40 = N45E diagonal, 40 lags)`));
c.push(P('The output `<base>_vario.csv` contains the columns **dir, sx, sy, sz, lag, h, npairs, gamma_model, gamma_mean, gamma_sd** and, when iperreal = 1, **gamma_r1 … gamma_rN**. gamma_model is the analytic semivariogram of the parameter file, gamma_mean and gamma_sd are the mean and standard deviation over the realizations. The program also prints the RMS difference between gamma_mean and gamma_model on the rising part of the variogram (γ_model < 0.95 of the sill).'));
c.push(...fig('figS1_variograms_yao', 'Figure 4 – Variograms computed with --vario for the Yao case (spherical, 20 realizations): single realizations (grey), mean ± 1 standard deviation and analytic model along the three principal directions.', 17));
const t4 = csv('figS1_yao_rms');
table([['Direction', 'RMS (mean vs model)', 'Max |error|'], ...t4.slice(1).map(r => [r[0], r[1], r[2]])], [4300, 2800, 2800], 'Table 2 – Error of the mean variogram of 20 realizations (rising part).').forEach(x => c.push(x));
const t10 = csv('fig10_variogram_speed');
c.push(...fig('fig10_variogram_speed', 'Figure 5 – Time per realization to compute 3 directions × 40 lags: Fortran module vs NumPy slicing.', 9.5));
table([['n (n³ nodes)', 'Fortran (s/realization)', 'NumPy (s/realization)', 'Speed-up'], ...t10.slice(1)], [2200, 2700, 2700, 2300], 'Table 3 – Variogram computation time (' + sys.cpu + ', ' + sys.cores + ' threads).').forEach(x => c.push(x));
c.push(P('The Fortran module is 40 to 130 times faster than the NumPy estimator on this machine (20 threads); at 64³ its time is close to the resolution of the program timer, so the largest speed-up is less reliable than the ones at 128³ and 256³. The maximum difference between the Fortran and the Python variograms is 5·10⁻⁷, which is the float32 rounding of the .bin file.'));

c.push(H('6. Validation', HeadingLevel.HEADING_1));
c.push(P('**Scope and honesty note.** SPECSIM3D was written for this project; the reference for all tests is the analytic covariance model (and, for 2D, the output of the original SPECSIM2.EXE), not an independent 3D spectral code. GSTools (version ' + (sys.gstools || 'n/a') + ') is used as an independent, widely used generator.'));
c.push(H('6.1 Reproduction of the model variograms', HeadingLevel.HEADING_2));
c.push(P('The Yao (1998) test case (64×64×16, unit variance, ranges 20/10/5, N45E) was simulated with 20 realizations for each model. GSTools was configured with the equivalent models: spherical len_scale = a, exponential len_scale = a/3, Gaussian len_scale = a·√(π/12), 1024 modes, rotation angle (90° − azimuth). The orientation of the GSTools models was verified to coincide with the analytic one. Variograms were computed with the same exact estimator for both generators.'));
c.push(...fig('fig04_models_vs_gstools', 'Figure 6 – Mean experimental variogram (20 realizations) of SPECSIM3D, GSTools 1.7.0 and (spherical model) GSTools modified by the author vs the model, for three models and three directions. Shaded bands: 5th–95th percentile of the variograms of the 20 individual realizations.', 17));
const t5 = csv('fig04_models_rms');
const rows5 = []; for (const m of ['spherical', 'exponential', 'gaussian']) for (const d of [...new Set(t5.slice(1).map(r => r[2]))]) {
  const a = t5.find(r => r[0] === m && r[1] === 'SPECSIM3D' && r[2] === d), b = t5.find(r => r[0] === m && r[1] === 'GSTools' && r[2] === d), g = t5.find(r => r[0] === m && r[1] === 'GSTools (modified)' && r[2] === d); rows5.push([m, d, a[3], b ? b[3] : 'n/a', g ? g[3] : '–']); }
table([['Model', 'Direction', 'SPECSIM3D', 'GSTools 1.7.0', 'GSTools modified'], ...rows5], [1700, 2700, 1700, 1700, 1700], 'Table 4 – RMS error of the mean variogram within the range. GSTools modified: version modified by the author for the 3D spherical model (exact radial sampling of the wave numbers, clipping of the spectral density); not applicable to the other models (–).').forEach(x => c.push(x));
c.push(P('SPECSIM3D reproduces the three models in all directions with errors of the order of 10⁻³, compatible with the Monte Carlo noise of 20 realizations. GSTools reproduces the exponential and Gaussian models well (errors around 0.01–0.025) but not the spherical one in 3D (errors around 0.12–0.15 and a spurious nugget-like behaviour at short lags). In GSTools 1.7.0 the wave numbers are drawn by inversion only where an analytic quantile function exists (Gaussian and exponential models in 1D and 2D); in 3D all three models, and the spherical model in any dimension, are sampled by MCMC (emcee, 50 walkers and 20 burn-in steps). For the spherical model the cause of the bias was traced to the numerical Hankel transform of the covariance: at large wave numbers it returns small negative values, which spectral_rad_pdf turns positive with an absolute value, giving a radial density with a non-integrable 1/k tail; the chain then drifts to very large wave numbers, the more so the longer the burn-in. Drawing the radii from the exact spherical spectrum removes the bias (RMS error 0.013 instead of 0.137 in a test with range 10 cells). The small errors of the 3D exponential and Gaussian models are consistent with the short burn-in of the same sampler. The problem was reported to the GSTools developers.'));
c.push(H('6.2 Different ranges, anisotropy ratios and rotations', HeadingLevel.HEADING_2));
c.push(P('Six anisotropy cases (isotropic, 2:1, 4:1, layered 5:1, extreme 10:1:20, and a rotated case with azimuth 30°, dip 20° and plunge 15°) times three models, plus a nested spherical + exponential model with nugget, were tested with 20 realizations each. For each case the mean variogram was compared with the analytic one at up to 300 arbitrary integer lag vectors (not only along the axes), so the full 3D structure is tested. The ratio RMS/noise compares the error with the expected Monte Carlo standard error of the mean variogram: values near 1 mean agreement within noise.'));
c.push(...fig('fig05_anisotropy_cases', 'Figure 7 – Mean experimental vs model semivariogram at arbitrary lag vectors, all cases and models.', 9.5));
const t8 = csv('fig05_cases');
table([['Case', 'Model', 'Lags', 'RMS', 'MC noise', 'RMS/noise', 'Bias'], ...t8.slice(1).map(r => [r[0], r[1], r[2], r[3], r[4], r[5], r[6]])], [3000, 1400, 800, 1000, 1200, 1300, 1200], 'Table 5 – Anisotropy and range tests (SPECSIM3D).').forEach(x => c.push(x));
c.push(P('RMS/noise ranges from 0.6 to 1.7 (the highest values occur in cases with few informative lags), so the deviations are of the size of the Monte Carlo noise. The small positive bias (≈ +0.005) is explained by the fact that all cases of this particular test used the same random seed, so errors are correlated; a dedicated test with 400 realizations and independent seeds gave no significant bias (|t| < 1.2).'));
c.push(H('6.3 Gaussianity', HeadingLevel.HEADING_2));
c.push(P('Random-phase spectral simulation yields a field that is Gaussian only approximately: each value is a sum of many independent contributions, so by the central limit theorem the marginal distribution approaches the normal one, but it is not exactly Gaussian. We tested, on a 16×16×8 grid (exponential model, ranges 6/3/2), the marginal distribution of 400 randomly chosen nodes across 2000 independent realizations (600 for GSTools): Shapiro–Wilk, skewness and excess kurtosis.'));
c.push(...fig('fig06_gaussianity', 'Figure 8 – Normal QQ plot, distribution of the Shapiro–Wilk p-values (uniform if Gaussian) and of the excess kurtosis.', 17));
const t7 = csv('fig06_gaussianity');
table([['Method', 'Realizations', 'p < 0.05', 'Skewness', 'Excess kurtosis', 'Node variance'], ...t7.slice(1).map(r => [r[0], r[1], r[3], r[4], r[5], r[6]])], [1800, 1600, 1500, 1500, 1800, 1700], 'Table 6 – Marginal gaussianity tests (expected: p<0.05 fraction 0.05, skewness 0, kurtosis 0, variance 1).').forEach(x => c.push(x));
c.push(P('Both generators pass. Additional tests (Mardia kurtosis on groups of four nodes and bivariate tail probabilities, test_gauss.py) also gave results consistent with a multivariate normal distribution for SPECSIM3D.'));
c.push(H('6.4 Check against the original SPECSIM2', HeadingLevel.HEADING_2));
c.push(P('With nz = 1 SPECSIM3D reproduces the 2D case. The parameter file supplied with SPECSIM2 (exponential model, nugget 0.3, sill 0.7, 128×128) was run with both programs (ranges converted to the GSLIB convention). Since the random number sequences differ in how they are consumed, the comparison is statistical: the variograms of the SPECSIM2.EXE output and of 20 SPECSIM3D realizations are compared with the model.'));
c.push(...fig('fig08_validation_2d', 'Figure 9 – 2D check: SPECSIM3D (nz = 1) and the original SPECSIM2.EXE output vs the model. The band is the 5th–95th percentile of the variograms of the 20 SPECSIM3D realizations.', 13.5));
const t11 = csv('fig08_validation_2d');
table([['Direction', 'RMS SPECSIM3D', 'RMS SPECSIM2 (1 realization)'], ...t11.slice(1)], [3000, 3200, 3700], 'Table 7 – 2D validation errors.').forEach(x => c.push(x));
c.push(H('6.5 Example realization', HeadingLevel.HEADING_2));
c.push(...fig('fig03_slices', 'Figure 10 – One realization of the Yao case: horizontal slice and two vertical slices.', 14));
c.push(H('6.6 Comparison with gstlearn', HeadingLevel.HEADING_2));
c.push(P('gstlearn (version 1.11.1) provides a spectral simulation, `simuSpectral`, based on a sum of ns cosines with random frequencies. In this version it accepts the exponential and Gaussian models but refuses the spherical, cubic, stable and linear ones ("not valid for SPECTRAL"), so no result is given for the spherical model. The test of Section 6.1 was repeated with ns = 1024 and ns = 10000 (script `tests/compare_gstlearn.py`; the Python function `simuFFT` of gstlearn was not used). The times were measured on the desktop of Section 7.1 (Intel Core i7-12700).'));
c.push(...fig('gstlearn_comparison', 'Figure 11 – Mean experimental variogram of 20 realizations, SPECSIM3D and gstlearn simuSpectral (ns = 1024 and 10000), exponential and Gaussian models. Shaded bands: 5th–95th percentile of the variograms of the individual realizations.', 17));
const tg = fs.readFileSync(path.join(FIG, 'tables', 'gstlearn_comparison_metrics.csv'), 'utf8').trim().split(/\r?\n/).map(l => l.split(';'));
table([['Model', 'Method', 'Time (s)', 'Mean variance', 'RMS N45E', 'RMS N45W', 'RMS vertical'], ...tg.slice(1).map(r => [r[0], r[1], r[2], r[3], r[4], r[7], r[10]])], [1300, 3400, 1000, 1300, 1000, 1000, 1200], 'Table 8 – SPECSIM3D and gstlearn simuSpectral: time for 20 realizations (64×64×16 grid), mean variance and RMS error of the mean variogram within the range in the three principal directions.').forEach(x => c.push(x));
c.push(P('With 1024 components gstlearn gives errors of 0.02–0.08 and a mean variance of 0.93–0.94, against 0.001–0.004 for SPECSIM3D; with 10000 components the errors decrease to 0.008–0.027 and the time rises to about 120 s, against 0.15–0.18 s for SPECSIM3D. The error of a random-frequency spectral sum decreases with the number of components, whereas the FFT method samples all frequencies of the periodic grid.'));

c.push(H('7. Performance and memory', HeadingLevel.HEADING_1));
c.push(P(`**Benchmark machine.** The scaling measurements of this manual (Figures 12–13, Tables 3 and 9–11) were run on ${scal.system.cpu}, ${scal.system.cores} logical cores (OpenMP threads: ${scal.system.omp_threads}), ${scal.system.ram_gb} GB RAM, ${scal.system.platform}; compiler: ${scal.system.compiler}. Timings depend strongly on the computer, in particular on its memory system (Section 7.1). Run \`python make_figures.py --full\` or \`run_full_validation.bat\` on your machine to regenerate them with its own CPU, RAM and compiler (reported in figures/system_info.json).`));
const t9 = csv('fig09_scaling');
c.push(...fig('fig09_scaling', 'Figure 12 – Time per realization, peak memory and variogram error as a function of the grid size, SPECSIM3D vs GSTools (spherical model, 1024 modes).', 17));
table([['n', 'Periodic grid', 'Spectrum (s)', 's/realization', 'Peak RSS (MB)', 'RMS error', 'GSTools s/real.'], ...t9.slice(1).map(r => [r[0], r[1], r[2], r[3], r[4], r[5], r[6] || '–'])], [800, 1700, 1300, 1500, 1500, 1300, 1800], 'Table 9 – Scaling of SPECSIM3D (spherical model, 20/10/5).').forEach(x => c.push(x));
c.push(P('A 256³ field is generated in about 1.8 s per realization (20 threads) with a peak memory of about 3.3 GB, because the periodic grid is 512³ (complex double precision). Memory is the practical limit: with ipad = 1 expect roughly 8 × n³ × 16 bytes plus work arrays. GSTools needs about 38 s per 128³ realization with 1024 modes, i.e. about 140 times more than SPECSIM3D, and the cost per realization grows with the number of modes and nodes.'));
c.push(H('7.1 Dependence on the number of threads', HeadingLevel.HEADING_2));
const t15 = csv('fig15_thread_scaling');
c.push(...fig('fig15_thread_scaling', 'Figure 13 – Time per realization and speed-up of the realizations and of the spectrum step against the number of OpenMP threads (64³, 128³ and 256³ output grids).', 17));
const sp = (n, th) => { const r = t15.slice(1).find(x => x[0] == n && x[2] == th); return r ? r[5] : '–'; };
c.push(P(`The program runs the FFT of every realization with OpenMP threads (set the number with the environment variable OMP_NUM_THREADS). On this machine the speed-up of a 256³ realization against one thread is ${sp(256, 2)} with 2 threads, ${sp(256, 8)} with 8 and ${sp(256, 16)} with 16; for 128³ it is ${sp(128, 8)} with 8 threads and ${sp(128, 16)} with 16. The speed-up is small because the transform works on a periodic grid much larger than the processor cache, so that memory bandwidth, and not arithmetic, probably limits the speed (not verified by profiling). Independent runs with different seeds can be started as separate processes when memory allows.`));



c.push(H('8. Conditional simulation', HeadingLevel.HEADING_1));
c.push(P('With icond = 1 in the parameter file the realizations are conditioned on a set of data; the conditioning block is read from the same file, directly after the flag. The unconditional spectral realizations are generated exactly as before and are then corrected by **kriging of the residual** (conditioning by kriging; Journel and Huijbregts, 1978; Chilès and Delfiner, 2012):'));
c.push(...code('z_c(x) = z*_SK(x; data) + [ z_s(x) - z*_SK(x; z_s at the data nodes) ]'));
c.push(P('z_s is the unconditional realization, and z*_SK is the simple kriging estimate computed with the model covariance of the parameter file. The same kriging weights are applied to the data and to the simulated values at the data nodes, so that the conditioned field honours the data and has, node by node, the simple kriging mean and the simple kriging variance. The weights are computed once per node and used for all realizations; each additional realization costs one inverse FFT and a weighted sum of at most *ndmax* values per node.'));
c.push(P('Yao (1998) conditions spectral simulations through the identification of the phases of the spectrum; SPECSIM3D keeps the exact spectral generator and conditions it afterwards by kriging of the residual, the approach also used by `CondSRF` of GSTools.'));
c.push(H('8.1 Conditioning block', HeadingLevel.HEADING_2));
c.push(P('The nine lines below follow the line icond = 1 of the parameter file. Because the block is part of the parameter file, the data file, the transform and the search settings used for a simulation are always recorded together with the model, the seed and the outputs.'));
c.push(...code(`datafile                  conditioning data, GSLIB format (title, nvar, nvar names, rows)
ixcol iycol izcol ivcol iwcol    columns of x, y, z, value, weight (0 = not used; izcol = 0 for 2D grids)
tmin tmax                 trimming limits on the value
itrans skmean             itrans = 0: values used as they are (simple kriging mean skmean);
                          itrans = 1: normal score transform + back transform
zmin zmax                 limits of the back transform (itrans = 1)
ltail ltpar               lower tail: 1 = linear, 2 = power with exponent ltpar
utail utpar               upper tail: 1 = linear, 2 = power, 4 = hyperbolic (exponent utpar)
radius_max radius_med radius_vert   search ellipsoid radii
sang1 sang2 sang3         search ellipsoid angles (azimuth, dip, plunge, same convention as the model)
ndmax                     maximum number of data per node (1 to 100)`));
c.push(P('Every line may be followed by a free comment. The data are assigned to the nearest grid node (as in GSLIB sgsim); data falling in the same node are averaged (weighted by iwcol), and data outside the grid are ignored. Therefore the realizations reproduce the data exactly at the nodes. Nodes with no data inside the search ellipsoid are not conditioned (they keep the unconditional value, which is the correct result of simple kriging with no data). **The variogram model in the main parameter file must be the model of the values used for conditioning**: the Gaussian values when itrans = 0, the normal scores (total sill 1) when itrans = 1.'));
c.push(H('8.2 Normal score transform and back transform', HeadingLevel.HEADING_2));
c.push(P('Spectral simulation produces Gaussian fields, so conditioning data that are not Gaussian must be transformed first. With itrans = 1 the program implements the GSLIB procedures: *nscore* ranks the data (optionally with declustering weights, column iwcol) and assigns the standard normal quantile of the mid-point cumulative probability; the realizations are conditioned in Gaussian space and every node is then back transformed with the rules of *backtr*: linear interpolation of the table between data, and the options ltail/utail with limits zmin/zmax for the tails (linear, power, hyperbolic). The transformation table is written to `<base>_nscore.trn` and the data actually used (with value and normal score) to `<base>_conddata.dat`. When the data are already Gaussian, set itrans = 0 (and skmean to their mean if it is not 0).'));
c.push(H('8.3 Neighbourhood search', HeadingLevel.HEADING_2));
c.push(P('For each node the program searches the data inside the search ellipsoid and keeps the *ndmax* closest ones in the scaled distance. To make the search independent of the number of data, the data are sorted into **super blocks** (as in the GSLIB superblock search): a regular block grid whose cell is half of the extent of the search ellipsoid; for a node only the blocks that the ellipsoid can reach are visited. The kriging system (at most ndmax × ndmax) is solved by Cholesky factorization, with a small regularization if the covariance matrix is nearly singular (a warning is printed). The loop over the grid nodes is parallelized with OpenMP. **Default values**: radii equal to the variogram ranges and ndmax = 16–24 give conditioning in a few seconds for grids of a few million nodes (Figure 17); larger values cost much more time (Figure 17b) and are rarely needed.'));
c.push(...code('specsim3d model.par'));
c.push(P('Ready-to-run examples (Gaussian data in 3D, lognormal data with normal score transform, a 2D case with nugget and off-grid data) are in the folder `par/conditional`, with a `README_conditional.txt` and the scripts `run_conditional_examples.bat`/`.sh`.'));
c.push(H('8.4 Validation', HeadingLevel.HEADING_2));
c.push(P('The test uses an exponential model (practical ranges 20/10/5, N45E, 64×64×16 nodes) and 150 data taken from an unconditional realization. The reference is the exact simple kriging with all the data, computed in NumPy. With 200 realizations the mean and the variance of the conditioned fields must differ from the kriging estimate and kriging variance only by the Monte Carlo error. SPECSIM3D was run with 24 neighbours, GSTools `CondSRF` (1.7.0, 1024 modes, simple kriging) uses all the data.'));
c.push(...fig('fig11_conditional_maps', 'Figure 14 – Conditional simulation: exact simple kriging, mean of 200 realizations, one realization and standard deviation of the realizations in the layer z = 8.', 17));
const t12 = csv('fig12_conditional_validation');
c.push(...fig('fig12_conditional_validation', 'Figure 15 – Mean and variance of the conditioned realizations against the exact simple kriging estimate and variance, SPECSIM3D and GSTools CondSRF.', 13));
table([['Code', 'Realizations', 'RMS mean diff.', 'MC error', 'RMS var. diff.', 'MC error', 'Var. slope', 'Max |data error|', 's/realization'], ...t12.slice(1).map(r => [r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[8], r[9]])], [1500, 1000, 1000, 900, 1000, 900, 900, 1300, 1200], 'Table 10 – Conditional realizations against exact simple kriging (MC error = expected difference from Monte Carlo sampling).').forEach(x => c.push(x));
c.push(P('Both codes agree with simple kriging within the Monte Carlo error and reproduce the data. The data error of SPECSIM3D is at the level of the float32 precision of the .bin output. The use of a limited neighbourhood (24 data) in SPECSIM3D does not produce a visible difference from the exact kriging in this test.'));
c.push(H('8.5 Normal score example', HeadingLevel.HEADING_2));
const t13 = csv('fig13_conditional_nscore');
c.push(P('Lognormal data (150 values taken from exp(1 + 0.8 g), g a Gaussian realization with the same model) are transformed to normal scores, the realizations are conditioned and back transformed (linear tails, zmin = 0). The cumulative distribution of the conditioned realizations follows that of the data and of the truth (Figure 16); the mean of the realizations (' + t13[4][1] + ') is that of the data (' + t13[3][1] + '), and the data are reproduced with a relative error of ' + t13[1][1] + ' (float32 output).'));
c.push(...fig('fig13_conditional_nscore', 'Figure 16 – Normal score transform: marginal distribution of truth, data and conditioned realizations, truth and one realization in the layer z = 8.', 17));
c.push(H('8.6 Performance', HeadingLevel.HEADING_2));
const t14 = csv('fig14_conditional_speed_grid');
c.push(...fig('fig14_conditional_speed', 'Figure 17 – Time of the conditional simulation as a function of the grid size (10 realizations, 16 neighbours) and of the maximum number of neighbours.', 16));
table([['Nodes', 'Data', 'Unconditional (s)', 'Conditional total (s)', 'Conditioning (s)', 'GSTools CondSRF (s)'], ...t14.slice(1).map(r => [r[0], r[1], r[2], r[3], r[4], r[5] === 'nan' ? '–' : r[5]])], [1400, 1000, 1700, 1900, 1700, 2000], 'Table 11 – Time for 10 realizations (' + sys.cpu + ', ' + sys.cores + ' threads, ' + sys.ram_gb + ' GB RAM, ' + sys.compiler + '; GSTools extrapolated from 3 realizations).').forEach(x => c.push(x));
c.push(P('On the two largest grids conditioning adds about 80 % to the time of the unconditional simulation when 16 neighbours are used; the times of the two smallest grids are dominated by the start-up of the program. GSTools CondSRF was run only on the two smaller grids, where its time for 10 realizations (extrapolated from 3) is one to two orders of magnitude larger.'));
c.push(H('8.7 Limitations of the conditional mode', HeadingLevel.HEADING_2));
['Simple kriging with known mean (0 or skmean): the mean is not estimated. Ordinary kriging and trend models are not implemented.',
 'Data are moved to the nearest grid node: the conditioning is exact at the nodes, not at the true data positions. Use a grid spacing small compared with the variogram range.',
 'The neighbourhood is limited (ndmax data in the search ellipsoid), so the covariance of the conditioned fields is approximated where the data are sparse. No octant search is available.',
 'Data errors (a nugget on the data) are not modelled: the nugget of the model is part of the covariance, and the data are honoured exactly.',
 'The normal score transform assumes a stationary univariate distribution; declustering weights must be supplied by the user.'].forEach(t => c.push(B(t)));

c.push(H('9. Limitations and known issues', HeadingLevel.HEADING_1));
['The program simulates stationary Gaussian fields, unconditional or conditioned on data by simple kriging of the residual (Section 8); no non-stationarity.',
 'The covariance must be positive definite on the periodic grid. The models implemented are valid in 3D; very small anisotropic ranges relative to the cell size may yield small negative spectral values, which are clipped and reported.',
 'Memory: with ipad = 1 the periodic grid may be up to 8× the output grid. For very large models use the compact grid (ipad = 0) only if periodicity is acceptable, or reduce the grid.',
 'The sign convention of the dip should be compared with the user\'s own software before use with dipping structures.',
 'The reference for validation is the analytic model; the original SPECSIM2 is 2D only, so for 3D the comparison is made with GSTools (version 1.7.0) and gstlearn (version 1.11.1), which do not reproduce all the models (Sections 6.1 and 6.6).',
 'Timings depend on the computer: all the measurements in this manual come from a desktop with an Intel Core i7-12700 (12 cores, 20 threads, 31.7 GB RAM, Windows 11).',
 'The program has not (yet) been compared with other 3D simulators (e.g. GSLIB sgsim, turning bands) on identical problems.'].forEach(t => c.push(B(t)));

c.push(H('10. Reproducing the figures', HeadingLevel.HEADING_1));
c.push(...code(`cd docs
python make_figures.py              # all figures + tables (needs numpy, scipy, matplotlib; gstools optional)
python make_figures.py --no-gstools # skip GSTools comparisons
python make_figures.py --full       # re-measure scaling (up to 256^3, ~3 GB RAM) on this machine
python make_figures.py --only figS1 fig07 fig12   # fig11-fig14: conditional simulation
node build_manual.js                # rebuild SPECSIM3D_Manual.docx (npm package "docx")`));
c.push(P('Figures are written to docs/figures (PNG, 300 dpi, and vector PDF) and the numerical tables to docs/figures/tables (CSV). The script finds the executable in the parent folder; on Windows use `make_figures.bat`.'));
c.push(H('10.1 Complete validation in one step', HeadingLevel.HEADING_2));
c.push(P('`run_full_validation.bat` (or `python run_full_validation.py`) performs the whole procedure automatically and stores everything in `validation/<date_time>/`. Run it from the Intel oneAPI command prompt, so that the ifx version is recorded, after `build_specsim3d.bat` and `setup_python_env.bat`, with the other heavy programs closed because the timings are part of the results. The steps are:'));
[ '**machine description**: CPU model, physical and logical cores, RAM, operating system, Fortran compiler version (ifx, ifort or gfortran), build date of the executable and versions of Python, numpy, scipy, matplotlib and GSTools, saved in `system_info.json` and in the summary;',
  '**tests that do not need GSTools**: `test_conditional`, `test_rotation`, `test_bias`;',
  '**tests against GSTools, run twice in separate processes**: with the released GSTools exactly as installed and with the author-modified GSTools (spherical 3D patch in `docs/gstools_patch`): `test_gauss`, `test_cases`, `compare_gstools`; the csv, png and json outputs of the two variants are kept in separate folders;',
  '**figures and benchmark**: `make_figures.py --full` (all figures at 300 dpi, tables, timings up to 256^3), or without `--full` with the option `--quick`;',
  '**report**: `summary.md` and `summary.json` with PASS/FAIL, time of every step and the machine description to be quoted in the Performance section of the paper and of this manual.'].forEach(t => c.push(B(t)));
c.push(P('Options: `--quick` (reduced benchmark), `--skip-tests`, `--skip-gstools`, `--skip-figures`. The full run takes from one to a few hours depending on the machine, mostly for the Gaussianity test with GSTools, which is executed twice. `run_all_tests.bat` runs only the quick tests (`test_scaling.py` is Linux-only: on Windows the scaling benchmark is part of `make_figures.py --full`). In the figures that compare variograms the line is the mean of the realizations and the shaded band is the 5th–95th percentile of the variograms of the individual realizations, so that the spread of a single realization around the model can be judged.'));

c.push(H('References', HeadingLevel.HEADING_1));
['Dietrich, C.R., Newsam, G.N. (1993). A fast and exact method for multidimensional Gaussian stochastic simulations. Water Resources Research 29(8), 2861–2869.',
 'Deutsch, C.V., Journel, A.G. (1998). GSLIB: Geostatistical Software Library and User\'s Guide, 2nd ed. Oxford University Press.',
 'Chilès, J.-P., Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty, 2nd ed. Wiley.',
 'Journel, A.G., Huijbregts, C.J. (1978). Mining Geostatistics. Academic Press.',
 'Le Ravalec, M., Noetinger, B., Hu, L.Y. (2000). The FFT moving average (FFT-MA) generator: an efficient numerical method for generating and conditioning Gaussian simulations. Mathematical Geology 32(6), 701–723.',
 'Müller, S., Schüler, L., Zech, A., Heße, F. (2022). GSTools v1.3: a toolbox for geostatistical modelling in Python. Geoscientific Model Development 15, 3161–3182.',
 'Pardo-Igúzquiza, E., Chica-Olmo, M. (1993). The Fourier integral method: an efficient spectral method for simulation of random fields. Mathematical Geology 25(2), 177–217.',
 'Pardo-Igúzquiza, E., Dowd, P.A. (2003). SPECSIM2: a program for spectral simulation of anisotropic two-dimensional random fields. Computers & Geosciences 29, 897–907.',
 'Yao, T. (1998). Conditional spectral simulation with phase identification. Mathematical Geology 30(3), 285–308.'].forEach(t => c.push(P(t, { indent: { left: 360, hanging: 360 } })));

const doc = new Document({
  creator: 'Stefano Bonduà', title: 'SPECSIM3D Manual',
  styles: { default: { document: { run: { font: 'Calibri', size: 21 } } },
    paragraphStyles: [
      { id: 'Heading1', name: 'Heading 1', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 32, bold: true, color: '1F3864' }, paragraph: { spacing: { before: 320, after: 140 }, outlineLevel: 0 } },
      { id: 'Heading2', name: 'Heading 2', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { size: 25, bold: true, color: '2E5597' }, paragraph: { spacing: { before: 220, after: 100 }, outlineLevel: 1 } }] },
  numbering: { config: [{ reference: 'b', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 270 } } } }] }] },
  sections: [{ properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1134, bottom: 1134, left: 1003, right: 1003 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ text: 'SPECSIM3D manual – Stefano Bonduà – page ', size: 17 }), new TextRun({ children: [PageNumber.CURRENT], size: 17 })] })] }) },
    children: c }]
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(path.join(__dirname, 'SPECSIM3D_Manual.docx'), b); console.log('manual written'); });
