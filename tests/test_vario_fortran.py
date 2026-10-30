"""Check of the Fortran variogram module (--vario) against a NumPy reference.

Five realizations of the Yao test case (64 x 64 x 16, spherical, ranges 20/10/5, N45E) are written as float32
(.bin) and the mean experimental variogram of the three directions (1,1,0), (1,-1,0), (0,0,1) computed by
the program is compared with the one computed by NumPy on the same .bin file. The expected maximum difference
is the float32 rounding of the file (about 5e-7). Result: results/vario_fortran_vs_numpy.csv.
Autore: Stefano Bonduà
"""
import csv
import os
import subprocess
import sys

import numpy as np

WD = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WD)
EXE = os.path.join(ROOT, "specsim3d.exe" if os.name == "nt" else "specsim3d")
WORK = os.path.join(WD, "work")
RES = os.path.join(WD, "results")
os.makedirs(WORK, exist_ok=True); os.makedirs(RES, exist_ok=True)
NX, NY, NZ, NREAL, KMAX = 64, 64, 16, 5, 40
DIRS = [(1, 1, 0), (1, -1, 0), (0, 0, 1)]
TOL = 2e-6


def vario_shift(f, step, kmax):
    """Mean NumPy variogram along k*step (f: nreal,x,y,z)."""
    out = []
    for k in range(1, kmax + 1):
        d = [k * s for s in step]
        if any(abs(di) >= n for di, n in zip(d, f.shape[1:])):
            break
        sa = [slice(None)] + [slice(max(0, -di), n - max(0, di)) for di, n in zip(d, f.shape[1:])]
        sb = [slice(None)] + [slice(max(0, di), n - max(0, -di)) for di, n in zip(d, f.shape[1:])]
        out.append(0.5 * np.mean((f[tuple(sa)] - f[tuple(sb)]) ** 2, axis=(1, 2, 3)))
    return np.array(out).T.mean(axis=0)


def main():
    if not os.path.exists(EXE):
        sys.exit("Executable not found: " + EXE)
    base = "t_vario"
    with open(os.path.join(WORK, base + ".par"), "w") as fh:
        fh.write("Parameters for SPECSIM3D\n************************\n\nSTART OF PARAMETERS:\n"
                 f"0 0 0\n1 1 1\n{NX} {NY} {NZ}\n0.0 1\n1 1.0 20 10 5 45 0 0\n31791191 {NREAL} 1\n{base}\n1\n3\n0\n")
    with open(os.path.join(WORK, base + "_v.par"), "w") as fh:
        fh.write(f"{len(DIRS)} 0\n" + "\n".join(f"{a} {b} {c} {KMAX}" for a, b, c in DIRS) + "\n")
    r = subprocess.run([EXE, base + ".par", "--vario", base + "_v.par"], cwd=WORK, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stdout + r.stderr)
    f = np.fromfile(os.path.join(WORK, base + ".bin"), dtype=np.float32).reshape(NREAL, NZ, NY, NX).transpose(0, 3, 2, 1).astype(np.float64)
    rows = list(csv.reader(open(os.path.join(WORK, base + "_vario.csv"))))[1:]
    out, worst = [], 0.0
    for d, step in enumerate(DIRS, 1):
        ref = np.array([float(x[8]) for x in rows if int(x[0]) == d])
        g = vario_shift(f, step, KMAX)[:len(ref)]
        diff = np.abs(g - ref)
        worst = max(worst, diff.max())
        out.append([d, *step, len(ref), f"{diff.max():.3e}"])
        print(f"direction {step}: {len(ref)} lags, max |Fortran - NumPy| = {diff.max():.3e}")
    with open(os.path.join(RES, "vario_fortran_vs_numpy.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["dir", "sx", "sy", "sz", "nlags", "max_abs_difference"]); w.writerows(out)
    print(f"maximum difference over all directions and lags: {worst:.3e} (tolerance {TOL:g})")
    for ext in (".bin", "_vario.csv"):
        os.remove(os.path.join(WORK, base + ext))
    if worst > TOL:
        sys.exit("FAILED: difference above tolerance")
    print("PASSED")


if __name__ == "__main__":
    main()
