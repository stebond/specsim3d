"""Check of the 3D angle convention (azimuth, dip, plunge) of SPECSIM3D.

The principal axes are chosen so that the directions are aligned with integer grid vectors:
azimuth 45, dip 35.2644 (tan = 1/sqrt2) -> major axis along (1,1,+-1).
For every integer direction the mean variogram (20 realizations, spherical) is compared with
the theoretical one for the ranges 20, 10, 5 and the range that describes it best is listed.
Autore: Stefano Bonduà
"""
import subprocess
import numpy as np

import os as _os
WD = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.path.dirname(WD)
EXE = _os.path.join(ROOT, "specsim3d.exe" if _os.name == "nt" else "specsim3d")
WORK = _os.path.join(WD, "work")      # generated .par/.bin/.out files
RES = _os.path.join(WD, "results")    # csv/png/json results
_os.makedirs(WORK, exist_ok=True); _os.makedirs(RES, exist_ok=True)
N = 64
NREAL = 20


def corr_sph(h, a):
    r = np.asarray(h, float) / a
    return np.where(r < 1, 1 - 1.5 * r + 0.5 * r**3, 0.0)


def vario_shift(f, step, kmax=60):
    out, lags = [], []
    for k in range(1, kmax + 1):
        d = [k * s for s in step]
        if any(abs(di) >= n for di, n in zip(d, f.shape)):
            break
        sa = tuple(slice(max(0, -di), n - max(0, di)) for di, n in zip(d, f.shape))
        sb = tuple(slice(max(0, di), n - max(0, -di)) for di, n in zip(d, f.shape))
        out.append(0.5 * np.mean((f[sa] - f[sb]) ** 2))
        lags.append(k * np.linalg.norm(step))
    return np.array(lags), np.array(out)


def run(az, dip, plg, tag):
    open(f"{WORK}/{tag}.par", "w").write(
        f"0 0 0\n1 1 1\n{N} {N} {N}\n0.0 1\n1 1.0 20 10 5 {az} {dip} {plg}\n777 {NREAL} 1\n{tag}\n")
    subprocess.run([EXE, f"{WORK}/{tag}.par"], cwd=WORK, check=True, capture_output=True)
    raw = np.loadtxt(f"{WORK}/{tag}.out", skiprows=2 + NREAL)
    return np.array([raw[:, r].reshape(N, N, N).transpose(2, 1, 0) for r in range(NREAL)])


dips = np.degrees(np.arctan(1 / np.sqrt(2)))
tests = [("az45 dip35.26 plg0", 45.0, dips, 0.0), ("az45 dip35.26 plg90", 45.0, dips, 90.0),
         ("az45 dip-35.26 plg0", 45.0, -dips, 0.0)]
dirs = {"(1,1,1)": (1, 1, 1), "(1,1,-1)": (1, 1, -1), "(1,-1,0)": (1, -1, 0), "(1,1,2)": (1, 1, 2), "(1,1,-2)": (1, 1, -2)}
for name, az, dip, plg in tests:
    fields = run(az, dip, plg, "rot_" + name.split()[0] + str(int(plg)) + ("m" if dip < 0 else "p"))
    print(f"\n=== {name}")
    for dname, step in dirs.items():
        gs_ = [vario_shift(f, step, 40) for f in fields]
        l = gs_[0][0]
        g = np.mean([x[1] for x in gs_], 0)
        best = None
        errs = {}
        for a in (20.0, 10.0, 5.0):
            m = l <= max(a * 1.3, 4)
            errs[a] = np.sqrt(np.mean((g[m] - (1 - corr_sph(l[m], a))) ** 2))
        best = min(errs, key=errs.get)
        print(f"  {dname:9s} RMS vs a=20: {errs[20.0]:.3f}  a=10: {errs[10.0]:.3f}  a=5: {errs[5.0]:.3f}  -> best range = {best:.0f}")
