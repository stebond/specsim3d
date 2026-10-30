"""Bias of the mean variogram of SPECSIM3D with many realizations (different seeds between runs).
For every realization the mean error (experimental variogram - model) over the lags within the
range in the three principal directions is computed; bias = mean over NREAL realizations, SE = sd/sqrt(NREAL).
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
NX, NY, NZ, NREAL = 64, 64, 16, 400
R = (20.0, 10.0, 5.0)


def corr(it, r):
    r = np.asarray(r, float)
    return np.where(r < 1, 1 - 1.5 * r + 0.5 * r**3, 0.0) if it == 1 else (np.exp(-3 * r) if it == 2 else np.exp(-3 * r**2))


def vario(f, step, kmax):                       # f: (nreal,x,y,z) -> (nreal,k), lags
    out, lags = [], []
    for k in range(1, kmax + 1):
        d = [k * s for s in step]
        sa = [slice(None)] + [slice(max(0, -di), n - max(0, di)) for di, n in zip(d, f.shape[1:])]
        sb = [slice(None)] + [slice(max(0, di), n - max(0, -di)) for di, n in zip(d, f.shape[1:])]
        out.append(0.5 * np.mean((f[tuple(sa)] - f[tuple(sb)]) ** 2, axis=(1, 2, 3)))
        lags.append(k * np.linalg.norm(step))
    return np.array(out).T, np.array(lags)


print(f"{'model':13s} {'ipad':>4s} {'bias':>9s} {'SE':>8s} {'t':>6s}   (mean error over lags <= range, 3 directions, {NREAL} realizations)")
for model, it in [("spherical", 1), ("exponential", 2), ("gaussian", 3)]:
    for ipad in (1, 0):
        tag = f"bias_{model}_{ipad}"
        open(f"{WORK}/{tag}.par", "w").write(f"0 0 0\n1 1 1\n{NX} {NY} {NZ}\n0.0 1\n{it} 1.0 {R[0]} {R[1]} {R[2]} 45 0 0\n{1000+it*10+ipad} {NREAL} {ipad}\n{tag}\n")
        subprocess.run([EXE, f"{WORK}/{tag}.par", "--bin"], cwd=WORK, check=True, capture_output=True)
        f = np.fromfile(f"{WORK}/{tag}.bin", dtype=np.float32).reshape(NREAL, NZ, NY, NX).astype(float).transpose(0, 3, 2, 1)
        errs = []
        for step, a in [((1, 1, 0), R[0]), ((1, -1, 0), R[1]), ((0, 0, 1), R[2])]:
            g, l = vario(f, step, int(a))
            m = l <= a
            errs.append((g[:, m] - (1 - corr(it, l[m] / a))).mean(1))
        e = np.mean(errs, axis=0)
        print(f"{model:13s} {ipad:4d} {e.mean():+9.4f} {e.std(ddof=1)/np.sqrt(NREAL):8.4f} {e.mean()/(e.std(ddof=1)/np.sqrt(NREAL)):6.2f}", flush=True)
