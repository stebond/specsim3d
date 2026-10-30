"""Scaling of SPECSIM3D (Fortran, OpenMP, 2 threads) vs GSTools (randomization, 1024 modes) on cubic grids.

Spherical model 20/10/5, azimuth 45, variance 1, unit spacing. Time, peak memory (RSS of the
process) and accuracy (RMS of the experimental variogram along the three principal axes within
the range, with respect to the model) are measured. GSTools is run only up to 128^3: its cost is
proportional to (nodes x modes) and at 256^3 it would take hours.
Autore: Stefano Bonduà
"""
import subprocess, re, time, json, sys
import numpy as np

import os as _os
WD = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.path.dirname(WD)
EXE = _os.path.join(ROOT, "specsim3d.exe" if _os.name == "nt" else "specsim3d")
WORK = _os.path.join(WD, "work")      # generated .par/.bin/.out files
RES = _os.path.join(WD, "results")    # csv/png/json results
_os.makedirs(WORK, exist_ok=True); _os.makedirs(RES, exist_ok=True)
SIZES = [32, 64, 128, 256]
RANGES = (20.0, 10.0, 5.0)


def corr_sph(h, a):
    r = np.asarray(h, float) / a
    return np.where(r < 1, 1 - 1.5 * r + 0.5 * r**3, 0.0)


def vario_shift(f, step, kmax):
    out, lags = [], []
    for k in range(1, kmax + 1):
        d = [k * s for s in step]
        if any(abs(di) >= n for di, n in zip(d, f.shape)):
            break
        sa = tuple(slice(max(0, -di), n - max(0, di)) for di, n in zip(d, f.shape))
        sb = tuple(slice(max(0, di), n - max(0, -di)) for di, n in zip(d, f.shape))
        out.append(0.5 * np.mean((f[sa] - f[sb]) ** 2, dtype=np.float64))
        lags.append(k * np.linalg.norm(step))
    return np.array(lags), np.array(out)


def acc(f):
    """RMS (3 axes) of the variogram within the range; f is a single field (nx,ny,nz)."""
    res = []
    for step, a in [((1, 1, 0), RANGES[0]), ((1, -1, 0), RANGES[1]), ((0, 0, 1), RANGES[2])]:
        l, g = vario_shift(f, step, int(a * 1.0))
        m = l <= a
        res.append(np.sqrt(np.mean((g[m] - (1 - corr_sph(l[m], a))) ** 2)))
    return float(np.mean(res))


out = {}
for n in SIZES:
    nreal = 1 if n >= 256 else 2
    tag = f"scal_{n}"
    open(f"{WORK}/{tag}.par", "w").write(f"0 0 0\n1 1 1\n{n} {n} {n}\n0.0 1\n1 1.0 {RANGES[0]} {RANGES[1]} {RANGES[2]} 45 0 0\n99 {nreal} 1\n{tag}\n")
    t0 = time.time()
    wrapper = ("import subprocess,resource,sys;p=subprocess.run(sys.argv[1:],capture_output=True,text=True);"
               "print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss);sys.stdout.write(p.stdout)")
    p = subprocess.run([sys.executable, "-I", "-c", wrapper, EXE, f"{WORK}/{tag}.par", "--bin"],
                       cwd=WORK, capture_output=True, text=True)
    wall = time.time() - t0
    first, _, rest = p.stdout.partition("\n")
    rss = int(first) / 1024.0
    p.stdout = rest
    grid = re.search(r"periodic grid:\s+(\d+)\s+(\d+)\s+(\d+)", p.stdout).groups()
    t_spec = float(re.search(r"covariance\+spectrum:\s+([\d.]+)", p.stdout).group(1))
    t_sim = float(re.search(r"realizations\s+:\s+([\d.]+)", p.stdout).group(1))
    a = np.fromfile(f"{WORK}/{tag}.bin", dtype=np.float32, count=n**3).reshape(n, n, n)   # memory order (z,y,x), transposed below
    f = a.reshape(n, n, n).transpose(2, 1, 0)                                               # (x,y,z) from (z,y,x)
    out[f"specsim3d_{n}"] = dict(n=n, periodic=list(map(int, grid)), wall_s=wall, spectrum_s=t_spec,
                                 per_realization_s=t_sim / nreal, rss_mb=rss, rms=acc(f), nreal=nreal)
    print("SPECSIM3D", out[f"specsim3d_{n}"], flush=True)
    subprocess.run(["rm", "-f", f"{WORK}/{tag}.bin"])

gs_code = r'''
import sys, time, json, resource, numpy as np, gstools as gs
n=int(sys.argv[1]); nreal=int(sys.argv[2])
mod=gs.Spherical(dim=3,var=1.0,len_scale=[20.,10.,5.],angles=[np.pi/4,0,0])
x=y=np.arange(n,dtype=float); z=np.arange(n,dtype=float)
t0=time.time(); srf=gs.SRF(mod,seed=99,mode_no=1024); f=srf.structured((x,y,z)); t1=time.time()
for i in range(1,nreal):
    f2=gs.SRF(mod,seed=99+i,mode_no=1024).structured((x,y,z))
t2=time.time()
np.save(sys.argv[3], f.astype(np.float32))
print(json.dumps(dict(per_realization_s=(t2-t0)/nreal, rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0)))
'''
open(f"{WORK}/_gs_scal.py", "w").write(gs_code)
for n in [32, 64, 128]:
    nreal = 2
    p = subprocess.run([sys.executable, "-I", f"{WORK}/_gs_scal.py", str(n), str(nreal), f"{WORK}/scal_gs_{n}.npy"], capture_output=True, text=True)
    r = json.loads(p.stdout.strip().splitlines()[-1])
    f = np.load(f"{WORK}/scal_gs_{n}.npy").astype(float)
    r.update(n=n, rms=acc(f), nreal=nreal)
    out[f"gstools_{n}"] = r
    print("GSTools  ", r, flush=True)
json.dump(out, open(f"{RES}/scaling.json", "w"), indent=1)
