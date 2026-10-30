"""Comparison of SPECSIM3D (spectral, FFT) with gstlearn simuSpectral (spectral method with random frequencies).

Same test case as compare_gstools.py (Yao, 1998): 64x64x16 grid, unit spacing, unit variance, practical ranges
20/10/5, major direction N45E, 20 realizations, exact experimental variograms along the three principal
directions, compared with the analytic variogram.
Only the exponential and the Gaussian models are compared: in gstlearn 1.11.1 simuSpectral refuses the
spherical (and cubic, stable, linear) covariance ("not valid for SPECTRAL"), so no result is reported for it.
simuSpectral is run with ns = 1024 and 10000 spectral components (ns is the number of random frequencies,
comparable to the 1024 modes of GSTools).
Needs:  pip install gstlearn      (tested with 1.11.1)
Autore: Stefano Bonduà
"""
import subprocess, time, json
import numpy as np
import gstlearn as gl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os as _os
WD = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.path.dirname(WD)
EXE = _os.path.join(ROOT, "specsim3d.exe" if _os.name == "nt" else "specsim3d")
WORK = _os.path.join(WD, "work")
RES = _os.path.join(WD, "results")
_os.makedirs(WORK, exist_ok=True); _os.makedirs(RES, exist_ok=True)
NX, NY, NZ = 64, 64, 16
NREAL = 20
SEED = 31791191
RANGES = (20.0, 10.0, 5.0)
MODELS = {"exponential": (2, gl.ECov.EXPONENTIAL), "gaussian": (3, gl.ECov.GAUSSIAN)}   # spherical: not supported by simuSpectral
NSPEC = (1024, 10000)                          # number of spectral components of simuSpectral


def corr(itype, h, a):
    r = np.asarray(h, float) / a
    if itype == 1:
        return np.where(r < 1, 1 - 1.5 * r + 0.5 * r**3, 0.0)
    if itype == 2:
        return np.exp(-3 * r)
    return np.exp(-3 * r**2)


def gamma_model(itype, h, a):
    return 1.0 - corr(itype, h, a)


def vario_shift(f, step, kmax):
    sx, sy, sz = step
    out, lags = [], []
    for k in range(1, kmax + 1):
        d = (k * sx, k * sy, k * sz)
        if any(abs(di) >= n for di, n in zip(d, f.shape)):
            break
        sa = tuple(slice(max(0, -di), n - max(0, di)) for di, n in zip(d, f.shape))
        sb = tuple(slice(max(0, di), n - max(0, -di)) for di, n in zip(d, f.shape))
        out.append(0.5 * np.mean((f[sa] - f[sb]) ** 2))
        lags.append(k * np.sqrt(sx**2 + sy**2 + sz**2))
    return np.array(lags), np.array(out)


DIRS = [("N45E (a=20)", (1, 1, 0), 0, 40), ("N45W (a=10)", (1, -1, 0), 1, 40), ("Vertical (a=5)", (0, 0, 1), 2, 10)]


def analyse(fields):
    res = []
    for name, step, ax, kmax in DIRS:
        gs_ = []
        for f in fields:
            l, g = vario_shift(f, step, kmax)
            gs_.append(g)
        res.append((l, np.mean(gs_, 0), np.array(gs_)))
    return res


def metrics(res, itype):
    rows = []
    for (name, step, ax, kmax), (l, gm, gall) in zip(DIRS, res):
        a = RANGES[ax]
        m = l <= a
        th = gamma_model(itype, l, a)
        rows.append((np.sqrt(np.mean((gm[m] - th[m]) ** 2)),
                     np.mean(np.sqrt(np.mean((gall[:, m] - th[m]) ** 2, axis=1))), gm[0] - th[0]))
    return rows


def run_specsim3d(itype, grid, nreal, ipad=1, tag="cl"):
    nx, ny, nz = grid
    with open(f"{WORK}/{tag}.par", "w") as fh:
        fh.write(f"0.5 0.5 0.5\n1 1 1\n{nx} {ny} {nz}\n0.0 1\n"
                 f"{itype} 1.0 {RANGES[0]} {RANGES[1]} {RANGES[2]} 45.0 0.0 0.0\n{SEED} {nreal} {ipad}\n{tag}\n")
    t0 = time.time()
    subprocess.run([EXE, f"{WORK}/{tag}.par", "--bin"], cwd=WORK, check=True, capture_output=True)
    dt = time.time() - t0
    a = np.fromfile(f"{WORK}/{tag}.bin", dtype=np.float32).reshape(nreal, nz, ny, nx).astype(float)
    return a.transpose(0, 3, 2, 1), dt


def run_gstlearn(ctype, grid, nreal, ns):
    nx, ny, nz = grid
    db = gl.DbGrid.create(nx=[nx, ny, nz], dx=[1.0, 1.0, 1.0], x0=[0.0, 0.0, 0.0])
    model = gl.Model.createFromParam(ctype, range=RANGES[0], sill=1.0, ranges=list(RANGES), angles=[45.0, 0.0, 0.0], flagRange=True)
    t0 = time.time()
    err = gl.simuSpectral(None, db, model, nbsimu=nreal, seed=SEED, ns=ns)
    dt = time.time() - t0
    assert err == 0, "gstlearn simuSpectral failed"
    names = db.getAllNames()[-nreal:]
    fields = np.array([np.array(db.getColumn(n)).reshape(nz, ny, nx).transpose(2, 1, 0) for n in names])
    return fields, dt


methods = {}
labels = ["SPECSIM3D padded (ipad=1)"] + [f"gstlearn spectral (ns={n})" for n in NSPEC]
for model, (itype, ctype) in MODELS.items():
    for label in labels:
        if label.startswith("SPECSIM3D"):
            fields, dt = run_specsim3d(itype, (NX, NY, NZ), NREAL)
        else:
            fields, dt = run_gstlearn(ctype, (NX, NY, NZ), NREAL, int(label.split("ns=")[1].rstrip(")")))
        res = analyse(fields)
        methods[(model, label)] = dict(res=res, metrics=metrics(res, itype), time=dt, var=float(fields.var(axis=(1, 2, 3)).mean()))
        print(f"{model:12s} {label:28s} t={dt:7.2f}s  var={methods[(model, label)]['var']:.3f}", flush=True)

lines = ["model;method;time_s;mean_variance;" +
         ";".join(f"RMS_mean_{d[0].split()[0]};RMS_single_{d[0].split()[0]};error_lag1_{d[0].split()[0]}" for d in DIRS)]
for (model, label), d in methods.items():
    row = [model, label, f"{d['time']:.2f}", f"{d['var']:.3f}"]
    for rms_m, rms_s, first in d["metrics"]:
        row += [f"{rms_m:.3f}", f"{rms_s:.3f}", f"{first:+.3f}"]
    lines.append(";".join(row))
open(f"{RES}/gstlearn_comparison_metrics.csv", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

fig, axs = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
colors = {labels[0]: "tab:blue", labels[1]: "tab:purple", labels[2]: "tab:brown"}
for i, (model, (itype, _)) in enumerate(MODELS.items()):
    for j, (name, step, ax_, kmax) in enumerate(DIRS):
        ax = axs[i, j]
        for label in labels:
            l, gm, gall = methods[(model, label)]["res"][j]
            lo, hi = np.percentile(gall, [5, 95], axis=0)
            ax.fill_between(l, lo, hi, color=colors[label], alpha=0.15, lw=0)
            ax.plot(l, gm, color=colors[label], lw=1.5, label=label)
        ax.plot(l, gamma_model(itype, l, RANGES[ax_]), "k-", lw=2.4, label="model")
        ax.set_title(f"{model} - {name}", fontsize=9); ax.grid(alpha=0.3); ax.set_ylim(0, 1.3)
axs[0, 0].legend(fontsize=7)
for ax in axs[-1]:
    ax.set_xlabel("lag")
for ax in axs[:, 0]:
    ax.set_ylabel("semivariogram (mean of 20, band 5-95 %)")
plt.tight_layout()
plt.savefig(f"{RES}/gstlearn_comparison.png", dpi=120)
json.dump({"gstlearn": getattr(gl, "__version__", "?"), "spherical": "not supported by simuSpectral"}, open(f"{RES}/gstlearn_version.json", "w"))
