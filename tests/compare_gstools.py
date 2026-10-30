"""Comparison of SPECSIM3D (Fortran) with GSTools (SRF, randomization method) in 3D.

Test case (Yao, 1998): 64x64x16 grid, unit spacing, unit variance, practical ranges 20/10/5,
major direction N45E. Three models: spherical, exponential, Gaussian.
Evaluation: exact experimental variograms (all pairs along grid shifts) in the three principal
directions, averaged over NREAL realizations and compared with the analytic variogram.
Autore: Stefano Bonduà
"""
import subprocess, time, json, os
import numpy as np
import gstools as gs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os as _os
WD = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.path.dirname(WD)
EXE = _os.path.join(ROOT, "specsim3d.exe" if _os.name == "nt" else "specsim3d")
WORK = _os.path.join(WD, "work")      # generated .par/.bin/.out files
RES = _os.path.join(WD, "results")    # csv/png/json results
_os.makedirs(WORK, exist_ok=True); _os.makedirs(RES, exist_ok=True)
NX, NY, NZ = 64, 64, 16
NREAL = 20
SEED = 31791191
RANGES = (20.0, 10.0, 5.0)
MODELS = {"spherical": 1, "exponential": 2, "gaussian": 3}


# ---------------- analytic model -------------------------------------------
def corr(itype, h, a):
    r = np.asarray(h, float) / a
    if itype == 1:
        return np.where(r < 1, 1 - 1.5 * r + 0.5 * r**3, 0.0)
    if itype == 2:
        return np.exp(-3 * r)
    return np.exp(-3 * r**2)


def gamma_model(itype, h, a):
    return 1.0 - corr(itype, h, a)


# ---------------- exact experimental variograms -------------------------------
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
    """fields: (nreal, nx, ny, nz) -> for each direction: lags, mean gamma, gamma of every realization."""
    res = []
    for name, step, ax, kmax in DIRS:
        gs_ = []
        for f in fields:
            l, g = vario_shift(f, step, kmax)
            gs_.append(g)
        res.append((l, np.mean(gs_, 0), np.array(gs_)))
    return res


def metrics(res, itype):
    """RMS error of the mean variogram and of the single realizations within the range."""
    rows = []
    for (name, step, ax, kmax), (l, gm, gall) in zip(DIRS, res):
        a = RANGES[ax]
        m = l <= a
        th = gamma_model(itype, l, a)
        rms_mean = np.sqrt(np.mean((gm[m] - th[m]) ** 2))
        rms_single = np.mean(np.sqrt(np.mean((gall[:, m] - th[m]) ** 2, axis=1)))
        first = (gm[0] - th[0])
        rows.append((rms_mean, rms_single, first))
    return rows


# ---------------- SPECSIM3D ------------------------------------------------------
def run_specsim3d(model, itype, ipad):
    tag = f"cmp_{model}_pad{ipad}"
    with open(f"{WORK}/{tag}.par", "w") as fh:
        fh.write(f"0.5 0.5 0.5\n1 1 1\n{NX} {NY} {NZ}\n0.0 1\n"
                 f"{itype} 1.0 {RANGES[0]} {RANGES[1]} {RANGES[2]} 45.0 0.0 0.0\n"
                 f"{SEED} {NREAL} {ipad}\n{tag}\n")
    t0 = time.time()
    subprocess.run([EXE, f"{WORK}/{tag}.par"], cwd=WORK, check=True, capture_output=True)
    dt = time.time() - t0
    raw = np.loadtxt(f"{WORK}/{tag}.out", skiprows=2 + NREAL)
    fields = np.array([raw[:, r].reshape(NZ, NY, NX).transpose(2, 1, 0) for r in range(NREAL)])
    return fields, dt


# ---------------- GSTools --------------------------------------------------------
def run_gstools(model, itype, mode_no=1024):
    if itype == 1:
        M, ls = gs.Spherical, np.array(RANGES)
    elif itype == 2:
        M, ls = gs.Exponential, np.array(RANGES) / 3.0
    else:
        M, ls = gs.Gaussian, np.array(RANGES) * np.sqrt(np.pi / 12.0)
    mod = M(dim=3, var=1.0, len_scale=list(ls), angles=[np.pi / 4, 0.0, 0.0])
    # check that the GSTools model coincides with the analytic model of SPECSIM3D
    for ax in range(3):
        h = np.linspace(0.5, 30, 30)
        assert np.allclose(mod.vario_axis(h, axis=ax), gamma_model(itype, h, RANGES[ax]), atol=1e-6), (model, ax)
    x, y, z = np.arange(NX, dtype=float), np.arange(NY, dtype=float), np.arange(NZ, dtype=float)
    t0 = time.time()
    fields = np.array([gs.SRF(mod, seed=SEED + i, mode_no=mode_no).structured((x, y, z)) for i in range(NREAL)])
    return fields, time.time() - t0


# ---------------- run -----------------------------------------------------
methods = {}
for model, itype in MODELS.items():
    for label, fn in [("SPECSIM3D legacy (ipad=0)", lambda m=model, t=itype: run_specsim3d(m, t, 0)),
                      ("SPECSIM3D padded (ipad=1)", lambda m=model, t=itype: run_specsim3d(m, t, 1)),
                      ("GSTools SRF (1024 modes)", lambda m=model, t=itype: run_gstools(m, t, 1024))]:
        fields, dt = fn()
        res = analyse(fields)
        methods[(model, label)] = dict(res=res, metrics=metrics(res, itype), time=dt,
                                       mean=float(fields.mean()), var=float(fields.var(axis=(1, 2, 3)).mean()))
        print(f"{model:13s} {label:28s} t={dt:6.1f}s  var={methods[(model,label)]['var']:.3f}", flush=True)

# ---------------- table --------------------------------------------------------
lines = ["model;method;time_s;mean_variance;" +
         ";".join(f"RMS_mean_{d[0].split()[0]};RMS_single_{d[0].split()[0]};error_lag1_{d[0].split()[0]}" for d in DIRS)]
for (model, label), d in methods.items():
    row = [model, label, f"{d['time']:.1f}", f"{d['var']:.3f}"]
    for rms_m, rms_s, first in d["metrics"]:
        row += [f"{rms_m:.3f}", f"{rms_s:.3f}", f"{first:+.3f}"]
    lines.append(";".join(row))
open(f"{RES}/gstools_comparison_metrics.csv", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

# ---------------- figure ---------------------------------------------------------
fig, axs = plt.subplots(3, 3, figsize=(14, 10), sharey=True)
colors = {"SPECSIM3D legacy (ipad=0)": "tab:green", "SPECSIM3D padded (ipad=1)": "tab:blue", "GSTools SRF (1024 modes)": "tab:red"}
for i, (model, itype) in enumerate(MODELS.items()):
    for j, (name, step, ax_, kmax) in enumerate(DIRS):
        ax = axs[i, j]
        for label in colors:
            l, gm, _ = methods[(model, label)]["res"][j]
            ax.plot(l, gm, color=colors[label], lw=1.5, label=label)
        l = methods[(model, "GSTools SRF (1024 modes)")]["res"][j][0]
        ax.plot(l, gamma_model(itype, l, RANGES[ax_]), "k-", lw=2.4, label="model")
        ax.set_title(f"{model} - {name}", fontsize=9)
        ax.grid(alpha=0.3)
        ax.set_ylim(0, 1.3)
axs[0, 0].legend(fontsize=7)
for ax in axs[2]:
    ax.set_xlabel("lag")
for ax in axs[:, 0]:
    ax.set_ylabel("semivariogram (mean of 20 realizations)")
plt.tight_layout()
plt.savefig(f"{RES}/gstools_comparison.png", dpi=120)
json.dump({f"{m}|{l}": dict(time=d["time"], var=d["var"]) for (m, l), d in methods.items()},
          open(f"{RES}/gstools_comparison_times.json", "w"), indent=1)
