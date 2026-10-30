"""Tests of SPECSIM3D (and of GSTools where applicable) with different ranges and anisotropy ratios.

For each case the experimental variogram (mean over NREAL realizations) is compared with the
analytic one at ~300 random lag vectors (not only along the axes), so that the whole 3D
structure (ranges, ratios, rotations) is verified. Metrics:
  RMS      : root-mean-square deviation of the mean variogram from the model
  noise    : expected standard error of the mean variogram (from the spread between realizations)
  RMS/noise ~ 1 -> deviation compatible with Monte Carlo noise alone; >> 1 -> systematic bias
  bias     : signed mean error
Autore: Stefano Bonduà
"""
import subprocess, time, csv, sys, os
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
SEED = 2024

# name: (nugget, [structures (sill, a1, a2, a3, az, dip, plg)], gstools_ok)
CASES = {
    "A iso 15/15/15": (0.0, [(1.0, 15, 15, 15, 0, 0, 0)], True),
    "B 2:1 20/10/10 az30": (0.0, [(1.0, 20, 10, 10, 30, 0, 0)], True),
    "C 4:1 20/5/5 az60": (0.0, [(1.0, 20, 5, 5, 60, 0, 0)], True),
    "D layered 10/10/2": (0.0, [(1.0, 10, 10, 2, 0, 0, 0)], True),
    "E extreme 40/4/2 az135": (0.0, [(1.0, 40, 4, 2, 135, 0, 0)], True),
    "F rotated 20/8/4 az30 dip20 plg15": (0.0, [(1.0, 20, 8, 4, 30, 20, 15)], False),
}
NESTED = ("G nested+nugget", 0.1, [(1, 0.5, 20, 10, 5, 45, 0, 0), (2, 0.4, 8, 8, 3, 0, 0, 0)])
MODELS = {"spherical": 1, "exponential": 2, "gaussian": 3}


# --------------------- analytic model (GSLIB rotation, rewritten in Python) ----------------------
def rotmat(a1, a2, a3, az, dip, plg):
    d = np.pi / 180
    alpha = (90 - az) * d if 0 <= az < 270 else (450 - az) * d
    beta, theta = -dip * d, plg * d
    sa, ca, sb, cb, st, ct = np.sin(alpha), np.cos(alpha), np.sin(beta), np.cos(beta), np.sin(theta), np.cos(theta)
    f1, f2 = a1 / a2, a1 / a3
    return np.array([[cb * ca, cb * sa, -sb],
                     [f1 * (-ct * sa + st * sb * ca), f1 * (ct * ca + st * sb * sa), f1 * st * cb],
                     [f2 * (st * sa + ct * sb * ca), f2 * (-st * ca + ct * sb * sa), f2 * ct * cb]])


def corr(itype, hr):
    if itype == 1:
        return np.where(hr < 1, 1 - 1.5 * hr + 0.5 * hr**3, 0.0)
    if itype == 2:
        return np.exp(-3 * hr)
    return np.exp(-3 * hr**2)


def cov_th(d, nugget, structs):
    """d: (n,3) lag vectors; structs: [(itype,sill,a1,a2,a3,az,dip,plg)]"""
    d = np.atleast_2d(d).astype(float)
    c = np.zeros(len(d))
    for it, sill, a1, a2, a3, az, dip, plg in structs:
        r = d @ rotmat(a1, a2, a3, az, dip, plg).T
        c += sill * corr(it, np.linalg.norm(r, axis=1) / a1)
    c += nugget * (np.linalg.norm(d, axis=1) == 0)
    return c


# --------------------- lag vectors ----------------------------------------------------------
rng = np.random.default_rng(7)


def pick_lags(nugget, structs, n=300):
    out = []
    tot = nugget + sum(s[1] for s in structs)
    while len(out) < n:
        d = np.array([rng.integers(-30, 31), rng.integers(-30, 31), rng.integers(-12, 13)])
        if not d.any():
            continue
        g = tot - cov_th(d, nugget, structs)[0]
        if 0.05 * tot < g < 0.97 * tot:           # informative lags (neither saturated nor null)
            out.append(d)
    return np.array(out)


def vario_vec(fields, d):
    """Variogram for the lag vector d=(dx,dy,dz) over all pairs; fields: (nreal,nx,ny,nz)."""
    sl_a, sl_b = [slice(None)], [slice(None)]
    for di, n in zip(d, fields.shape[1:]):
        sl_a.append(slice(max(0, -di), n - max(0, di)))
        sl_b.append(slice(max(0, di), n - max(0, -di)))
    diff = fields[tuple(sl_a)] - fields[tuple(sl_b)]
    return 0.5 * np.mean(diff**2, axis=(1, 2, 3))      # per realization


# --------------------- generators -------------------------------------------------------------
def run_specsim3d(tag, nugget, structs_with_type, ipad=1):
    lines = [f"0 0 0", "1 1 1", f"{NX} {NY} {NZ}", f"{nugget} {len(structs_with_type)}"]
    for it, sill, a1, a2, a3, az, dip, plg in structs_with_type:
        lines.append(f"{it} {sill} {a1} {a2} {a3} {az} {dip} {plg}")
    lines += [f"{SEED} {NREAL} {ipad}", tag]
    open(f"{WORK}/{tag}.par", "w").write("\n".join(lines) + "\n")
    t0 = time.time()
    subprocess.run([EXE, f"{WORK}/{tag}.par", "--bin"], cwd=WORK, check=True, capture_output=True)
    dt = time.time() - t0
    a = np.fromfile(f"{WORK}/{tag}.bin", dtype=np.float32).reshape(NREAL, NZ, NY, NX).astype(float)
    return a.transpose(0, 3, 2, 1), dt


def run_gstools(tag, it, a1, a2, a3, az, d_check=None):
    cache = f"{WORK}/cases_gs_{tag}.npy"
    if it == 1:
        M, ls = gs.Spherical, np.array([a1, a2, a3], float)
    elif it == 2:
        M, ls = gs.Exponential, np.array([a1, a2, a3], float) / 3.0
    else:
        M, ls = gs.Gaussian, np.array([a1, a2, a3], float) * np.sqrt(np.pi / 12.0)
    mod = M(dim=3, var=1.0, len_scale=list(ls), angles=[np.radians(90 - az), 0.0, 0.0])
    if d_check is not None:   # is the orientation of the GSTools model the analytic one?
        hiso = np.linalg.norm(mod.isometrize(d_check.T.astype(float)), axis=0)
        g_gs = mod.variogram(hiso)
        g_th = 1 - cov_th(d_check, 0.0, [(it, 1.0, a1, a2, a3, az, 0, 0)])
        assert np.allclose(g_gs, g_th, atol=1e-6), f"{tag}: GSTools orientation != analytic, max {np.abs(g_gs-g_th).max():.3g}"
    if os.path.exists(cache):
        return np.load(cache), 0.0
    x, y, z = np.arange(NX, dtype=float), np.arange(NY, dtype=float), np.arange(NZ, dtype=float)
    t0 = time.time()
    out = np.array([gs.SRF(mod, seed=SEED + i, mode_no=1024).structured((x, y, z)) for i in range(NREAL)])
    dt = time.time() - t0
    np.save(cache, out)
    return out, dt


def evaluate(fields, nugget, structs, lags):
    tot = nugget + sum(s[1] for s in structs)
    g_th = tot - cov_th(lags, nugget, structs)
    G = np.array([vario_vec(fields, d) for d in lags])        # (nlag, nreal)
    g_mean, g_se = G.mean(1), G.std(1, ddof=1) / np.sqrt(G.shape[1])
    err = g_mean - g_th
    return dict(rms=np.sqrt(np.mean(err**2)), noise=np.sqrt(np.mean(g_se**2)), bias=err.mean(),
                var=fields.var(axis=(1, 2, 3)).mean(), g_th=g_th, g_mean=g_mean, tot=tot)


# --------------------- run -------------------------------------------------------------
rows, scat = [], {}


def record(case, model, method, ev, dt):
    rows.append(dict(case=case, model=model, method=method, time_s=round(dt, 1), RMS=ev["rms"], MC_noise=ev["noise"],
                     RMS_over_noise=ev["rms"] / ev["noise"], bias=ev["bias"], mean_variance=ev["var"], sill=ev["tot"]))
    scat[(case, model, method)] = (ev["g_th"], ev["g_mean"])
    print(f"{case:34s} {model:13s} {method:10s} RMS={ev['rms']:.4f} noise={ev['noise']:.4f} ratio={ev['rms']/ev['noise']:.2f} bias={ev['bias']:+.4f}", flush=True)


for case, (nug, st, gs_ok) in CASES.items():
    for model, it in MODELS.items():
        structs = [(it,) + s for s in st]
        lags = pick_lags(nug, structs)
        tag = "case_" + case.split()[0] + "_" + model
        f, dt = run_specsim3d(tag, nug, structs)
        record(case, model, "SPECSIM3D", evaluate(f, nug, structs, lags), dt)
        if gs_ok and it in (2, 3, 1):
            s0 = st[0]
            fg, dtg = run_gstools(tag, it, s0[1], s0[2], s0[3], s0[4], d_check=lags[:50])
            record(case, model, "GSTools", evaluate(fg, nug, structs, lags), dtg)

name, nug, structs = NESTED
lags = pick_lags(nug, structs)
f, dt = run_specsim3d("case_G_nested", nug, structs)
record(name, "sph+exp", "SPECSIM3D", evaluate(f, nug, structs, lags), dt)

with open(f"{RES}/anisotropy_cases.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter=";")
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})

# figure: mean vs theoretical variogram (scatter) for each method
fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
for ax, method in zip(axs, ["SPECSIM3D", "GSTools"]):
    for (case, model, m), (gth, gm) in scat.items():
        if m != method:
            continue
        ax.plot(gth, gm, ".", ms=3, alpha=0.5, color={"spherical": "tab:green", "exponential": "tab:blue", "gaussian": "tab:orange"}.get(model, "tab:purple"))
    ax.plot([0, 1], [0, 1], "k-", lw=1)
    ax.set_title(method + " - mean variogram vs model (all cases)")
    ax.set_xlabel("theoretical gamma")
    ax.grid(alpha=0.3)
axs[0].set_ylabel("experimental gamma (mean of the realizations)")
plt.tight_layout()
plt.savefig(f"{RES}/anisotropy_cases.png", dpi=120)
