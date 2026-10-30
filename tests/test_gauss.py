"""Gaussianity tests of SPECSIM3D vs GSTools (SRF, randomization method).

NREAL independent realizations are used on a small grid (16x16x8): the values of the same node in
different realizations are independent samples, so the tests (Shapiro-Wilk, skewness, kurtosis)
are valid and are not affected by the spatial correlation.
 1) marginal: for each node Shapiro-Wilk, skewness, excess kurtosis, mean and variance;
    under H0 the p-values are uniform (expected fraction p < 0.05 about 5 %).
 2) multivariate: Mardia kurtosis on vectors of p = 4 neighbouring nodes (expected p(p+2) = 24).
 3) bivariate tail: P(z1 > q, z2 > q) observed vs that of a bivariate normal with the same
    empirical correlation, for pairs at a fixed distance.
Autore: Stefano Bonduà
"""
import os, subprocess, time
import numpy as np
import gstools as gs
from scipy import stats
from scipy.stats import multivariate_normal
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
NX, NY, NZ = 16, 16, 8
NREAL = 2000
SEED = 4242
RANGES = (6.0, 3.0, 2.0)
MODELS = {"spherical": 1, "exponential": 2, "gaussian": 3}
rng = np.random.default_rng(1)


def fields_specsim3d(model, itype):
    tag = f"gauss_{model}"
    open(f"{WORK}/{tag}.par", "w").write(
        f"0 0 0\n1 1 1\n{NX} {NY} {NZ}\n0.0 1\n{itype} 1.0 {RANGES[0]} {RANGES[1]} {RANGES[2]} 45 0 0\n{SEED} {NREAL} 1\n{tag}\n")
    t0 = time.time()
    subprocess.run([EXE, f"{WORK}/{tag}.par", "--bin"], cwd=WORK, check=True, capture_output=True)
    dt = time.time() - t0
    a = np.fromfile(f"{WORK}/{tag}.bin", dtype=np.float32).reshape(NREAL, NZ, NY, NX).astype(float)
    return a.transpose(0, 3, 2, 1), dt           # (nreal, x, y, z)


def fields_gstools(model, itype):
    cache = f"{WORK}/gauss_gs_{model}.npy"
    if os.path.exists(cache):
        return np.load(cache), 0.0
    if itype == 1:
        M, ls = gs.Spherical, np.array(RANGES)
    elif itype == 2:
        M, ls = gs.Exponential, np.array(RANGES) / 3.0
    else:
        M, ls = gs.Gaussian, np.array(RANGES) * np.sqrt(np.pi / 12.0)
    mod = M(dim=3, var=1.0, len_scale=list(ls), angles=[np.pi / 4, 0.0, 0.0])
    x, y, z = np.arange(NX, dtype=float), np.arange(NY, dtype=float), np.arange(NZ, dtype=float)
    t0 = time.time()
    out = np.array([gs.SRF(mod, seed=SEED + i, mode_no=1024).structured((x, y, z)) for i in range(NREAL)])
    dt = time.time() - t0
    np.save(cache, out)
    return out, dt


def marginal_tests(f):
    v = f.reshape(NREAL, -1)                          # (nreal, nodes)
    mean, var = v.mean(0), v.var(0, ddof=1)
    skew, kurt = stats.skew(v, axis=0), stats.kurtosis(v, axis=0)
    sw = np.array([stats.shapiro(v[:, n]).pvalue for n in range(v.shape[1])])
    ks = np.array([stats.kstest(v[:, n], "norm").pvalue for n in range(v.shape[1])])
    return dict(mean=mean, var=var, skew=skew, kurt=kurt, sw=sw, ks=ks)


def mardia(f, ntriples=60, p=4):
    """Mardia kurtosis on p neighbouring nodes; returns z-scores (expected 0, sd 1)."""
    n = NREAL
    zs = []
    for _ in range(ntriples):
        x0, y0, z0 = rng.integers(0, NX - 3), rng.integers(0, NY - 3), rng.integers(0, NZ - 2)
        pts = [(x0, y0, z0), (x0 + 2, y0, z0), (x0, y0 + 2, z0), (x0 + 1, y0 + 1, z0 + 1)][:p]
        X = np.stack([f[:, a, b, c] for a, b, c in pts], axis=1)
        Xc = X - X.mean(0)
        S = np.cov(Xc.T, bias=True)
        d2 = np.einsum("ij,jk,ik->i", Xc, np.linalg.inv(S), Xc)
        b = np.mean(d2**2)
        zs.append((b - p * (p + 2)) / np.sqrt(8 * p * (p + 2) / n))
    return np.array(zs)


def tail_bivariate(f, lag, q):
    pairs = []
    sl = lambda a, d, n: slice(0, n - d) if d >= 0 else slice(-d, n)
    sl2 = lambda a, d, n: slice(d, n) if d >= 0 else slice(0, n + d)
    A = f[:, sl(0, lag[0], NX), :, :][:, :, sl(0, lag[1], NY), :][:, :, :, sl(0, lag[2], NZ)]
    B = f[:, sl2(0, lag[0], NX), :, :][:, :, sl2(0, lag[1], NY), :][:, :, :, sl2(0, lag[2], NZ)]
    a, b = A.ravel(), B.ravel()
    rho = np.corrcoef(a, b)[0, 1]
    obs = np.mean((a > q) & (b > q))
    th = multivariate_normal(mean=[0, 0], cov=[[1, rho], [rho, 1]]).cdf([-q, -q])  # P(z1>q, z2>q) by symmetry
    return rho, obs, th


rows, store = [], {}
for model, itype in MODELS.items():
    for label, fn in [("SPECSIM3D", fields_specsim3d), ("GSTools", fields_gstools)]:
        f, dt = fn(model, itype)
        m = marginal_tests(f)
        z = mardia(f)
        t_x = tail_bivariate(f, (2, 0, 0), 1.0)
        t_x2 = tail_bivariate(f, (2, 0, 0), 1.645)
        t_z = tail_bivariate(f, (0, 0, 1), 1.0)
        store[(model, label)] = (f, m)
        rows.append(dict(model=model, method=label, time_s=round(dt, 1),
                         mean_nodes=m["mean"].mean(), var_nodes=m["var"].mean(), sd_var_nodes=m["var"].std(),
                         skew_mean=m["skew"].mean(), kurt_mean=m["kurt"].mean(),
                         SW_frac_p005=np.mean(m["sw"] < 0.05), KS_frac_p005=np.mean(m["ks"] < 0.05),
                         Mardia_z_mean=z.mean(), Mardia_frac_abs_z_gt196=np.mean(np.abs(z) > 1.96),
                         tail_x_rho=t_x[0], tail_x_obs=t_x[1], tail_x_theory=t_x[2],
                         tail1645_obs=t_x2[1], tail1645_theory=t_x2[2],
                         tail_z_rho=t_z[0], tail_z_obs=t_z[1], tail_z_theory=t_z[2]))
        print(model, label, "done", flush=True)

import csv
with open(f"{RES}/gaussianity.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter=";")
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
for r in rows:
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})

# figure: QQ plot of the standardized values (exponential) and histogram of the Shapiro p-values
fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
cols = {"SPECSIM3D": "tab:blue", "GSTools": "tab:red"}
for label in cols:
    f, m = store[("exponential", label)]
    v = (f.reshape(NREAL, -1)[:, ::64]).ravel()
    v = (v - v.mean()) / v.std()
    osm, osr = stats.probplot(v, dist="norm", fit=False)
    sub = slice(None, None, 25)
    axs[0].plot(osm[sub], osr[sub], ".", ms=2, color=cols[label], label=label)
    axs[1].hist(m["sw"], bins=20, range=(0, 1), histtype="step", color=cols[label], lw=1.6, label=label)
    axs[2].hist(m["kurt"], bins=30, histtype="step", color=cols[label], lw=1.6, label=label)
axs[0].plot([-4, 4], [-4, 4], "k-", lw=1)
axs[0].set_title("QQ plot (exponential)"); axs[0].set_xlabel("normal quantiles"); axs[0].set_ylabel("observed quantiles")
axs[1].set_title("Shapiro-Wilk p-value per node (uniform if Gaussian)"); axs[1].set_xlabel("p")
axs[2].set_title("excess kurtosis per node (0 if Gaussian)")
axs[0].legend(); axs[1].legend()
plt.tight_layout()
plt.savefig(f"{RES}/gaussianity.png", dpi=120)
