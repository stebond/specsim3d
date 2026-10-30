#!/usr/bin/env python3
"""
make_figures.py - regenerates every figure (and the numeric tables) of the SPECSIM3D manual.

Usage (from any folder):
    python make_figures.py                 # all figures, GSTools comparisons included if gstools is installed
    python make_figures.py --no-gstools    # only SPECSIM3D figures (fast, ~1 min)
    python make_figures.py --full          # also re-measure scaling (32^3..256^3) and variogram speed up to 256^3
    python make_figures.py --only fig07 fig05   # only some figures (fig01..fig10 = figures of the paper, figS1 = manual only)

Requirements: numpy, scipy, matplotlib (psutil optional, gstools optional) and the compiled
specsim3d executable in the parent folder (build_specsim3d.bat / gfortran, see README.txt).

Everything is written to docs/figures (PNG, 300 dpi, and vector PDF) and docs/figures/tables (CSV).
File names follow the figure numbers of the paper (Computers & Geosciences manuscript):
    fig01 geometry, fig02 algorithm, fig03 slices, fig04 models vs GSTools, fig05 anisotropy cases,
    fig06 Gaussianity, fig07 legacy vs padded grid, fig08 2D validation, fig09 scaling, fig10 variogram speed,
    fig11-fig14 conditional simulation (maps, validation against kriging, normal score, speed);
    figS1 (variograms of the Yao case) appears only in the manual.
Machine, compiler and package versions are stored in docs/figures/system_info.json and printed in
the text of the paper (not in the figures).

Author: Stefano Bonduà (DICAM, University of Bologna)
"""
import argparse, csv, json, os, platform, re, shutil, subprocess, sys, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

DOCS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(DOCS)
FIG = os.path.join(DOCS, "figures")
TAB = os.path.join(FIG, "tables")
WORK = os.path.join(FIG, "_work")
EXE = os.path.join(ROOT, "specsim3d.exe" if os.name == "nt" else "specsim3d")
ORIG2D = os.path.join(DOCS, "data", "specsim2_original_output.out")
DPI = 300
C_SPEC, C_GS, C_LEG, C_MOD = "#1f77b4", "#d62728", "#2ca02c", "black"
C_GSM = "#ff7f0e"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "axes.titlesize": 9.5,
                     "figure.dpi": 100, "savefig.dpi": DPI, "savefig.bbox": "tight", "pdf.fonttype": 42})

try:
    import gstools as gs
    HAVE_GS = True
except Exception:
    HAVE_GS = False
try:
    from scipy import stats
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

RNG = np.random.default_rng(12345)


# ============================================================================ utilities
def system_info():
    info = {"date": time.strftime("%Y-%m-%d %H:%M"), "platform": platform.platform(),
            "python": sys.version.split()[0], "numpy": np.__version__, "cores": os.cpu_count(),
            "omp_threads": os.environ.get("OMP_NUM_THREADS", str(os.cpu_count()))}
    cpu = platform.processor()
    if os.path.exists("/proc/cpuinfo"):
        m = re.search(r"model name\s*:\s*(.+)", open("/proc/cpuinfo").read())
        cpu = m.group(1).strip() if m else cpu
    if os.name == "nt":
        try:
            o = subprocess.run(["powershell", "-NoProfile", "-Command",
                                "(Get-CimInstance Win32_Processor | Select-Object -First 1).Name"],
                               capture_output=True, text=True, timeout=30).stdout.strip()
            cpu = o or cpu
        except Exception:
            pass
    info["cpu"] = cpu or "unknown"
    ram = None
    try:
        import psutil
        ram = psutil.virtual_memory().total / 2**30
    except Exception:
        if os.path.exists("/proc/meminfo"):
            ram = int(re.search(r"MemTotal:\s+(\d+)", open("/proc/meminfo").read()).group(1)) / 2**20
    if not ram and os.name == "nt":
        try:
            o = subprocess.run(["powershell", "-NoProfile", "-Command",
                                "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"],
                               capture_output=True, text=True, timeout=30).stdout.strip()
            ram = int(o) / 2**30 if o.isdigit() else None
        except Exception:
            ram = None
    info["ram_gb"] = round(ram, 1) if ram else None
    comp = "unknown (edit system_info.json or set SPECSIM3D_COMPILER)"
    for cmd in (["gfortran", "--version"], ["ifx", "/QV"], ["ifx", "--version"], ["ifort", "--version"]):
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
            txt = (out.stdout + out.stderr).strip().splitlines()
            if out.returncode == 0 and txt:
                comp = txt[0]
                break
        except Exception:
            pass
    for fn in ("build_info.txt", "compiler.txt"):    # written by build_specsim3d.bat / by run_full_validation.py
        cf = os.path.join(ROOT, fn)
        if os.path.exists(cf):
            txt = " ".join(l.strip() for l in open(cf, encoding="utf-8", errors="replace") if l.strip())
            if txt and not txt.startswith("not found"):
                comp = txt                              # the record of the build beats whatever compiler is in PATH now
                break
    info["compiler"] = os.environ.get("SPECSIM3D_COMPILER", comp)
    if HAVE_GS:
        info["gstools"] = gs.__version__
    try:
        import scipy
        info["scipy"] = scipy.__version__
    except Exception:
        pass
    info["matplotlib"] = matplotlib.__version__
    return info


SYS = None


def footer(fig, extra=""):
    """No text is drawn under the figures any more: the hardware and compiler used for the benchmark
    are given in the text of the paper and in docs/figures/system_info.json, and the notes (extra) in the captions."""
    return


BAND_NOTE = "shaded bands: 5th-95th percentile of the individual variograms"


def band(ax, h, per, color, alpha=0.18, label=None):
    """Envelope (5th-95th percentile) of the variograms of the individual realizations.
    per: array (n_lags, nreal) or (nreal, n_lags) -- pass it as (nreal, n_lags)."""
    lo, hi = np.percentile(per, [5, 95], axis=0)
    ax.fill_between(h, lo, hi, color=color, alpha=alpha, lw=0, label=label)


def save(fig, name):
    path = os.path.join(FIG, name + ".png")
    fig.savefig(path)
    fig.savefig(os.path.join(FIG, name + ".pdf"))   # vector version (journal artwork requirement)
    plt.close(fig)
    print("  wrote", os.path.relpath(path, DOCS), flush=True)


def write_table(name, header, rows):
    os.makedirs(TAB, exist_ok=True)
    with open(os.path.join(TAB, name + ".csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def make_par(base, grid=(64, 64, 16), nugget=0.0, structs=((1, 1.0, 20, 10, 5, 45, 0, 0),), seed=31791191,
             nreal=20, ipad=1, spacing=(1, 1, 1)):
    lines = ["Parameters for SPECSIM3D", "************************", "", "START OF PARAMETERS:", "0 0 0", " ".join(map(str, spacing)), " ".join(map(str, grid)), f"{nugget} {len(structs)}"]
    lines += [" ".join(map(str, s)) for s in structs]
    lines += [f"{seed} {nreal} {ipad}", base]
    path = os.path.join(WORK, base + ".par")
    open(path, "w").write("\n".join(lines) + "\n")
    return path


def make_vario_par(base, dirs, iper=0):
    path = os.path.join(WORK, base + "_vp.par")
    open(path, "w").write(f"{len(dirs)} {iper}\n" + "\n".join(" ".join(map(str, d)) for d in dirs) + "\n")
    return path


def run_specsim(par, flags=("--bin",), vario=None):
    cmd = [EXE, par] + list(flags) + (["--vario", vario] if vario else [])
    t0 = time.time()
    out = subprocess.run(cmd, cwd=WORK, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stdout + out.stderr)
    return out.stdout, time.time() - t0


def load_bin(base, grid, nreal):
    nx, ny, nz = grid
    a = np.fromfile(os.path.join(WORK, base + ".bin"), dtype=np.float32).reshape(nreal, nz, ny, nx)
    return a.astype(float).transpose(0, 3, 2, 1)                       # (nreal, x, y, z)


def read_vario(base):
    rows = list(csv.DictReader(open(os.path.join(WORK, base + "_vario.csv"))))
    out = {}
    for r in rows:
        d = int(r["dir"])
        out.setdefault(d, {k: [] for k in ("lag", "h", "n", "model", "mean", "sd", "per")})
        o = out[d]
        o["lag"].append(int(r["lag"])); o["h"].append(float(r["h"])); o["n"].append(int(r["npairs"]))
        o["model"].append(float(r["gamma_model"])); o["mean"].append(float(r["gamma_mean"])); o["sd"].append(float(r["gamma_sd"]))
        o["per"].append([float(v) for k, v in r.items() if k.startswith("gamma_r")])
    for d in out:
        for k in out[d]:
            out[d][k] = np.array(out[d][k])
    return out


def vario_shift(f, step, kmax):
    """Reference (numpy) variogram of an (nreal,x,y,z) array along k*step: (nreal,k) and lags."""
    out, lags = [], []
    for k in range(1, kmax + 1):
        d = [k * s for s in step]
        if any(abs(di) >= n for di, n in zip(d, f.shape[1:])):
            break
        sa = [slice(None)] + [slice(max(0, -di), n - max(0, di)) for di, n in zip(d, f.shape[1:])]
        sb = [slice(None)] + [slice(max(0, di), n - max(0, -di)) for di, n in zip(d, f.shape[1:])]
        out.append(0.5 * np.mean((f[tuple(sa)] - f[tuple(sb)]) ** 2, axis=(1, 2, 3)))
        lags.append(k * np.linalg.norm(step))
    return np.array(out).T, np.array(lags)


def corr(it, r):
    r = np.asarray(r, float)
    if it == 1:
        return np.where(r < 1, 1 - 1.5 * r + 0.5 * r**3, 0.0)
    return np.exp(-3 * r) if it == 2 else np.exp(-3 * r**2)


MODELS = {"spherical": 1, "exponential": 2, "gaussian": 3}
RANGES = (20.0, 10.0, 5.0)
DIRS = [("N45E (a = 20)", (1, 1, 0), 0, 40), ("N45W (a = 10)", (1, -1, 0), 1, 40), ("vertical (a = 5)", (0, 0, 1), 2, 10)]
YAO_DIRS = [(1, 1, 0, 40), (1, -1, 0, 40), (0, 0, 1, 10)]


def gs_model(it, ranges=RANGES, az=45.0):
    a = np.array(ranges, float)
    if it == 1:
        M, ls = gs.Spherical, a
    elif it == 2:
        M, ls = gs.Exponential, a / 3.0
    else:
        M, ls = gs.Gaussian, a * np.sqrt(np.pi / 12.0)
    return M(dim=3, var=1.0, len_scale=list(ls), angles=[np.radians(90 - az), 0.0, 0.0])


def gs_fields(it, grid, nreal, seed, mode_no=1024):
    mod = gs_model(it)
    x, y, z = (np.arange(n, dtype=float) for n in grid)
    return np.array([gs.SRF(mod, seed=seed + i, mode_no=mode_no).structured((x, y, z)) for i in range(nreal)])


# ============================================================================ figures
def fig02_algorithm():
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 60); ax.grid(False)

    def box(x, y, w, h, txt, fc="#e8f1fa", ec="#1f77b4"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=1.2))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=8)

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10, color="0.25", lw=1.1))

    box(2, 42, 24, 13, "Parameter file\n(grid, nugget, nested\nstructures, ranges, angles)")
    box(34, 42, 28, 13, "Covariance C(h)\nsampled on a periodic grid\n(minimum-image lags)")
    box(70, 42, 28, 13, "3D FFT\nspectral density S(f)\n(real, S >= 0 check)")
    box(70, 20, 28, 13, "Amplitude sqrt(S)\n+ random phases with\nHermitian symmetry")
    box(34, 20, 28, 13, "Inverse 3D FFT\nreal field z(x)\n(cropped to n_x n_y n_z)")
    box(2, 20, 24, 13, "Outputs\n.out (GSLIB)  .bin  .vtk\n(ParaView, binary)")
    box(34, 2, 28, 12, "Experimental variograms\n(Fortran, OpenMP)\n--vario  ->  _vario.csv", fc="#eaf6ea", ec="#2ca02c")
    arrow(26, 48.5, 34, 48.5); arrow(62, 48.5, 70, 48.5); arrow(84, 42, 84, 33)
    arrow(70, 26.5, 62, 26.5); arrow(34, 26.5, 26, 26.5); arrow(48, 20, 48, 14)
    ax.text(50, 58.5, "SPECSIM3D workflow (repeated for every realization from the amplitude spectrum step)", ha="center", fontsize=9, weight="bold")
    save(fig, "fig02_algorithm")


def fig01_geometry():
    def rotmat(a1, a2, a3, az, dip, plg):                         # same formulas as specsim3d.f90 (GSLIB setrot)
        d = np.pi / 180
        alpha = (90 - az) * d if 0 <= az < 270 else (450 - az) * d
        beta, theta = -dip * d, plg * d
        sa, ca, sb, cb, st, ct = np.sin(alpha), np.cos(alpha), np.sin(beta), np.cos(beta), np.sin(theta), np.cos(theta)
        return np.array([[cb * ca, cb * sa, -sb], [-ct * sa + st * sb * ca, ct * ca + st * sb * sa, st * cb],
                         [st * sa + ct * sb * ca, -st * ca + ct * sb * sa, ct * cb]])   # unit principal directions (rows)
    fig = plt.figure(figsize=(8.4, 3.9))
    for k, (az, dip, plg, title) in enumerate([(45, 0, 0, "azimuth 45, dip 0, plunge 0"),
                                               (30, 20, 15, "azimuth 30, dip 20, plunge 15")]):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        U = rotmat(20, 10, 5, az, dip, plg)
        a = np.array([20, 10, 5.0])
        u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:15j]
        e = np.stack([a[0] * np.cos(u) * np.sin(v), a[1] * np.sin(u) * np.sin(v), a[2] * np.cos(v)], -1)
        p = e @ U
        ax.plot_wireframe(p[..., 0], p[..., 1], p[..., 2], color="0.75", lw=0.4)
        for i, (col, lab) in enumerate([("#d62728", "$a_{hmax}$"), ("#1f77b4", "$a_{hmin}$"), ("#2ca02c", "$a_{vert}$")]):
            vec = U[i] * a[i] * 1.25
            ax.plot([0, vec[0]], [0, vec[1]], [0, vec[2]], color=col, lw=2)
            ax.text(*(vec * 1.08), lab, color=col, fontsize=9)
        ax.plot([0, 0], [0, 28], [0, 0], "k--", lw=0.8); ax.text(0, 30, 0, "N", fontsize=8)
        ax.plot([0, 28], [0, 0], [0, 0], "k--", lw=0.8); ax.text(30, 0, 0, "E", fontsize=8)
        ax.set_xlim(-22, 22); ax.set_ylim(-22, 22); ax.set_zlim(-22, 22)
        ax.set_box_aspect((1, 1, 1)); ax.set_title(title, fontsize=9)
        ax.set_xlabel("x (East)"); ax.set_ylabel("y (North)"); ax.set_zlabel("z (up)")
        ax.view_init(24, -58)
    fig.suptitle("Anisotropy ellipsoid (ranges 20 / 10 / 5) and angle convention", fontsize=10)
    save(fig, "fig01_geometry")


def fig03_slices():
    base, grid = "f03", (64, 64, 16)
    run_specsim(make_par(base, grid, nreal=1, seed=31791191), ("--bin",))
    f = load_bin(base, grid, 1)[0]
    fig = plt.figure(figsize=(8.2, 6.4))
    gsp = fig.add_gridspec(2, 2, height_ratios=[3, 1.15], hspace=0.35, wspace=0.18)
    vmin, vmax = -3, 3
    ax = fig.add_subplot(gsp[0, :])
    im = ax.imshow(f[:, :, 8].T, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax, extent=(0, 64, 0, 64))
    ax.set_title("Horizontal slice z = 8 (major axis N45E, ranges 20 / 10 / 5, spherical)"); ax.set_xlabel("x"); ax.set_ylabel("y"); ax.grid(False)
    ax2 = fig.add_subplot(gsp[1, 0]); ax2.imshow(f[:, 32, :].T, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax, extent=(0, 64, 0, 16), aspect="equal")
    ax2.set_title("Vertical slice y = 32"); ax2.set_xlabel("x"); ax2.set_ylabel("z"); ax2.grid(False)
    ax3 = fig.add_subplot(gsp[1, 1]); ax3.imshow(f[32, :, :].T, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax, extent=(0, 64, 0, 16), aspect="equal")
    ax3.set_title("Vertical slice x = 32"); ax3.set_xlabel("y"); ax3.set_ylabel("z"); ax3.grid(False)
    fig.colorbar(im, ax=[ax, ax2, ax3], shrink=0.7, label="simulated value (unit variance)")
    save(fig, "fig03_slices")


def figS1_variograms_yao():
    base = "f04"
    par = make_par(base, nreal=20)
    vp = make_vario_par(base, YAO_DIRS, iper=1)
    run_specsim(par, ("--bin",), vp)
    v = read_vario(base)
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3))
    rows = []
    for k, ((name, step, ax_i, kmax), ax) in enumerate(zip(DIRS, axs), start=1):
        d = v[k]
        for r in range(d["per"].shape[1]):
            ax.plot(d["h"], d["per"][:, r], color="0.78", lw=0.7, label="realizations" if r == 0 else None)
        ax.fill_between(d["h"], d["mean"] - d["sd"], d["mean"] + d["sd"], color=C_SPEC, alpha=0.2, label="mean $\\pm$ 1 sd")
        ax.plot(d["h"], d["mean"], color=C_SPEC, lw=1.6, label="SPECSIM3D mean (20)")
        ax.plot(d["h"], d["model"], color=C_MOD, lw=2.0, label="model")
        ax.set_title(name); ax.set_xlabel("lag distance h"); ax.set_ylim(0, 1.45)
        m = d["model"] < 0.95
        rows.append([name, round(float(np.sqrt(np.mean((d["mean"][m] - d["model"][m]) ** 2))), 4), round(float(np.max(np.abs(d["mean"][m] - d["model"][m]))), 4)])
    axs[0].set_ylabel("semivariogram"); axs[0].legend(fontsize=7, loc="lower right")
    fig.suptitle("Spherical model, 64x64x16 grid, 20 realizations (variograms computed by SPECSIM3D --vario)", fontsize=9.5, y=1.03)
    write_table("figS1_yao_rms", ["direction", "rms_mean_vs_model", "max_abs_error"], rows)
    save(fig, "figS1_variograms_yao")


def gs_fields_modified(it, grid, nreal, seed):
    """GSTools fields with the author-modified spherical 3D sampling (separate process, see gs_modified_fields.py)."""
    out = os.path.join(WORK, f"gsmod_{it}.npy")
    subprocess.run([sys.executable, os.path.join(DOCS, "gs_modified_fields.py"), str(it), *map(str, grid), str(nreal), str(seed), out], check=True)
    return np.load(out)


def fig04_models(use_gs):
    grid, nreal = (64, 64, 16), 20
    fig, axs = plt.subplots(3, 3, figsize=(10.5, 8.8), sharey=True, gridspec_kw={"hspace": 0.38})
    rows = []
    for i, (mname, it) in enumerate(MODELS.items()):
        base = f"f05_{mname}"
        run_specsim(make_par(base, structs=((it, 1.0) + RANGES + (45, 0, 0),), nreal=nreal), ("--bin",), make_vario_par(base, YAO_DIRS, iper=1))
        v = read_vario(base)
        fg = gs_fields(it, grid, nreal, 31791191) if use_gs else None
        fm = gs_fields_modified(it, grid, nreal, 31791191) if use_gs and it == 1 else None   # only the spherical model is modified
        for j, (dname, step, ax_i, kmax) in enumerate(DIRS):
            ax = axs[i, j]
            d = v[j + 1]
            ax.plot(d["h"], d["model"], color=C_MOD, lw=2.3, label="model")
            band(ax, d["h"], d["per"].T, C_SPEC)
            ax.plot(d["h"], d["mean"], color=C_SPEC, lw=1.4, label="SPECSIM3D (mean of 20)")
            m = d["h"] <= RANGES[ax_i]
            rows.append([mname, "SPECSIM3D", dname, round(float(np.sqrt(np.mean((d["mean"][m] - d["model"][m]) ** 2))), 4)])
            if use_gs:
                g, l = vario_shift(fg, step, kmax)
                th = 1 - corr(it, l / RANGES[ax_i])
                band(ax, l, g, C_GS)
                ax.plot(l, g.mean(0), color=C_GS, lw=1.4, label="GSTools 1.7.0 (mean of 20)")
                mm = l <= RANGES[ax_i]
                rows.append([mname, "GSTools", dname, round(float(np.sqrt(np.mean((g.mean(0)[mm] - th[mm]) ** 2))), 4)])
                if fm is not None:
                    gm, lm = vario_shift(fm, step, kmax)
                    band(ax, lm, gm, C_GSM)
                    ax.plot(lm, gm.mean(0), color=C_GSM, lw=1.4, ls="--", label="GSTools modified (mean of 20)")
                    rows.append([mname, "GSTools (modified)", dname, round(float(np.sqrt(np.mean((gm.mean(0)[mm] - th[mm]) ** 2))), 4)])
            ax.set_title(f"{mname} - {dname}", fontsize=8.5); ax.set_ylim(0, 1.3)
            ax.set_xlabel("lag distance h")
    axs[0, 0].legend(fontsize=6.5, loc="lower right")
    for a in axs[:, 0]:
        a.set_ylabel("semivariogram")
    write_table("fig04_models_rms", ["model", "method", "direction", "rms_within_range"], rows)
    footer(fig, BAND_NOTE + " (20 realizations)")
    save(fig, "fig04_models_vs_gstools")


def fig07_legacy_vs_padded():
    grid, nreal = (64, 64, 16), 40
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.5))
    rows = []
    for ipad, col, lab in [(0, C_LEG, "legacy grid (ipad = 0)"), (1, C_SPEC, "padded grid (ipad = 1)")]:
        base = f"f06_g{ipad}"
        run_specsim(make_par(base, structs=((3, 1.0) + RANGES + (45, 0, 0),), nreal=nreal, ipad=ipad), ("--bin",), make_vario_par(base, [(1, 1, 0, 57)]))
        d = read_vario(base)[1]
        axs[0].plot(d["h"], d["mean"], color=col, lw=1.4, label=lab)
        if ipad == 0:
            axs[0].plot(d["h"], d["model"], color=C_MOD, lw=2, label="model")
            hmax = d["h"].max()
        m = d["h"] > 25
        rows.append([lab, "gaussian N45E, h > 25: mean gamma", round(float(d["mean"][m].mean()), 4)])
    axs[0].set_title("Gaussian model, N45E direction"); axs[0].set_xlabel("lag distance h"); axs[0].set_ylabel("semivariogram"); axs[0].legend(fontsize=7)
    axs[0].set_ylim(0.8, 1.15); axs[0].set_xlim(0, 60)
    nr2 = 300
    for ipad, col, lab in [(0, C_LEG, "legacy (ipad = 0)"), (1, C_SPEC, "padded (ipad = 1)")]:
        base = f"f06_s{ipad}"
        run_specsim(make_par(base, structs=((1, 1.0) + RANGES + (45, 0, 0),), nreal=nr2, ipad=ipad, seed=777), ("--bin",))
        f = load_bin(base, grid, nr2)
        var = f.reshape(nr2, -1).var(axis=1)
        mean = f.reshape(nr2, -1).mean(axis=1)
        axs[1].hist(var, bins=np.linspace(0.6, 1.5, 46), color=col, alpha=0.65, label=f"{lab}: sd = {var.std():.3f}")
        axs[2].hist(mean, bins=np.linspace(-0.6, 0.6, 49), color=col, alpha=0.65, label=f"{lab}: sd = {mean.std():.3f}")
        rows.append([lab, "spherical, 300 realizations: sd of sample variance", round(float(var.std()), 5)])
        rows.append([lab, "spherical, 300 realizations: sd of sample mean", round(float(mean.std()), 5)])
    axs[1].set_title("Sample variance of each realization"); axs[1].set_xlabel("variance"); axs[1].set_ylabel("count"); axs[1].legend(fontsize=7)
    axs[2].set_title("Sample mean of each realization"); axs[2].set_xlabel("mean"); axs[2].legend(fontsize=7)
    write_table("fig07_legacy_vs_padded", ["case", "quantity", "value"], rows)
    save(fig, "fig07_legacy_vs_padded")


def fig06_gaussianity(use_gs):
    grid, nreal, nreal_gs, nnodes = (16, 16, 8), 2000, 600, 400
    base = "f07"
    run_specsim(make_par(base, grid=grid, structs=((2, 1.0, 6, 3, 2, 45, 0, 0),), nreal=nreal, seed=4242), ("--bin",))
    fields = {"SPECSIM3D": load_bin(base, grid, nreal)}
    if use_gs:
        mod = gs_model(2, ranges=(6, 3, 2))
        x, y, z = (np.arange(n, dtype=float) for n in grid)
        fields["GSTools"] = np.array([gs.SRF(mod, seed=4242 + i, mode_no=1024).structured((x, y, z)) for i in range(nreal_gs)])
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.5))
    rows = []
    for label, f in fields.items():
        col = C_SPEC if label == "SPECSIM3D" else C_GS
        n = f.shape[0]
        v = f.reshape(n, -1)
        nodes = RNG.choice(v.shape[1], nnodes, replace=False)
        sk = stats.skew(v[:, nodes], axis=0); ku = stats.kurtosis(v[:, nodes], axis=0)
        sw = np.array([stats.shapiro(v[:, k]).pvalue for k in nodes])
        z = v[:, nodes[::20]].ravel(); z = (z - z.mean()) / z.std()
        osm, osr = stats.probplot(z, dist="norm", fit=False)
        axs[0].plot(osm[::7], osr[::7], ".", ms=2, color=col, label=f"{label} (n = {n})")
        axs[1].hist(sw, bins=20, range=(0, 1), histtype="step", color=col, lw=1.5, label=label)
        axs[2].hist(ku, bins=25, histtype="step", color=col, lw=1.5, label=f"{label}: mean {ku.mean():+.3f}")
        rows.append([label, n, nnodes, round(float(np.mean(sw < 0.05)), 4), round(float(sk.mean()), 4), round(float(ku.mean()), 4), round(float(v.var(axis=0).mean()), 4)])
    axs[0].plot([-4, 4], [-4, 4], "k-", lw=0.9); axs[0].set_title("Normal QQ plot (exponential model)"); axs[0].set_xlabel("theoretical quantiles"); axs[0].set_ylabel("observed quantiles"); axs[0].legend(fontsize=7)
    axs[1].set_title("Shapiro-Wilk p-value per node"); axs[1].set_xlabel("p"); axs[1].legend(fontsize=7)
    axs[2].set_title("Excess kurtosis per node"); axs[2].legend(fontsize=7)
    write_table("fig06_gaussianity", ["method", "n_realizations", "n_nodes_tested", "fraction_p_below_0.05", "mean_skewness", "mean_excess_kurtosis", "mean_node_variance"], rows)
    save(fig, "fig06_gaussianity")


def fig05_cases():
    grid = (64, 64, 16)
    cases = {"A iso 15/15/15": (0.0, [(1.0, 15, 15, 15, 0, 0, 0)]), "B 2:1 20/10/10": (0.0, [(1.0, 20, 10, 10, 30, 0, 0)]),
             "C 4:1 20/5/5": (0.0, [(1.0, 20, 5, 5, 60, 0, 0)]), "D layered 10/10/2": (0.0, [(1.0, 10, 10, 2, 0, 0, 0)]),
             "E extreme 40/4/2": (0.0, [(1.0, 40, 4, 2, 135, 0, 0)]), "F rotated 20/8/4 (az30 dip20 plg15)": (0.0, [(1.0, 20, 8, 4, 30, 20, 15)])}
    cand = np.column_stack([RNG.integers(-30, 31, 1500), RNG.integers(-30, 31, 1500), RNG.integers(-12, 13, 1500)])
    cand = cand[np.any(cand != 0, axis=1)]
    vp = make_vario_par("f08", [(a, b, c, 1) for a, b, c in cand])
    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    rows = []
    cols = {1: "tab:green", 2: "tab:blue", 3: "tab:orange"}
    nested = ("G nested + nugget 0.1 (sph + exp)", 0.1, [(1, 0.5, 20, 10, 5, 45, 0, 0), (2, 0.4, 8, 8, 3, 0, 0, 0)])
    jobs = []
    for cname, (nug, st) in cases.items():
        for mname, it in MODELS.items():
            jobs.append((cname, mname, it, nug, [(it,) + s for s in st]))
    jobs.append((nested[0], "nested", 0, nested[1], nested[2]))
    for cname, mname, it, nug, structs in jobs:
        base = "f08_run"
        run_specsim(make_par(base, grid=grid, nugget=nug, structs=tuple(structs), nreal=20, seed=2024), ("--bin",), vp)
        v = read_vario(base)
        gm = np.array([v[d]["mean"][0] for d in sorted(v) if v[d]["n"][0] > 0]); gt = np.array([v[d]["model"][0] for d in sorted(v) if v[d]["n"][0] > 0])
        se = np.array([v[d]["sd"][0] for d in sorted(v) if v[d]["n"][0] > 0]) / np.sqrt(20)
        tot = nug + sum(s[1] for s in structs)
        keep = (gt > 0.05 * tot) & (gt < 0.97 * tot)
        idx = np.flatnonzero(keep)[:300]
        ax.plot(gt[idx], gm[idx], ".", ms=2.2, alpha=0.5, color=cols.get(it, "tab:purple"))
        rms = float(np.sqrt(np.mean((gm[idx] - gt[idx]) ** 2))); noise = float(np.sqrt(np.mean(se[idx] ** 2)))
        rows.append([cname, mname, len(idx), round(rms, 4), round(noise, 4), round(rms / noise, 2), round(float(np.mean(gm[idx] - gt[idx])), 4)])
    ax.plot([0, 1.05], [0, 1.05], "k-", lw=0.9)
    for mname, it in MODELS.items():
        ax.plot([], [], ".", color=cols[it], label=mname)
    ax.plot([], [], ".", color="tab:purple", label="nested + nugget")
    ax.set_xlabel("model semivariogram"); ax.set_ylabel("mean experimental semivariogram (20 realizations)")
    ax.set_title("300 arbitrary lag vectors per case, 6 anisotropy cases + nested"); ax.legend(fontsize=7)
    write_table("fig05_cases", ["case", "model", "n_lags", "rms", "mc_noise", "rms_over_noise", "bias"], rows)
    save(fig, "fig05_anisotropy_cases")


def peak_rss_run(cmd, cwd):
    """Run cmd, return (stdout, wall seconds, peak RSS in MB or None)."""
    t0 = time.time()
    if os.name != "nt":
        wrapper = ("import subprocess,resource,sys;p=subprocess.run(sys.argv[1:],capture_output=True,text=True);"
                   "print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss);sys.stdout.write(p.stdout)")
        p = subprocess.run([sys.executable, "-I", "-c", wrapper] + cmd, cwd=cwd, capture_output=True, text=True)
        first, _, rest = p.stdout.partition("\n")
        return rest, time.time() - t0, int(first) / 1024.0
    try:
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, text=True)
        try:
            import psutil
            ps, peak = psutil.Process(proc.pid), 0
        except Exception:
            ps, peak = None, 0
        while proc.poll() is None:
            if ps is not None:
                try:
                    mi = ps.memory_info()
                    peak = max(peak, getattr(mi, "peak_wset", 0), mi.rss)
                except Exception:
                    ps = None
            time.sleep(0.05)
        out = proc.stdout.read()
        if peak == 0:                                  # no psutil: ask Windows for the peak working set of the finished process
            peak = _win_peak_wset(proc)
        return out, time.time() - t0, (peak / 2**20 if peak else None)
    except Exception:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
        return p.stdout, time.time() - t0, None


def _win_peak_wset(proc):
    """Peak working set (bytes) of a finished child process on Windows, with ctypes (no extra package)."""
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
        pmc = PMC(); pmc.cb = ctypes.sizeof(PMC)
        ok = ctypes.windll.psapi.GetProcessMemoryInfo(wintypes.HANDLE(int(proc._handle)), ctypes.byref(pmc), pmc.cb)
        return pmc.PeakWorkingSetSize if ok else 0
    except Exception:
        return 0


def measure_scaling(use_gs, sizes):
    res = {}
    for n in sizes:
        nreal = 1 if n >= 256 else 2
        base = f"sc_{n}"
        par = make_par(base, grid=(n, n, n), nreal=nreal, seed=99)
        out, wall, rss = peak_rss_run([EXE, par, "--bin"], WORK)
        grid = [int(x) for x in re.search(r"periodic grid:\s+(\d+)\s+(\d+)\s+(\d+)", out).groups()]
        t_sim = float(re.search(r"realizations\s+:\s+([\d.]+)", out).group(1))
        t_spec = float(re.search(r"covariance\+spectrum:\s+([\d.]+)", out).group(1))
        f = load_bin(base, (n, n, n), nreal)[0:1]
        errs = []
        for step, a in [((1, 1, 0), RANGES[0]), ((1, -1, 0), RANGES[1]), ((0, 0, 1), RANGES[2])]:
            g, l = vario_shift(f, step, int(a))
            errs.append(np.sqrt(np.mean((g[0] - (1 - corr(1, l / a))) ** 2)))
        res[f"specsim3d_{n}"] = dict(n=n, periodic=grid, wall_s=wall, spectrum_s=t_spec, per_realization_s=t_sim / nreal, rss_mb=rss, rms=float(np.mean(errs)), nreal=nreal)
        os.remove(os.path.join(WORK, base + ".bin"))
        print(f"    SPECSIM3D n={n}: {t_sim / nreal:.2f} s/realization", flush=True)
    if use_gs:
        for n in [s for s in sizes if s <= 128]:
            mod = gs_model(1)
            x = np.arange(n, dtype=float)
            t0 = time.time()
            f = gs.SRF(mod, seed=99, mode_no=1024).structured((x, x, x))
            gs.SRF(mod, seed=100, mode_no=1024).structured((x, x, x))
            tt = (time.time() - t0) / 2
            errs = []
            for step, a in [((1, 1, 0), RANGES[0]), ((1, -1, 0), RANGES[1]), ((0, 0, 1), RANGES[2])]:
                g, l = vario_shift(f[None], step, int(a))
                errs.append(np.sqrt(np.mean((g[0] - (1 - corr(1, l / a))) ** 2)))
            res[f"gstools_{n}"] = dict(n=n, per_realization_s=tt, rms=float(np.mean(errs)))
            print(f"    GSTools   n={n}: {tt:.2f} s/realization", flush=True)
    return res


def fig09_scaling(use_gs, full):
    cache = os.path.join(FIG, "scaling_cached.json")
    if full:
        print("  (re-measuring scaling; this takes a few minutes)", flush=True)
        data = measure_scaling(use_gs, [32, 64, 128, 256])
        data["system"] = SYS
        json.dump(data, open(cache, "w"), indent=1)
    else:
        data = json.load(open(cache))
        print("  (using cached scaling measurements in figures/scaling_cached.json; --full re-measures on this machine)")
    sysinfo = data.get("system", {"cpu": "see cached note", "cores": "?", "omp_threads": "?", "ram_gb": "?", "compiler": "?"})
    n_s = sorted(int(k.split("_")[1]) for k in data if k.startswith("specsim3d_"))
    n_g = sorted(int(k.split("_")[1]) for k in data if k.startswith("gstools_"))
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.5))
    ns_t = [n for n in n_s if data[f"specsim3d_{n}"]["per_realization_s"] > 0]      # a time of 0 (below the 0.01 s resolution of the timer) cannot be drawn on a log axis
    axs[0].loglog(ns_t, [data[f"specsim3d_{n}"]["per_realization_s"] for n in ns_t], "o-", color=C_SPEC, label="SPECSIM3D")
    if n_g:
        axs[0].loglog(n_g, [data[f"gstools_{n}"]["per_realization_s"] for n in n_g], "s-", color=C_GS, label="GSTools SRF (1024 modes)")
    axs[0].set_xlabel("grid side n (n$^3$ nodes)"); axs[0].set_ylabel("seconds per realization"); axs[0].legend(fontsize=7)
    mem = [data[f"specsim3d_{n}"].get("rss_mb") for n in n_s]
    if all(m is not None for m in mem):
        axs[1].loglog(n_s, mem, "o-", color=C_SPEC)
    axs[1].set_xlabel("grid side n"); axs[1].set_ylabel("peak memory (MB)")
    axs[2].loglog(n_s, [data[f"specsim3d_{n}"]["rms"] for n in n_s], "o-", color=C_SPEC, label="SPECSIM3D")
    if n_g:
        axs[2].loglog(n_g, [data[f"gstools_{n}"]["rms"] for n in n_g], "s-", color=C_GS, label="GSTools (spherical)")
    axs[2].set_xlabel("grid side n"); axs[2].set_ylabel("RMS variogram error\n(1 realization, 3 axes)"); axs[2].legend(fontsize=7)
    for a in axs:                                       # plain tick labels at the simulated grid sides (no overlapping minor labels)
        a.set_xticks(n_s); a.set_xticklabels([str(n) for n in n_s]); a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    fig.subplots_adjust(wspace=0.4)
    write_table("fig09_scaling", ["n", "periodic_grid", "spectrum_s", "s_per_realization", "peak_rss_mb", "rms_one_realization", "gstools_s_per_realization", "gstools_rms"],
                [[n, "x".join(map(str, data[f"specsim3d_{n}"]["periodic"])), round(data[f"specsim3d_{n}"]["spectrum_s"], 2), round(data[f"specsim3d_{n}"]["per_realization_s"], 3),
                  None if data[f"specsim3d_{n}"].get("rss_mb") is None else round(data[f"specsim3d_{n}"]["rss_mb"]), round(data[f"specsim3d_{n}"]["rms"], 4),
                  round(data[f"gstools_{n}"]["per_realization_s"], 2) if f"gstools_{n}" in data else "", round(data[f"gstools_{n}"]["rms"], 3) if f"gstools_{n}" in data else ""] for n in n_s])
    save(fig, "fig09_scaling")


def fig10_vario_speed(full):
    sizes = [64, 128] + ([256] if full else [])
    vp = make_vario_par("f10", [(1, 1, 0, 40), (1, -1, 0, 40), (0, 0, 1, 40)])
    rows = []
    for n in sizes:
        base = f"f10_{n}"
        nreal = {64: 20, 128: 6}.get(n, 1)              # several realizations: the timer of the program has 0.01 s resolution
        out, _ = run_specsim(make_par(base, grid=(n, n, n), nreal=nreal, seed=5), ("--bin",), vp)
        t_f = float(re.search(r"time variograms\s+:\s+([\d.]+)", out).group(1)) / nreal
        f = load_bin(base, (n, n, n), nreal)
        t0 = time.time()
        for step in [(1, 1, 0), (1, -1, 0), (0, 0, 1)]:
            vario_shift(f, step, 40)
        t_p = (time.time() - t0) / nreal
        rows.append([n, round(t_f, 4), round(t_p, 4), round(t_p / max(t_f, 1e-4), 1)])
        os.remove(os.path.join(WORK, base + ".bin"))
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    ax.loglog([r[0] for r in rows], [r[1] for r in rows], "o-", color=C_SPEC, label="Fortran (--vario, OpenMP)")
    ax.loglog([r[0] for r in rows], [r[2] for r in rows], "s-", color="0.35", label="numpy (array slicing)")
    ax.set_xlabel("grid side n (n$^3$ nodes)"); ax.set_ylabel("seconds per realization\n(3 directions x 40 lags)"); ax.legend(fontsize=7)
    ax.set_ylim(top=max(r[2] for r in rows) * 6)                  # room for the speed-up labels above the highest point
    ax.set_xticks([r[0] for r in rows]); ax.set_xticklabels([str(r[0]) for r in rows]); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    for r in rows:
        ax.annotate(f"x{r[3]:.0f}", (r[0], r[2]), textcoords="offset points", xytext=(0, 5), ha="center", fontsize=7)
    write_table("fig10_variogram_speed", ["n", "fortran_s_per_realization", "numpy_s_per_realization", "speedup"], rows)
    footer(fig)
    save(fig, "fig10_variogram_speed")


def fig08_validation_2d():
    if not os.path.exists(ORIG2D):
        print("  skipped (docs/data/specsim2_original_output.out not found)")
        return
    d = np.loadtxt(ORIG2D)
    z = d[:, 2].reshape(128, 128)                                    # x slow, y fast (SPECSIM2 ordering)
    def vs(f, step, kmax):
        out = []
        for k in range(1, kmax + 1):
            dx, dy = k * step[0], k * step[1]
            a = f[max(0, -dx):f.shape[0] - max(0, dx), max(0, -dy):f.shape[1] - max(0, dy)]
            b = f[max(0, dx):f.shape[0] - max(0, -dx), max(0, dy):f.shape[1] - max(0, -dy)]
            out.append(0.5 * np.mean((a - b) ** 2))
        return np.array(out)
    base = "f11"
    structs = ((2, 0.7, 30.0, 30.0, 30.0, 45.0, 0, 0),)
    run_specsim(make_par(base, grid=(128, 128, 1), nugget=0.3, structs=structs, nreal=20, seed=12345, ipad=1), ("--bin",), make_vario_par(base, [(1, 0, 0, 40), (0, 1, 0, 40)], iper=1))
    v = read_vario(base)
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.3), sharey=True)
    h = np.arange(1, 41)
    gmod = 0.3 + 0.7 * (1 - np.exp(-h / 10.0))
    rows = []
    for k, (ax, name, step) in enumerate(zip(axs, ["x direction", "y direction"], [(1, 0), (0, 1)]), start=1):
        band(ax, v[k]["h"], v[k]["per"].T, C_SPEC, alpha=0.2, label="SPECSIM3D 5th-95th percentile")
        ax.plot(v[k]["h"], v[k]["mean"], color=C_SPEC, lw=1.5, label="SPECSIM3D, nz = 1 (mean of 20)")
        g0 = vs(z, step, 40)
        ax.plot(h, g0, color=C_GS, lw=1.2, label="SPECSIM2.EXE (1 realization)")
        ax.plot(h, gmod, color=C_MOD, lw=2, label="model: 0.3 + 0.7 (1 - exp(-h/10))")
        ax.set_title(name); ax.set_xlabel("lag distance h")
        rows.append([name, round(float(np.sqrt(np.mean((v[k]["mean"] - gmod) ** 2))), 4), round(float(np.sqrt(np.mean((g0 - gmod) ** 2))), 4)])
    axs[0].set_ylabel("semivariogram"); axs[0].legend(fontsize=7, loc="lower right")
    fig.suptitle("2D check against the original SPECSIM2 output (exponential, nugget 0.3, sill 0.7, 128 x 128)", fontsize=9.5, y=1.03)
    write_table("fig08_validation_2d", ["direction", "rms_specsim3d_vs_model", "rms_specsim2_vs_model"], rows)
    save(fig, "fig08_validation_2d")


# ============================================================================ conditional simulation
CSTRUCT = ((2, 1.0, 20.0, 10.0, 5.0, 45.0, 0, 0),)      # exponential, practical ranges 20/10/5, N45E (Yao geometry)
CGRID = (64, 64, 16)


def write_gslib_points(name, pts, names=("x", "y", "z", "value")):
    path = os.path.join(WORK, name)
    with open(path, "w") as fh:
        fh.write("conditioning data\n%d\n" % len(names) + "\n".join(names) + "\n")
        for p in pts:
            fh.write(" ".join("%.8g" % v for v in p) + "\n")
    return path


def make_cond_par(name, datafile, ncol=4, izc=3, itrans=0, skmean=0.0, zlim=(0.0, 100.0), ltail=(1, 0.0),
                  utail=(1, 0.0), radii=(30.0, 15.0, 8.0), angles=(45.0, 0.0, 0.0), ndmax=24, iwt=0):
    lines = [os.path.basename(datafile), f"1 2 {izc} 4 {iwt}", "-1.0e21 1.0e21", f"{itrans} {skmean}",
             f"{zlim[0]} {zlim[1]}", f"{ltail[0]} {ltail[1]}", f"{utail[0]} {utail[1]}",
             " ".join(map(str, radii)), " ".join(map(str, angles)), str(ndmax)]
    path = os.path.join(WORK, name)
    open(path, "w").write("\n".join(lines) + "\n")
    return path


def run_cond(base, grid, structs, nreal, seed, cfile, nugget=0.0):
    par = make_par(base, grid=grid, structs=structs, nreal=nreal, seed=seed, ipad=1, nugget=nugget)
    out, t = run_specsim(par, ("--bin", "--cond", os.path.basename(cfile)))
    m = re.search(r"time conditioning:\s+([\d.]+) s", out)
    tc = float(m.group(1)) if m else float("nan")
    return load_bin(base, grid, nreal), t, tc, out


def sk_reference(grid, structs, xyz, val, nugget=0.0):
    """Exact simple kriging (mean 0, all data) with the GSLIB-type covariance of SPECSIM3D (azimuth rotation only)."""
    it, sill, a1, a2, a3, az = structs[0][:6]
    c0 = nugget + sum(s[1] for s in structs)
    azr = np.radians(az)
    u = np.array([np.sin(azr), np.cos(azr), 0.0]); w = np.array([np.cos(azr), -np.sin(azr), 0.0])

    def C(d):
        h = np.sqrt((d @ u / a1) ** 2 + (d @ w / a2) ** 2 + (d[..., 2] / a3) ** 2)
        return sill * corr(it, h)
    K = C(xyz[:, None, :] - xyz[None, :, :]) + nugget * np.eye(len(xyz))
    gx, gy, gz = np.meshgrid(*(np.arange(n, dtype=float) for n in grid), indexing="ij")
    P = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]
    k = C(xyz[None, :, :] - P[:, None, :])
    W = np.linalg.solve(K + 1e-10 * np.eye(len(xyz)), k.T).T
    est = (W @ val).reshape(grid)
    var = (c0 - (W * k).sum(1)).reshape(grid)
    return est, var


def gs_cond(grid, xyz, val, nreal, seed, it=2):
    mod = gs_model(it)
    kr = gs.krige.Simple(mod, cond_pos=tuple(xyz.T), cond_val=val, mean=0.0)
    cs = gs.CondSRF(kr, seed=seed, mode_no=1024)
    ax = tuple(np.arange(n, dtype=float) for n in grid)
    t0 = time.time()
    f = []
    for i in range(nreal):
        f.append(cs(ax, seed=seed + i, mesh_type="structured"))
    return np.array(f), (time.time() - t0) / nreal


def cond_case(nd=150, nreal=200):
    """Truth realization, conditioning data (nd random nodes), SPECSIM3D conditional realizations and exact SK."""
    run_specsim(make_par("ctruth", grid=CGRID, structs=CSTRUCT, nreal=1, seed=2024), ("--bin",))
    truth = load_bin("ctruth", CGRID, 1)[0]
    rng = np.random.default_rng(7)
    idx = rng.choice(np.prod(CGRID), nd, replace=False)
    ijk = np.array(np.unravel_index(idx, CGRID)).T
    xyz = ijk.astype(float)
    val = truth[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
    df = write_gslib_points("cdata.dat", np.c_[xyz, val])
    cf = make_cond_par("ccond.cpar", df)
    sims, t, tc, out = run_cond("ccond", CGRID, CSTRUCT, nreal, 777, cf)
    est, var = sk_reference(CGRID, CSTRUCT, xyz, val)
    return dict(truth=truth, xyz=xyz, val=val, sims=sims, est=est, var=var, t=t, tc=tc, nreal=nreal)


_CC = {}


def get_cc():
    if "c" not in _CC:
        _CC["c"] = cond_case()
    return _CC["c"]


def fig11_cond_maps():
    cc = get_cc()
    k = 8
    sims, est = cc["sims"], cc["est"]
    inl = np.abs(cc["xyz"][:, 2] - k) < 0.5
    fig = plt.figure(figsize=(10.6, 3.0))
    gsp = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 0.05, 1, 0.05], wspace=0.55)
    axs = [fig.add_subplot(gsp[0, 0]), fig.add_subplot(gsp[0, 1]), fig.add_subplot(gsp[0, 2]), None, fig.add_subplot(gsp[0, 4])]
    cax1, cax2 = fig.add_subplot(gsp[0, 3]), fig.add_subplot(gsp[0, 5])
    ext = (-0.5, CGRID[0] - 0.5, -0.5, CGRID[1] - 0.5)
    panels = [("(a) exact simple kriging", est[:, :, k]), ("(b) mean of %d realizations" % cc["nreal"], sims.mean(0)[:, :, k]),
              ("(c) one realization", sims[0][:, :, k])]
    for ax, (t, a) in zip(axs[:3], panels):
        im = ax.imshow(a.T, origin="lower", extent=ext, cmap="RdBu_r", vmin=-3, vmax=3)
        ax.scatter(cc["xyz"][inl, 0], cc["xyz"][inl, 1], c=cc["val"][inl], cmap="RdBu_r", vmin=-3, vmax=3,
                   edgecolors="k", linewidths=0.8, s=26)
        ax.set_title(t, fontsize=8.5); ax.set_xlabel("x"); ax.grid(False)
    axs[0].set_ylabel("y")
    cb1 = fig.colorbar(im, cax=cax1); cb1.set_label("value"); cb1.ax.yaxis.set_label_position("left")
    im2 = axs[4].imshow(sims.std(0, ddof=1)[:, :, k].T, origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=1)
    axs[4].scatter(cc["xyz"][inl, 0], cc["xyz"][inl, 1], c="w", edgecolors="k", linewidths=0.8, s=14)
    axs[4].set_title("(d) std. deviation of realizations", fontsize=8.5); axs[4].set_xlabel("x"); axs[4].grid(False)
    axs[4].tick_params(labelleft=False)
    fig.colorbar(im2, cax=cax2).set_label("std. deviation")
    fig.suptitle("Conditional simulation, horizontal slice z = 8 (circles: conditioning data within the layer)", fontsize=9.5, y=1.02)
    save(fig, "fig11_conditional_maps")


def fig12_cond_validation(use_gs):
    cc = get_cc()
    sims, est, var, nr = cc["sims"], cc["est"], cc["var"], cc["nreal"]
    sets = [("SPECSIM3D", sims, C_SPEC, cc["t"] / nr)]
    if use_gs:
        g, tg = gs_cond(CGRID, cc["xyz"], cc["val"], 100, 4242)
        sets.append(("GSTools CondSRF", g, C_GS, tg))
    fig, axs = plt.subplots(2, len(sets), figsize=(3.6 * len(sets) + 0.4, 6.2), squeeze=False)
    rows = []
    rng = np.random.default_rng(3)
    sub = rng.choice(est.size, 4000, replace=False)
    for j, (name, f, col, tt) in enumerate(sets):
        n = len(f)
        m = f.mean(0); v = f.var(0, ddof=1)
        se = np.sqrt(var / n)
        zs = (m - est) / np.maximum(se, 1e-9)
        free = var > 1e-6
        ax = axs[0, j]
        ax.plot(est.ravel()[sub], m.ravel()[sub], ".", ms=2, color=col, alpha=0.5)
        lim = [est.min(), est.max()]; ax.plot(lim, lim, "k-", lw=1)
        ax.set_xlabel("exact simple kriging estimate"); ax.set_ylabel("mean of realizations"); ax.set_title(f"{name} ({n} realizations)", fontsize=9)
        ax = axs[1, j]
        ax.plot(var.ravel()[sub], v.ravel()[sub], ".", ms=2, color=col, alpha=0.5)
        lim = [0, 1.05]; ax.plot(lim, lim, "k-", lw=1)
        ax.set_xlabel("exact simple kriging variance"); ax.set_ylabel("variance of realizations")
        rms_m = float(np.sqrt(np.mean((m - est) ** 2)))
        rms_se = float(np.sqrt(np.mean(se ** 2)))
        rms_v = float(np.sqrt(np.mean((v - var) ** 2)))
        rms_vse = float(np.sqrt(np.mean(2 * var ** 2 / (n - 1))))       # sampling error of a variance for Gaussian data
        slope = float(np.sum(v * var) / np.sum(var * var))
        zstd = float(np.std(zs[free]))
        ii = cc["xyz"].astype(int)
        dm = float(abs(f[:, ii[:, 0], ii[:, 1], ii[:, 2]] - cc["val"]).max())
        rows.append([name, n, round(rms_m, 4), round(rms_se, 4), round(rms_v, 4), round(rms_vse, 4), round(slope, 4), round(zstd, 3),
                     "%.2e" % dm, round(tt, 3)])
        axs[0, j].text(0.04, 0.93, f"rms diff = {rms_m:.3f}\n(Monte Carlo error {rms_se:.3f})", transform=axs[0, j].transAxes, fontsize=7.5, va="top")
        axs[1, j].text(0.04, 0.93, f"rms diff = {rms_v:.3f}\n(Monte Carlo error {rms_vse:.3f})\nregression slope = {slope:.3f}",
                       transform=axs[1, j].transAxes, fontsize=7.5, va="top")
    write_table("fig12_conditional_validation", ["code", "nreal", "rms_mean_diff", "mc_error_mean", "rms_var_diff", "mc_error_var",
                                                 "var_slope", "std_of_standardized_mean_diff", "max_abs_data_error", "seconds_per_realization"], rows)
    fig.suptitle("Conditional realizations against exact simple kriging (64 x 64 x 16 grid, 150 data)", fontsize=9.5, y=1.0)
    fig.tight_layout()
    save(fig, "fig12_conditional_validation")


def fig13_cond_nscore():
    cc = get_cc()
    truth = np.exp(1.0 + 0.8 * cc["truth"])                      # lognormal "truth" with the same spatial structure
    xyz = cc["xyz"]
    ijk = xyz.astype(int)
    val = truth[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
    df = write_gslib_points("cdata_ln.dat", np.c_[xyz, val])
    zhi = float(val.max() * 2)
    cf = make_cond_par("ccond_ln.cpar", df, itrans=1, zlim=(0.0, zhi), ltail=(1, 0.0), utail=(1, 0.0))
    sims, t, tc, out = run_cond("ccond_ln", CGRID, CSTRUCT, 50, 99, cf)
    dmax = float(abs(sims[:, ijk[:, 0], ijk[:, 1], ijk[:, 2]] - val).max() / val.max())
    fig = plt.figure(figsize=(10.2, 3.1))
    ax = fig.add_subplot(1, 3, 1)
    for a, lab, c, ls in [(truth.ravel(), "truth (all nodes)", "k", "-"), (val, f"{len(val)} data", C_GS, "--"),
                          (sims.ravel(), "conditional realizations", C_SPEC, "-")]:
        s = np.sort(a); ax.plot(s, (np.arange(len(s)) + 0.5) / len(s), color=c, ls=ls, lw=1.4, label=lab)
    ax.set_xscale("log"); ax.set_xlim(0.1, 100); ax.set_xlabel("value"); ax.set_ylabel("cumulative probability"); ax.legend(fontsize=7, loc="upper left")
    ax.set_title("(a) marginal distribution", fontsize=9)
    k = 8
    inl = np.abs(xyz[:, 2] - k) < 0.5
    vmax = float(np.percentile(truth, 99))
    for n, (t_, a) in enumerate([("(b) truth, z = 8", truth[:, :, k]), ("(c) one realization, z = 8", sims[0][:, :, k])], start=2):
        ax = fig.add_subplot(1, 3, n)
        im = ax.imshow(a.T, origin="lower", extent=(-0.5, 63.5, -0.5, 63.5), cmap="viridis", vmin=0, vmax=vmax)
        ax.scatter(xyz[inl, 0], xyz[inl, 1], c=val[inl], cmap="viridis", vmin=0, vmax=vmax, edgecolors="w", linewidths=0.8, s=26)
        ax.set_title(t_, fontsize=9); ax.set_xlabel("x"); ax.grid(False)
        if n == 3:
            ax.set_ylabel("y")
    fig.colorbar(im, ax=fig.axes[1:], fraction=0.02, pad=0.01).set_label("value")
    fig.suptitle("Normal score transform, conditioning and back transform (lognormal data)", fontsize=9.5, y=1.03)
    write_table("fig13_conditional_nscore", ["quantity", "value"],
                [["max |realization - data| / max data at data nodes", "%.2e" % dmax],
                 ["mean truth", round(float(truth.mean()), 3)], ["mean data", round(float(val.mean()), 3)], ["mean realizations", round(float(sims.mean()), 3)],
                 ["variance truth", round(float(truth.var()), 3)], ["variance data", round(float(val.var()), 3)], ["variance realizations", round(float(sims.var()), 3)]])
    save(fig, "fig13_conditional_nscore")


def fig14_cond_speed(use_gs):
    rng = np.random.default_rng(11)
    sizes = [(32, 32, 8), (64, 64, 16), (128, 128, 32), (256, 256, 64)]
    rows = []
    for g in sizes:
        n = int(np.prod(g)); nd = max(20, n // 500)
        xyz = np.c_[rng.uniform(0, g[0] - 1, nd), rng.uniform(0, g[1] - 1, nd), rng.uniform(0, g[2] - 1, nd)]
        df = write_gslib_points("cspeed.dat", np.c_[xyz, rng.normal(size=nd)])
        cf = make_cond_par("cspeed.cpar", df, ndmax=16, radii=(20.0, 10.0, 5.0))
        par0 = make_par("cs0", grid=g, structs=CSTRUCT, nreal=10, seed=5)
        _, t_unc = run_specsim(par0, ("--bin",))
        out, t_c = run_specsim(make_par("cs1", grid=g, structs=CSTRUCT, nreal=10, seed=5), ("--bin", "--cond", "cspeed.cpar"))
        tc = float(re.search(r"time conditioning:\s+([\d.]+) s", out).group(1))
        tgs = float("nan")
        if use_gs and n <= 64 * 64 * 16:
            _, tgs = gs_cond(g, xyz, rng.normal(size=nd), 3, 1)
            tgs *= 10                                                     # extrapolated to 10 realizations
        rows.append([n, nd, round(t_unc, 2), round(t_c, 2), round(tc, 2), round(tgs, 2)])
    nds = [4, 8, 16, 32, 64]
    g = (128, 128, 32); n = int(np.prod(g)); nd = n // 250
    xyz = np.c_[rng.uniform(0, g[0] - 1, nd), rng.uniform(0, g[1] - 1, nd), rng.uniform(0, g[2] - 1, nd)]
    df = write_gslib_points("cspeed2.dat", np.c_[xyz, rng.normal(size=nd)])
    rows2 = []
    for m in nds:
        cf = make_cond_par("cspeed2.cpar", df, ndmax=m, radii=(60.0, 60.0, 30.0))
        tt = []
        for _rep in range(3):
            out, _ = run_specsim(make_par("cs2", grid=g, structs=CSTRUCT, nreal=10, seed=5), ("--bin", "--cond", "cspeed2.cpar"))
            tt.append(float(re.search(r"time conditioning:\s+([\d.]+) s", out).group(1)))
        rows2.append([m, min(tt)])
    fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.2))
    r = np.array(rows, float)
    axs[0].loglog(r[:, 0], r[:, 2], "s--", color="0.4", label="SPECSIM3D, unconditional")
    axs[0].loglog(r[:, 0], r[:, 3], "o-", color=C_SPEC, label="SPECSIM3D, conditional (total)")
    axs[0].loglog(r[:, 0], r[:, 4], "^:", color=C_SPEC, label="of which conditioning")
    ok = ~np.isnan(r[:, 5])
    if ok.any():
        axs[0].loglog(r[ok, 0], r[ok, 5], "d-", color=C_GS, label="GSTools CondSRF (10 realizations)")
    axs[0].set_xlabel("number of grid nodes (1 datum per 500 nodes)"); axs[0].set_ylabel("time for 10 realizations (s)")
    axs[0].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, frameon=False); axs[0].set_title("(a) grid size, 16 neighbours", fontsize=9)
    r2 = np.array(rows2)
    axs[1].plot(r2[:, 0], r2[:, 1], "o-", color=C_SPEC)
    axs[1].set_xlabel("maximum number of neighbours (ndmax)"); axs[1].set_ylabel("conditioning time (s)")
    axs[1].set_title("(b) 128 x 128 x 32 grid, %d data, best of 3" % nd, fontsize=9)
    axs[0].set_ylim(top=200)
    fig.tight_layout(); fig.subplots_adjust(bottom=0.36)
    footer(fig)
    write_table("fig14_conditional_speed_grid", ["nodes", "ndata", "t_uncond_s", "t_cond_total_s", "t_conditioning_s", "t_gstools_cond_10real_s"], rows)
    write_table("fig14_conditional_speed_ndmax", ["ndmax", "t_conditioning_s"], rows2)
    save(fig, "fig14_conditional_speed")


# ============================================================================ main

def fig15_thread_scaling(full):
    """Strong scaling with the number of OpenMP threads (unconditional simulation, spherical model, 20/10/5)."""
    ncpu = os.cpu_count() or 1
    threads = sorted({t for t in (1, 2, 4, 8, 16, 32, 64) if t <= ncpu} | {ncpu})
    sizes = [(64, 3), (128, 3), (256, 1)] if full else [(64, 3), (128, 2)]     # (output grid n, realizations)
    rows, res = [], {}
    for n, nreal in sizes:
        for th in threads:
            base = f"th_{n}_{th}"
            par = make_par(base, grid=(n, n, n), nreal=nreal, seed=99)
            best = None
            for rep in range(2 if n < 256 else 1):                                # best of 2 against OS noise
                env = dict(os.environ, OMP_NUM_THREADS=str(th))
                t0 = time.time()
                out = subprocess.run([EXE, par, "--bin"], cwd=WORK, capture_output=True, text=True, env=env).stdout
                wall = time.time() - t0
                t_sim = float(re.search(r"realizations\s+:\s+([\d.]+)", out).group(1)) / nreal
                t_spec = float(re.search(r"covariance\+spectrum:\s+([\d.]+)", out).group(1))
                if best is None or t_sim < best[0]:
                    best = (t_sim, t_spec, wall)
            for ext in (".bin",):
                try:
                    os.remove(os.path.join(WORK, base + ext))
                except OSError:
                    pass
            res[(n, th)] = best
            print(f"    n={n} threads={th}: {best[0]:.3f} s/realization, spectrum {best[1]:.2f} s", flush=True)
    for n, _ in sizes:
        t1 = res[(n, 1)][0]; s1 = res[(n, 1)][1]
        for th in threads:
            t, sp, wall = res[(n, th)]
            rows.append([n, 2 * n, th, round(sp, 3), round(t, 4), round(t1 / t, 2), round(t1 / t / th, 2), round(s1 / sp, 2)])
    write_table("fig15_thread_scaling", ["n", "periodic_n", "threads", "spectrum_s", "s_per_realization", "speedup", "efficiency", "spectrum_speedup"], rows)
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.5))
    cols = {64: "tab:green", 128: C_SPEC, 256: "tab:purple"}
    for n, _ in sizes:
        t = np.array([res[(n, th)][0] for th in threads]); sp = np.array([res[(n, th)][1] for th in threads])
        axs[0].plot(threads, t, "o-", color=cols[n], label=f"{n}$^3$ ({2 * n}$^3$ periodic)")
        axs[1].plot(threads, t[0] / t, "o-", color=cols[n], label=f"{n}$^3$")
        axs[2].plot(threads, sp[0] / sp, "o-", color=cols[n], label=f"{n}$^3$")
    for ax, ttl in zip(axs, ["Time per realization", "Speed-up of the realizations", "Speed-up of the spectrum step"]):
        ax.set_xscale("log", base=2); ax.set_xticks(threads); ax.set_xticklabels([str(t) for t in threads]); ax.set_title(ttl)
        ax.set_xlabel("OpenMP threads")
    axs[0].set_yscale("log"); axs[0].set_ylabel("s per realization")
    for ax in axs[1:]:
        ax.plot(threads, threads, "k--", lw=0.9, label="ideal"); ax.set_ylabel("speed-up vs 1 thread")
    axs[1].legend(fontsize=7); axs[0].legend(fontsize=7)
    footer(fig, f"{ncpu} logical cores available; best of 2 runs (1 for 256$^3$)")
    save(fig, "fig15_thread_scaling")

FIGS = {"fig01": lambda a: fig01_geometry(), "fig02": lambda a: fig02_algorithm(), "fig03": lambda a: fig03_slices(),
        "fig04": lambda a: fig04_models(a.use_gs), "fig05": lambda a: fig05_cases(), "fig06": lambda a: fig06_gaussianity(a.use_gs),
        "fig07": lambda a: fig07_legacy_vs_padded(), "fig08": lambda a: fig08_validation_2d(), "fig09": lambda a: fig09_scaling(a.use_gs, a.full),
        "fig10": lambda a: fig10_vario_speed(a.full), "figS1": lambda a: figS1_variograms_yao(),
        "fig11": lambda a: fig11_cond_maps(), "fig12": lambda a: fig12_cond_validation(a.use_gs),
        "fig13": lambda a: fig13_cond_nscore(), "fig14": lambda a: fig14_cond_speed(a.use_gs),
        "fig15": lambda a: fig15_thread_scaling(a.full)}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-gstools", action="store_true", help="skip the GSTools comparisons")
    ap.add_argument("--full", action="store_true", help="re-measure scaling and variogram speed up to 256^3 (minutes, ~3 GB RAM)")
    ap.add_argument("--only", nargs="*", help="figure ids to produce (e.g. fig07 fig05; figS1 = supplementary figure of the manual)")
    args = ap.parse_args()
    args.use_gs = HAVE_GS and not args.no_gstools
    if not os.path.exists(EXE):
        sys.exit(f"Executable not found: {EXE}\nBuild it first (build_specsim3d.bat on Windows, or: gfortran -O3 -fopenmp -o specsim3d specsim3d.f90).")
    if not HAVE_SCIPY:
        print("warning: scipy not available, fig06 will fail")
    os.makedirs(WORK, exist_ok=True); os.makedirs(TAB, exist_ok=True)
    SYS = system_info()
    json.dump(SYS, open(os.path.join(FIG, "system_info.json"), "w"), indent=1)
    print("System:", json.dumps(SYS))
    print("GSTools comparisons:", "ON" if args.use_gs else "OFF")
    todo = args.only or list(FIGS)
    t00 = time.time()
    for key in todo:
        t0 = time.time()
        print(f"[{key}]", flush=True)
        FIGS[key](args)
        print(f"  done in {time.time() - t0:.1f} s", flush=True)
    shutil.rmtree(WORK, ignore_errors=True)
    print(f"All done in {time.time() - t00:.0f} s. Figures in {FIG}")
