"""GSTools test: anisotropic 3D simulation (Yao 1998 case) and variogram analysis.

Model: spherical, variance 1, ranges 20/10/5 (major N45E, intermediate N45W, vertical),
64 x 64 x 16 grid, unit spacing, seed 31791191, 20 realizations.
Autore: Stefano Bonduà
"""
import time
import numpy as np
import gstools as gs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os as _os
OUT = _os.path.dirname(_os.path.abspath(__file__))
RES = _os.path.join(OUT, "results"); _os.makedirs(RES, exist_ok=True)
NX, NY, NZ = 64, 64, 16
NREAL = 20
SEED = 31791191
MODE_NO = 1024

# --- model -------------------------------------------------------------
# angles: rotation in the xy plane (counterclockwise from x = East). N45E with x = East, y = North
# is 45 degrees from the x axis.
model = gs.Spherical(dim=3, var=1.0, len_scale=[20.0, 10.0, 5.0], angles=[np.pi / 4, 0.0, 0.0])
print("len_scale:", model.len_scale, "| anis:", model.anis, "| integral_scale:", model.integral_scale)

# check that the GSTools len_scale = range (the spherical variogram reaches the sill at r = range)
r_test = np.array([10.0, 20.0, 25.0])
print("variogram along axis 0 at r=10,20,25:", model.vario_axis(r_test, axis=0))
sph = lambda h, a: np.where(h < a, 1.5 * h / a - 0.5 * (h / a) ** 3, 1.0)
print("theoretical spherical (range 20):     ", sph(r_test, 20.0))

# --- simulations ---------------------------------------------------------
x, y, z = np.arange(NX, dtype=float), np.arange(NY, dtype=float), np.arange(NZ, dtype=float)
fields = []
t0 = time.time()
for i in range(NREAL):
    srf = gs.SRF(model, seed=SEED + i, mode_no=MODE_NO)
    f = srf.structured((x, y, z))
    fields.append(f)
print(f"{NREAL} realizations {NX}x{NY}x{NZ} in {time.time() - t0:.1f} s")
fields = np.array(fields)
print("global mean: %.3f | mean variance: %.3f" % (fields.mean(), fields.var(axis=(1, 2, 3)).mean()))


# --- exact variograms by grid shift (all pairs) -----------------------------
def vario_shift(f, step):
    """Variogram along the integer direction 'step'=(sx,sy,sz); lag k -> k*|step|."""
    sx, sy, sz = step
    nx, ny, nz = f.shape
    lags, gam = [], []
    k = 1
    while True:
        dx, dy, dz = k * sx, k * sy, k * sz
        if abs(dx) >= nx or abs(dy) >= ny or abs(dz) >= nz:
            break
        sl_a = tuple(slice(max(0, -d), n - max(0, d)) for d, n in zip((dx, dy, dz), (nx, ny, nz)))
        sl_b = tuple(slice(max(0, d), n - max(0, -d)) for d, n in zip((dx, dy, dz), (nx, ny, nz)))
        gam.append(0.5 * np.mean((f[sl_a] - f[sl_b]) ** 2))
        lags.append(k * np.sqrt(sx**2 + sy**2 + sz**2))
        k += 1
    return np.array(lags), np.array(gam)


dirs = {
    "N45E (major axis, range 20)": ((1, 1, 0), 0),
    "N45W (intermediate axis, range 10)": ((1, -1, 0), 1),
    "Vertical (minor axis, range 5)": ((0, 0, 1), 2),
}
MAXLAG = {0: 40, 1: 40, 2: 10}

# cross-check: gstools vario_estimate (structured grid, Cartesian axes) vs own function
bins = np.arange(0.5, 11.5, 1.0)
bc, gv = gs.vario_estimate((x, y, z), fields[0], bins, direction=((0, 0, 1),), mesh_type="structured")
lz, gz = vario_shift(fields[0], (0, 0, 1))
print("check vertical gstools vs shift (first 5 lags):", np.round(gv[:5], 4), np.round(gz[:5], 4))

fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
stats = {}
for ax, (name, (step, axis)) in zip(axs, dirs.items()):
    allg = []
    for f in fields:
        l, g = vario_shift(f, step)
        allg.append(g)
        m = l <= MAXLAG[axis]
        ax.plot(l[m], g[m], color="0.7", lw=0.8)
    allg = np.array(allg)
    mean_g = allg.mean(axis=0)
    m = l <= MAXLAG[axis]
    th = model.vario_axis(l, axis=axis)
    ax.plot(l[m], mean_g[m], "b-", lw=1.8, label="mean of %d realizations" % NREAL)
    ax.plot(l[m], th[m], "k-", lw=2.4, label="model")
    ax.set_title(name, fontsize=10)
    ax.set_xlabel("lag")
    ax.set_ylim(0, 1.5)
    ax.grid(alpha=0.3)
    # RMS deviation within the range (where the comparison is meaningful)
    rng = {0: 20, 1: 10, 2: 5}[axis]
    mm = l <= rng
    stats[name] = (np.sqrt(np.mean((mean_g[mm] - th[mm]) ** 2)), np.max(np.abs(mean_g[mm] - th[mm])))
axs[0].set_ylabel("semivariogram")
axs[0].legend()
plt.tight_layout()
plt.savefig(f"{RES}/gstools_variograms_3d.png", dpi=130)

print("\nDeviation of the mean of the realizations from the model within the range (RMS, max):")
for k, (a, b) in stats.items():
    print(f"  {k}: RMS={a:.3f}  max={b:.3f}")

# --- one slice for visual inspection ------------------------------------
fig, ax = plt.subplots(figsize=(5, 4.5))
im = ax.imshow(fields[0][:, :, 8].T, origin="lower", cmap="viridis")
ax.set_title("Realization 1, slice z=8")
plt.colorbar(im)
plt.tight_layout()
plt.savefig(f"{RES}/slice_xy.png", dpi=110)
