import time
import numpy as np
import gstools as gs

x, y, z = np.arange(64.0), np.arange(64.0), np.arange(16.0)


def vario_shift(f, step, kmax):
    sx, sy, sz = step
    out = []
    for k in range(1, kmax + 1):
        d = (k * sx, k * sy, k * sz)
        sa = tuple(slice(max(0, -di), n - max(0, di)) for di, n in zip(d, f.shape))
        sb = tuple(slice(max(0, di), n - max(0, -di)) for di, n in zip(d, f.shape))
        out.append(0.5 * np.mean((f[sa] - f[sb]) ** 2))
    return np.array(out)


# 1) check of the estimator: gstools with the smallest angular tolerance
srf = gs.SRF(gs.Spherical(dim=3, var=1, len_scale=[20, 10, 5], angles=[np.pi / 4, 0, 0]), seed=1, mode_no=1024)
f = srf.structured((x, y, z))
bins = np.arange(0.5, 10.5, 1.0)
bc, gv = gs.vario_estimate((x, y, z), f, bins, direction=((0, 0, 1),), angles_tol=1e-3, mesh_type="structured")
print("check vertical (angles_tol=1e-3) gstools:", np.round(gv[:6], 4))
print("                              own shift  :", np.round(vario_shift(f, (0, 0, 1), 6), 4))

# 2) bias at short lag / sill as a function of model and number of modes (4 realizations each)
cases = [
    ("Spherical modes=1024", lambda: gs.Spherical, 1024),
    ("Spherical modes=8192", lambda: gs.Spherical, 8192),
    ("Exponent modes=1024", lambda: gs.Exponential, 1024),
    ("Gaussian  modes=1024", lambda: gs.Gaussian, 1024),
]
for label, getm, nm in cases:
    M = getm()
    model = M(dim=3, var=1, len_scale=[20, 10, 5], angles=[np.pi / 4, 0, 0])
    gz1, gd1, gsill, var = [], [], [], []
    t0 = time.time()
    for i in range(4):
        fl = gs.SRF(model, seed=100 + i, mode_no=nm).structured((x, y, z))
        gz1.append(vario_shift(fl, (0, 0, 1), 1)[0])
        g = vario_shift(fl, (1, 1, 0), 30)
        gd1.append(g[0])
        gsill.append(g[-10:].mean())
        var.append(fl.var())
    th_z = model.vario_axis(np.array([1.0]), axis=2)[0]
    th_d = model.vario_axis(np.array([np.sqrt(2)]), axis=0)[0]
    print(f"{label}: gamma_vert(1)={np.mean(gz1):.3f} (model {th_z:.3f}) | "
          f"gamma_N45E(1.41)={np.mean(gd1):.3f} (model {th_d:.3f}) | "
          f"plateau N45E lags 21-30={np.mean(gsill):.3f} | var={np.mean(var):.3f} | {time.time() - t0:.0f}s")
