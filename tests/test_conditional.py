"""Regression test of the conditional simulation (--cond) of SPECSIM3D.

Checks, on a small 3D case with all the data in the neighbourhood:
 1) the data are reproduced at their nodes;
 2) the mean and the variance of the realizations equal the exact simple kriging estimate and variance
    within the Monte Carlo error;
 3) the normal score transform table matches the standard normal quantiles (scipy) and a back transformed
    lognormal case reproduces the data.
Run:  python test_conditional.py      (needs numpy, scipy and the compiled specsim3d executable)
Author: Stefano Bonduà
"""
import os, subprocess, sys
import numpy as np
from scipy.stats import norm

WD = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WD)
EXE = os.path.join(ROOT, "specsim3d.exe" if os.name == "nt" else "specsim3d")
WORK = os.path.join(WD, "work")
os.makedirs(WORK, exist_ok=True)
NX, NY, NZ, NREAL = 30, 30, 10, 400


def run(par, *flags):
    out = subprocess.run([EXE, par, "--bin", *flags], cwd=WORK, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stdout + out.stderr)
    return out.stdout


def write_data(name, pts, names):
    with open(os.path.join(WORK, name), "w") as fh:
        fh.write("test data\n%d\n" % len(names) + "\n".join(names) + "\n")
        for p in pts:
            fh.write(" ".join("%.8g" % v for v in p) + "\n")


def test_gaussian():
    rng = np.random.default_rng(3)
    nd = 40
    idx = rng.choice(NX * NY * NZ, nd, replace=False)
    i, j, k = np.unravel_index(idx, (NX, NY, NZ))
    v = rng.normal(size=nd)
    write_data("tc_data.dat", np.c_[i, j, k, v], ["x", "y", "z", "v"])
    open(os.path.join(WORK, "tc.par"), "w").write(
        f"0 0 0\n1 1 1\n{NX} {NY} {NZ}\n0.0 1\n1 1.0 12 8 6 30 0 0\n12345 {NREAL} 1\ntc\n")
    open(os.path.join(WORK, "tc.cpar"), "w").write(
        "tc_data.dat\n1 2 3 4 0\n-1e21 1e21\n0 0.0\n-5 5\n1 0\n1 0\n100 100 100\n0 0 0\n100\n")
    run("tc.par", "--cond", "tc.cpar")
    sims = np.fromfile(os.path.join(WORK, "tc.bin"), dtype="<f4").reshape(NREAL, NZ, NY, NX)
    sims = np.moveaxis(sims, 0, -1).transpose(2, 1, 0, 3)               # (x, y, z, realization)
    # exact simple kriging with the spherical model, major axis azimuth 30 degrees
    az = np.radians(30.0)
    u = np.array([np.sin(az), np.cos(az), 0.0]); w = np.array([np.cos(az), -np.sin(az), 0.0])

    def cov(d):
        h = np.sqrt((d @ u / 12.0) ** 2 + (d @ w / 8.0) ** 2 + (d[..., 2] / 6.0) ** 2)
        return np.where(h < 1, 1 - 1.5 * h + 0.5 * h ** 3, 0.0)
    xyz = np.c_[i, j, k].astype(float)
    K = cov(xyz[:, None] - xyz[None]) + 1e-10 * np.eye(nd)
    g = np.meshgrid(*(np.arange(n, dtype=float) for n in (NX, NY, NZ)), indexing="ij")
    P = np.c_[g[0].ravel(), g[1].ravel(), g[2].ravel()]
    kv = cov(xyz[None] - P[:, None])
    W = np.linalg.solve(K, kv.T).T
    est = (W @ v).reshape(NX, NY, NZ); var = (1 - (W * kv).sum(1)).reshape(NX, NY, NZ)
    err_data = abs(sims[i, j, k, :] - v[:, None]).max()
    rms_mean = np.sqrt(np.mean((sims.mean(-1) - est) ** 2)); mc_mean = np.sqrt(var.mean() / NREAL)
    ratio_var = sims.var(-1, ddof=1).mean() / var.mean()
    print(f"data reproduction max error {err_data:.2e}; RMS(mean - SK) {rms_mean:.4f} (Monte Carlo {mc_mean:.4f}); "
          f"mean variance ratio {ratio_var:.4f}")
    assert err_data < 1e-5
    assert rms_mean < 1.3 * mc_mean
    assert abs(ratio_var - 1) < 0.02


def test_nscore():
    rng = np.random.default_rng(5)
    nd = 200
    idx = rng.choice(64 * 64, nd, replace=False)
    x, y = (idx % 64).astype(float), (idx // 64).astype(float)
    v = np.exp(rng.normal(1.0, 0.8, nd))
    write_data("tn_data.dat", np.c_[x, y, v], ["x", "y", "v"])
    open(os.path.join(WORK, "tn.par"), "w").write("0 0 0\n1 1 1\n64 64 1\n0.0 1\n2 1.0 15 15 15 0 0 0\n777 10 1\ntn\n")
    open(os.path.join(WORK, "tn.cpar"), "w").write(
        "tn_data.dat\n1 2 0 3 0\n0 1e21\n1 0.0\n0.0 200.0\n1 0\n1 0\n30 30 30\n0 0 0\n16\n")
    run("tn.par", "--cond", "tn.cpar")
    table = np.loadtxt(os.path.join(WORK, "tn_nscore.trn"), skiprows=4)
    err_q = abs(table[:, 1] - norm.ppf((np.arange(nd) + 0.5) / nd)).max()
    sims = np.fromfile(os.path.join(WORK, "tn.bin"), dtype="<f4").reshape(10, 64, 64)
    err_d = abs(sims[:, idx // 64, idx % 64] - v[None]).max() / v.max()
    print(f"nscore table vs normal quantiles {err_q:.1e}; data reproduction after back transform {err_d:.1e}")
    assert err_q < 1e-6
    assert err_d < 1e-5


def test_parameter_file_format():
    """The conditioning block and the output codes inside the parameter file give exactly the same
    realizations as the legacy command line (--bin --cond condfile)."""
    test_gaussian()                                   # writes tc.par, tc.cpar and the legacy result tc.bin
    base = open(os.path.join(WORK, "tc.par")).read().splitlines()
    new = base[:-1] + ["tc_new", "3", "3", "2", "4", "1"] + open(os.path.join(WORK, "tc.cpar")).read().splitlines()
    open(os.path.join(WORK, "tc_new.par"), "w").write("\n".join(new) + "\n")
    out = subprocess.run([EXE, "tc_new.par"], cwd=WORK, capture_output=True, text=True)
    assert out.returncode == 0, out.stdout + out.stderr
    a = open(os.path.join(WORK, "tc.bin"), "rb").read()
    b = open(os.path.join(WORK, "tc_new.bin"), "rb").read()
    assert a == b, "parameter-file conditioning differs from the command-line conditioning"
    assert os.path.exists(os.path.join(WORK, f"tc_new_{NREAL}.vtk")), "output code 2 (VTK of every realization) missing"
    # code 4: one VTK XML ImageData file; the arrays must equal the .bin (read with the standard-library parser)
    raw = open(os.path.join(WORK, "tc_new.vti"), "rb").read()
    head = raw[:raw.index(b"<AppendedData")].decode()
    assert 'WholeExtent="0 %d 0 %d 0 %d"' % (NX - 1, NY - 1, NZ - 1) in head and head.count("<DataArray") == NREAL
    pos = raw.index(b"_", raw.index(b"<AppendedData")) + 1
    nb = 4 * NX * NY * NZ
    ref = np.fromfile(os.path.join(WORK, "tc_new.bin"), dtype="<f4").reshape(NREAL, NZ * NY * NX)
    for r in (0, NREAL // 2, NREAL - 1):
        off = pos + r * (nb + 8)
        assert int.from_bytes(raw[off:off + 8], "little") == nb
        assert np.array_equal(np.frombuffer(raw[off + 8:off + 8 + nb], dtype="<f4"), ref[r]), "vti array %d differs from .bin" % r
    # unconditional file with the new format and a bad output code must stop with a message
    bad = base[:-1] + ["tc_bad", "1", "9", "0"]
    open(os.path.join(WORK, "tc_bad.par"), "w").write("\n".join(bad) + "\n")
    out = subprocess.run([EXE, "tc_bad.par"], cwd=WORK, capture_output=True, text=True)
    assert out.returncode != 0 and "output code" in (out.stdout + out.stderr)
    print("parameter file format (conditioning block, output codes): identical to the legacy run")


if __name__ == "__main__":
    if not os.path.exists(EXE):
        sys.exit("Executable not found: build it first (see README.txt)")
    test_gaussian(); test_nscore(); test_parameter_file_format()
    print("conditional tests passed")
