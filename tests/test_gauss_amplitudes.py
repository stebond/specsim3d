"""Check of the Gaussian-coefficients option (igauss = 1 / --gauss).

16 x 16 x 16 compact periodic grid (exponential model, practical range 5, 300 realizations):
  * fixed amplitudes (default): the sample variance of every realization is the same (Parseval);
  * Gaussian coefficients: the sample variance fluctuates, with the same mean;
  * the mean variogram along x is the same in both modes within Monte Carlo noise;
  * the pooled values have a kurtosis excess close to zero for the Gaussian option;
  * the fourth value of the parameter line and --gauss give identical realizations.
Autore: Stefano Bonduà
"""
import os
import subprocess
import sys

import numpy as np

WD = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WD)
EXE = os.path.join(ROOT, "specsim3d.exe" if os.name == "nt" else "specsim3d")
WORK = os.path.join(WD, "work")
os.makedirs(WORK, exist_ok=True)
N, NREAL = 16, 300


def run(base, fourth="", flag=None, seed=4711):
    with open(os.path.join(WORK, base + ".par"), "w") as fh:
        fh.write("Parameters for SPECSIM3D\n************************\n\nSTART OF PARAMETERS:\n"
                 f"0 0 0\n1 1 1\n{N} {N} {N}\n0.0 1\n2 1.0 5 5 5 0 0 0\n{seed} {NREAL} 0 {fourth}\n{base}\n1\n3\n0\n")
    cmd = [EXE, base + ".par"] + ([flag] if flag else [])
    r = subprocess.run(cmd, cwd=WORK, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stdout + r.stderr)
    f = np.fromfile(os.path.join(WORK, base + ".bin"), dtype=np.float32).reshape(NREAL, N, N, N).astype(np.float64)
    return f, r.stdout


def vario_x(f, kmax=5):
    return np.array([0.5 * np.mean((f[:, :, :, k:] - f[:, :, :, :-k]) ** 2) for k in range(1, kmax + 1)])


def main():
    if not os.path.exists(EXE):
        sys.exit("Executable not found: " + EXE)
    fa, _ = run("t_gauss_a")
    fb, out = run("t_gauss_b", fourth="1")
    fc, _ = run("t_gauss_c", flag="--gauss")
    va = fa.reshape(NREAL, -1).var(axis=1)
    vb = fb.reshape(NREAL, -1).var(axis=1)
    ok = True

    def check(name, cond, txt):
        nonlocal ok
        ok &= bool(cond)
        print(f"{'PASS' if cond else 'FAIL'}  {name}: {txt}")

    check("fixed amplitudes", va.std() < 1e-5, f"sd of the sample variance {va.std():.2e} (deterministic)")
    check("Gaussian coefficients", vb.std() > 0.03, f"sd of the sample variance {vb.std():.3f}")
    se = vb.std() / np.sqrt(NREAL)
    check("same mean variance", abs(vb.mean() - va.mean()) < 4 * se,
          f"{va.mean():.4f} (fixed) vs {vb.mean():.4f} +/- {se:.4f} (Gaussian)")
    d = np.abs(vario_x(fa) - vario_x(fb)).max()
    check("same mean variogram", d < 0.04, f"max difference along x, lags 1-5: {d:.4f}")
    z = fb.ravel() / fb.std()
    kur = np.mean(z ** 4) - 3.0
    check("Gaussian kurtosis", abs(kur) < 0.1, f"pooled kurtosis excess {kur:.3f}")
    check("parameter file value = --gauss", np.array_equal(fb, fc), "identical realizations")
    check("message", "Gaussian" in out, "program reports the option")
    print("\nALL TESTS PASSED" if ok else "\nSOME TESTS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
