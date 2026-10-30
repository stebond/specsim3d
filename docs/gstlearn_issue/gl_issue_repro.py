"""Minimal reproducer: gstlearn simuSpectral rejects the 3D spherical covariance.
Tested with gstlearn 1.11.1 (pip, Linux x86_64, Python 3.13)."""
import gstlearn as gl

print("gstlearn", gl.__version__ if hasattr(gl, "__version__") else "?")
grid = gl.DbGrid.create(nx=[10, 10, 4], dx=[1., 1., 1.], x0=[0., 0., 0.])

for name in ("SPHERICAL", "EXPONENTIAL", "GAUSSIAN", "CUBIC", "MATERN", "CAUCHY", "STABLE", "LINEAR", "NUGGET"):
    ct = getattr(gl.ECov, name)
    try:
        model = gl.Model.createFromParam(ct, range=20., sill=1., param=1.5 if name in ("MATERN", "CAUCHY", "STABLE") else 1.,
                                         ranges=[20., 10., 5.], angles=[45., 0., 0.], flagRange=True)
        err = gl.simuSpectral(None, grid, model, nbsimu=1, seed=123)
    except Exception as e:
        err = f"exception: {e}"
    print(f"{name:12s} -> {'OK' if err == 0 else 'FAILED'}")
