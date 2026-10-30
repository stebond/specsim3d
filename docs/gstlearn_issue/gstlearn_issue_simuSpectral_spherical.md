# `simuSpectral` rejects the 3D Spherical covariance ("not valid for SPECTRAL") although its spectral density is known and non-negative

*(Written as a question / feature request: the refusal may be intended, in that case please document it and, if possible, consider adding the model.)*

## Summary
`gl.simuSpectral` aborts for `ECov.SPHERICAL` on a 3D grid with the message

```
The structure (Spherical) is not valid for SPECTRAL
(Space=RN NDim=3 NVar=1)
Error: .../src/Calculators/ACalculator.cpp@33: Check has failed. Calculation aborted
```

The same call works for Exponential, Gaussian, Matern, Cauchy and Nugget. The spherical model is the most used covariance in geostatistics, so this blocks the spectral method for most practical 3D cases.

## Why I think it should be possible
In 3D the spherical covariance C(h) = 1 - 1.5 h/a + 0.5 (h/a)^3 (h < a), 0 otherwise, is positive definite and its spectral density is non-negative and available in closed form:

f(k) ∝ [ j1(k a) / (k a) ]^2 , with j1 the spherical Bessel function of order 1,

so the radial distribution of the wave numbers has the cdf (6/π) ∫0^x j1(t)^2 dt, which has a closed form (sine integral, sin, cos), and can be sampled by inversion. (I used this expression to correct the 3D spherical sampling of GSTools; the corrected fields reproduce the model within the Monte Carlo noise.)

## Environment
- gstlearn 1.11.1 (`pip install gstlearn`), Linux x86_64, Python 3.13
- Same behaviour expected on Windows (not tested)

## Minimal reproducer
```python
import gstlearn as gl

print("gstlearn", gl.__version__ if hasattr(gl, "__version__") else "?")
grid = gl.DbGrid.create(nx=[10, 10, 4], dx=[1., 1., 1.], x0=[0., 0., 0.])

for name in ("SPHERICAL", "EXPONENTIAL", "GAUSSIAN", "CUBIC", "MATERN", "CAUCHY", "STABLE", "LINEAR", "NUGGET"):
    ct = getattr(gl.ECov, name)
    try:
        model = gl.Model.createFromParam(ct, range=20., sill=1.,
                                         param=1.5 if name in ("MATERN", "CAUCHY", "STABLE") else 1.,
                                         ranges=[20., 10., 5.], angles=[45., 0., 0.], flagRange=True)
        err = gl.simuSpectral(None, grid, model, nbsimu=1, seed=123)
    except Exception as e:
        err = f"exception: {e}"
    print(f"{name:12s} -> {'OK' if err == 0 else 'FAILED'}")
```

## Observed result
```
SPHERICAL    -> FAILED   (not valid for SPECTRAL)
EXPONENTIAL  -> OK
GAUSSIAN     -> OK
CUBIC        -> FAILED   (not valid for SPECTRAL)
MATERN       -> OK
CAUCHY       -> OK
STABLE       -> FAILED   (not valid for SPECTRAL)
LINEAR       -> FAILED   (not valid for SPECTRAL)
NUGGET       -> OK
```

## Expected result
For `SPHERICAL` (3D, anisotropic, ranges 20/10/5, angles 45/0/0) `simuSpectral` returns 0 and the realizations reproduce the model: unit variance and an experimental variogram equal, on average over realizations, to 1 - C(h) (e.g. at lag 1 along the N45E direction: 1 - C(√2) = 0.106 with a = 20).

## Check that the other models work (sanity of my set-up)
Same grid but 64 x 64 x 16, 10 realizations, ranges 20/10/5, angles 45/0/0, practical ranges (`flagRange=True`). Mean experimental variogram (all pairs) along N45E, lags k·√2:

| model | lag k=1 | k=8 |
|---|---|---|
| Exponential, theory | 0.191 | 0.817 |
| Exponential, `simuSpectral` | 0.188 | 0.815 |
| Gaussian, theory | 0.015 | 0.617* |
| Gaussian, `simuSpectral` | 0.016 | 0.639 |

(*theory for the Gaussian model 1 - exp(-3 (h/a)^2) with h = 8·√2 and a = 20.) So the method is accurate when it is available (the 0.02 difference at k=8 for the Gaussian model is compatible with the Monte Carlo noise of only 10 realizations); it is, however, slow with the default settings (about 145 s for 10 realizations on 65,536 nodes, 2 cores). If there is a parameter for the number of frequencies that I missed, a pointer in the documentation would help.

## Request
1. Either add the spherical covariance (and possibly Cubic) to the models accepted by `simuSpectral`, using the closed-form radial cdf above; or
2. document why it is excluded and which covariances are valid for the spectral method.

Thank you for the library.
