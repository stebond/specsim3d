"""Runtime patch for GSTools 1.7.0: exact radial spectral cdf/ppf for the 3D Spherical model.

Usage:  import gstools_spherical3d_patch   # before creating Spherical(dim=3) models / SRF objects

Fixes the spurious nugget of RandMeth for gs.Spherical(dim=3): the numerical hankel transform gives small
negative spectral densities at large k, spectral_rad_pdf took abs() of them and the MCMC sampler (emcee)
drifted to huge wave numbers. This patch (1) gives Spherical(dim=3) an analytic radial cdf and ppf, so
sampling="auto" uses inversion and no MCMC, and (2) clips instead of abs() in spectral_rad_pdf, which also
improves other models that still use MCMC (e.g. Spherical 2D, Stable).
Author: Stefano Bonduà (analysis and patch prepared with Claude)
"""
import numpy as np
from scipy import special as sps
import gstools
from gstools.covmodel import base as _base, tools as _tools


def _sph3d_rad_cdf(x):
    """(6/pi) * integral_0^x j1(t)^2 dt  with j1 the spherical Bessel function of order 1."""
    x = np.asarray(x, dtype=np.double)
    res = np.empty_like(x)
    small = x < 1.0
    xs = x[small]
    res[small] = (xs**3 / 27.0 - xs**5 / 225.0 + xs**7 / 3675.0 - 4.0 * xs**9 / 382725.0
                  + 2.0 * xs**11 / 7203735.0 - xs**13 / 184459275.0 + xs**15 / 12314176875.0)
    xl = x[~small]
    x2 = 2.0 * xl
    res[~small] = (sps.sici(x2)[0] / 3.0 + np.cos(x2) / (6.0 * xl) + np.sin(x2) / (3.0 * xl**2)
                   + np.cos(x2) / (6.0 * xl**3) - 0.5 / xl - 1.0 / (6.0 * xl**3))
    return 6.0 / np.pi * res


_X = np.concatenate(([0.0], np.geomspace(1e-3, 1e6, 40000)))
_F = _sph3d_rad_cdf(_X)


def _spectral_rad_cdf(self, r):
    if self.dim != 3:
        return None
    return _sph3d_rad_cdf(np.abs(np.asarray(r, dtype=np.double)) * self.len_rescaled / 2.0)


def _spectral_rad_ppf(self, u):
    if self.dim != 3:
        return None
    u = np.asarray(u, dtype=np.double)
    x = np.interp(u, _F, _X)
    x = np.where(u > _F[-1], 3.0 / (np.pi * np.maximum(1.0 - u, 1e-300)), x)
    return 2.0 * x / self.len_rescaled


def _spectral_rad_pdf(model, r):
    r = np.asarray(np.abs(r), dtype=np.double)
    if model.dim > 1:
        r_gz = np.logical_not(np.isclose(r, 0))
        res = np.zeros_like(r, dtype=np.double)
        res[r_gz] = _tools.rad_fac(model.dim, r[r_gz]) * model.spectral_density(r[r_gz])
    else:
        res = _tools.rad_fac(model.dim, r) * model.spectral_density(r)
    res[np.logical_not(np.isfinite(res))] = 0.0
    return np.maximum(res, 0.0)


gstools.Spherical.spectral_rad_cdf = _spectral_rad_cdf
gstools.Spherical.spectral_rad_ppf = _spectral_rad_ppf
gstools.Spherical._has_cdf = lambda self: self.dim == 3
gstools.Spherical._has_ppf = lambda self: self.dim == 3
_tools.spectral_rad_pdf = _spectral_rad_pdf
_base.spectral_rad_pdf = _spectral_rad_pdf
