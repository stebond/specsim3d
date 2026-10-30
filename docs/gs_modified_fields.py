#!/usr/bin/env python3
"""Helper of make_figures.py: GSTools fields with the author-modified spherical 3D sampling.
Runs in its own process so that the patch never affects the other (released-version) GSTools results.
Usage: gs_modified_fields.py <model_id> <nx> <ny> <nz> <nreal> <seed> <out.npy>"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "gstools_patch"))
sys.path.insert(0, HERE)
import gstools_spherical3d_patch  # noqa: F401  (applies the patch on import)
import make_figures as mf
it, nx, ny, nz, nreal, seed = (int(v) for v in sys.argv[1:7])
np.save(sys.argv[7], mf.gs_fields(it, (nx, ny, nz), nreal, seed))
