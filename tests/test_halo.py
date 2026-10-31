'''Analytic SHM halo integrals (shm_halo) at vmin = vesc - ve, where the closed form of eta_0
changes branch: eta_i must be continuous there and not drop to 0 because of rounding.'''
import ctypes

import numpy as np
import pytest

import rapidd

_shm_halo = rapidd.core._crapidd.shm_halo
_shm_halo.argtypes = [ctypes.c_double] * 4 + [ctypes.c_int]
_shm_halo.restype = ctypes.c_double

rng = np.random.default_rng(2026)
# (vesc, v0, ve) in km/s: paper / MadDM defaults and random values
HALOS = [(544., 220., 244.), (544., 238., 250.5), (544., 220., 232.), (600., 230., 240.)] \
    + [tuple(rng.uniform([450., 150., 200.], [650., 280., 280.])) for _ in range(500)]


def shm_halo(vmin, vesc, v0, ve, i):
    return _shm_halo(vmin, vesc, v0, ve, i)


@pytest.mark.parametrize("i", (0, 1, 2))
def test_continuous_at_vesc_minus_ve(i):
    bad = []
    for vesc, v0, ve in HALOS:
        vm = vesc - ve
        at = shm_halo(vm, vesc, v0, ve, i)
        for nb in (vm * (1 - 1e-12), vm * (1 + 1e-12)):
            if not np.isclose(at, shm_halo(nb, vesc, v0, ve, i), rtol=1e-8, atol=0):
                bad.append((vesc, v0, ve, at))
                break
    assert not bad, bad[:5]


def test_eta0_against_numerical_integral():
    '''eta_0 = int_{vmin} f(v)/v dv of the boosted, truncated Maxwellian (paper halo), also at
    the branch point vmin = vesc - ve = 300 km/s'''
    vesc, v0, ve = 544., 220., 244.
    from scipy.special import erf
    nesc = erf(vesc / v0) - 2 * vesc / v0 * np.exp(-(vesc / v0)**2) / np.sqrt(np.pi)
    v = np.linspace(1e-3, vesc + ve, 200001)
    # angular integral of exp(-|v + ve|^2/v0^2) over cos < cmax (|v + ve| < vesc), done exactly
    cmax = np.clip((vesc**2 - v**2 - ve**2) / (2 * v * ve), -1., 1.)
    a = 2 * v * ve / v0**2
    ang = np.exp(-(v**2 + ve**2) / v0**2) * (np.exp(a) - np.exp(-a * cmax)) / a
    fv = 2 * np.pi * v**2 * ang / (nesc * np.pi**1.5 * v0**3)
    for vmin in (100., 299., 300., 301., 600., 750.):
        sel = v >= vmin
        g = fv[sel] / v[sel]
        ref = np.sum((g[1:] + g[:-1]) / 2 * np.diff(v[sel]))
        np.testing.assert_allclose(shm_halo(vmin, vesc, v0, ve, 0), ref, rtol=1e-3)
