'''DM-electron scattering in silicon (crystal response W_1 from the QCDark2 dielectric function).

Reference: tests/data/si_qcdark2_reference.txt, dR/dE of the QCDark2 rate formula (arXiv:2603.12326, Eqs. 6 and 10)
computed from the original Si_comp.h5 by lib/crystal_responses/generate/qcdark2_to_rapidd.py.
'''
import os
import warnings

import numpy as np
import pytest
from scipy.integrate import simpson
from scipy.special import erf

import rapidd
from rapidd.electron import (coeff_from_sigma_e, difrate_electron, difrate_dQ, total_rate_electron,
                             read_ionization_table, UnsupportedCrystalCouplingError, CrystalTruncationWarning,
                             crystal_edge_fractions, _crystal_warned)

from conftest import DATA_DIR

REF = os.path.join(DATA_DIR, "si_qcdark2_reference.txt")
ALPHA, ME_EV, C_KMS = 1. / 137.03599908, 5.1099894e5, 299792.458     # as in QCDark2
CM_INV_EV, KG_EV, YEAR_S = 50677.3093773, 5.609588e35, 3.15576e7
QC2 = {'v0': 238., 'vEarth': 250.2, 'vEscape': 544., 'rhoX': 0.3, 'sigma_e': 1e-37, 'M_cell': 5.2164e10}
Q_EV = (0.01 + 0.02 * np.arange(rapidd.CRYSTAL_NQ)) * ALPHA * ME_EV
E_EV = 0.1 * np.arange(rapidd.CRYSTAL_NE)
DAMICM_EFF = {2: 0.38, 3: 0.65, 4: 0.79}


def load_reference():
    names = open(REF).readlines()[2].split()[1:]
    return names, np.loadtxt(REF)


def load_tables():
    d = rapidd.crystal_resp_dir
    return tuple(np.loadtxt(os.path.join(d, f), comments="//").T                 # (q, E)
                 for f in ("Si_W1.txt", "Si_W1_screened.txt", "Si_W1_nolfe.txt"))


def eta_mb(q, E, m):
    '''QCDark2 get_eta_MB (inverse mean speed, units of 1/c), reimplemented'''
    v0, ve, vesc = QC2['v0'] / C_KMS, QC2['vEarth'] / C_KMS, QC2['vEscape'] / C_KMS
    vmin = q[:, None] / (2 * m) + E[None, :] / q[:, None]
    a = -4 * ve * np.exp(-(vesc / v0)**2) + np.sqrt(np.pi) * v0 * (erf((vmin + ve) / v0) - erf((vmin - ve) / v0))
    b = -2 * (ve + vesc - vmin) * np.exp(-(vesc / v0)**2) + np.sqrt(np.pi) * v0 * (erf(vesc / v0) - erf((vmin - ve) / v0))
    val = np.where(vmin < vesc - ve, a, np.where(vmin < vesc + ve, b, 0.))
    K = v0**3 * (-2 * np.pi * (vesc / v0) * np.exp(-(vesc / v0)**2) + np.pi**1.5 * erf(vesc / v0))
    return v0**2 * np.pi / (2 * ve * K) * val


def rate_from_w1(W, m_ev, light, n_sub=1, M_cell=QC2['M_cell']):
    '''
    dR/dE [events/kg/year/eV] of the crystal rate n/(128 pi m^2 m_e^2 M_cell) int dq q eta R_1 W_1 / Delta E with
    R_1 = 16 pi m^2 m_e^2 sigma_e/mu^2 F_DM^2, written independently of the C library: Simpson's rule on the q grid
    (n_sub = 1, the QCDark2 quadrature) or on n_sub sub-intervals per cell with W_1 interpolated log-log in q
    (linearly where it vanishes, the RAPIDD quadrature).
    '''
    q = np.linspace(Q_EV[0], Q_EV[-1], (Q_EV.size - 1) * n_sub + 1)
    Wf = np.zeros((q.size, W.shape[1]))
    Wf[::n_sub] = W
    W0, W1 = W[:-1], W[1:]
    for j in range(1, n_sub):
        t = j / n_sub
        qq = Q_EV[:-1] + t * (Q_EV[1:] - Q_EV[:-1])
        pw = (np.log(qq / Q_EV[:-1]) / np.log(Q_EV[1:] / Q_EV[:-1]))[:, None]
        with np.errstate(all='ignore'):
            Wf[j::n_sub] = np.where((W0 > 0) & (W1 > 0), W0 * (W1 / W0)**pw, W0 + t * (W1 - W0))
    mu = m_ev * ME_EV / (m_ev + ME_EV)
    F2 = (ALPHA * ME_EV / q)**4 if light else np.ones_like(q)
    rho = QC2['rhoX'] * 1e9 / CM_INV_EV**3
    sig = QC2['sigma_e'] * CM_INV_EV**2
    I = simpson(q[:, None] * F2[:, None] * eta_mb(q, E_EV, m_ev) * Wf, x=q, axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.where(E_EV > 0, (rho / m_ev) * sig / (8 * mu**2) * I / M_cell / E_EV, 0.)
    return r * CM_INV_EV * C_KMS * 1e5 * YEAR_S * KG_EV


def recoil_spectrum(E, dRdE, n_max=10):
    '''QCDark2 recoil_spectrum (rate per number of e-h pairs), reimplemented: rectangle rule, dE = 0.1 eV'''
    Et, p = read_ionization_table()
    sel = (E >= Et.min()) & (E <= Et.max())
    return np.array([np.sum(dRdE[sel] * 0.1 * np.interp(E[sel], Et, p[n])) for n in range(n_max)])


@pytest.fixture(scope="module")
def qcdark2_halo(tmp_path_factory):
    '''Halo of QCDark2 (v0 = 238, v_E = 250.2, vesc = 544 km/s) and the silicon tables'''
    path = str(tmp_path_factory.mktemp("halo") / "halo_qcdark2.dat")
    rapidd.calc_halo(path, vesc=QC2['vEscape'], v0=QC2['v0'], ve=QC2['vEarth'], i=2)
    rapidd.read_halo(path)
    rapidd.read_crystal_responses()
    rapidd.reset_coefficients()
    yield
    rapidd.reset_coefficients()
    rapidd.read_halo(rapidd.halo_path)


def rapidd_rate(m_ev, light, screening, M_cell=QC2['M_cell'], op=1):
    '''RAPIDD Si dR/dE in events/kg/year/eV, rescaled to the QCDark2 M_cell (c_op from sigma_e, c_1 by default)'''
    m = m_ev * 1e-9
    rapidd.reset_coefficients()
    rapidd.set_any_Ecoeff(coeff_from_sigma_e(QC2['sigma_e'], m), op, interaction="long-range" if light else "contact")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CrystalTruncationWarning)
        r = difrate_electron(QC2['rhoX'], m, jchi=0.5, target="Si", screening=screening)
    rapidd.reset_coefficients()
    M_rapidd = rapidd.CRYSTAL_CELL_MASS_AU["Si"] * rapidd.amu * 1e9
    return r['E'] * 1e3, r['dRdE'] * 1e-3 * 365.25 * M_rapidd / M_cell


@pytest.mark.parametrize("table", ["none", "RPA", "noLFE"])
def test_tables_reproduce_qcdark2(table):
    '''
    W_1 tables in the RAPIDD rate formula with the QCDark2 quadrature = QCDark2 formula on the original eps
    (Si_comp.h5 unscreened and RPA screened, Si_nolfe.h5 unscreened)
    '''
    names, ref = load_reference()
    W = dict(zip(("none", "RPA", "noLFE"), load_tables()))[table]
    for col, name in enumerate(names[1:], start=1):
        m, med, scr = name.split("_")
        if scr != table:
            continue
        r = rate_from_w1(W, float(m[1:]), med == "light")
        R = ref[:, col]
        sel = R > 1e-8 * R.max()
        assert np.allclose(r[sel], R[sel], rtol=1e-5), name


@pytest.mark.parametrize("m_ev", [1e6, 1e7, 1e9])
@pytest.mark.parametrize("light", [False, True])
def test_kernel_vs_independent_quadrature(qcdark2_halo, m_ev, light):
    '''
    C kernel = independent reimplementation of the same quadrature (4 log-log sub-intervals per q cell): c_1 with
    the screened and unscreened tables, c_4 (R_1 = jchi(jchi+1)/4 c_4^2 = 3/16 of c_1^2 for the same coefficient)
    with the table without local field effects, for both values of screening
    '''
    W1, W1s, W1n = load_tables()
    for screening, W, op, f in ((True, W1s, 1, 1.), (False, W1, 1, 1.), (True, W1n, 4, 3. / 16.),
                                (False, W1n, 4, 3. / 16.)):
        E, r = rapidd_rate(m_ev, light, screening, op=op)
        ref = f * rate_from_w1(W, m_ev, light, n_sub=4)[1:]
        sel = ref > 1e-4 * ref.max()
        assert np.allclose(r[sel], ref[sel], rtol=2e-3), (m_ev, light, screening, op)


@pytest.mark.parametrize("m_ev", [1e6, 3e6, 1e7, 1e8, 1e9])
@pytest.mark.parametrize("light", [False, True])
@pytest.mark.parametrize("screening", [True, False])
def test_rate_vs_qcdark2(qcdark2_halo, m_ev, light, screening):
    '''
    Total rate and rate per e-h pair vs the QCDark2 formula. The two differ only by the interpolation of W_1 between
    the q points (QCDark2: Simpson on the grid), up to ~20% in single bins near the gap; with 32 sub-intervals the
    RAPIDD quadrature changes by < 1e-3.
    '''
    names, ref = load_reference()
    R = ref[:, names.index("m%.0e_%s_%s" % (m_ev, "light" if light else "heavy", "RPA" if screening else "none"))]
    E, r = rapidd_rate(m_ev, light, screening)
    tol_tot = 0.03 if m_ev < 2e6 else 0.005
    assert np.sum(r) / np.sum(R[1:]) == pytest.approx(1., abs=tol_tot)
    RQ, RQ_ref = recoil_spectrum(E, r), recoil_spectrum(E_EV, R)
    for Q in range(1, 6):
        if RQ_ref[Q - 1] > 1e-3 * RQ_ref.max():
            assert RQ[Q - 1] / RQ_ref[Q - 1] == pytest.approx(1., abs=0.05), Q


def test_dRdQ_normalization(qcdark2_halo):
    '''sum_Q dR/dQ (Q <= 20) = rate with 1.1 eV <= Delta E <= 50 eV (normalized probabilities)'''
    rapidd.set_any_Ecoeff(1., 1)
    r = difrate_electron(0.3, 0.01, target="Si")
    dq = difrate_dQ(r, Q_max=20)
    E = np.round(r['E'] * 1e3, 5)
    sel = E >= 1.1
    assert np.sum(dq['dRdQ']) == pytest.approx(np.sum(r['dRdE'][sel]) * 1e-4, rel=1e-6)
    assert list(dq['Q']) == list(range(1, 21))
    assert total_rate_electron(None, 0.01, rate=r) == pytest.approx(np.sum(r['dRdE']) * 1e-4, rel=2e-3)


def test_screening_per_coupling(qcdark2_halo):
    '''
    RPA screening acts on the density couplings (c_1, c_11), not on the spin couplings (c_4, c_6, c_9, c_10), as in
    arXiv:2402.06817 for non-spin-polarised materials
    '''
    m = 0.05

    def rate(coeffs, screening=True):
        rapidd.reset_coefficients()
        for (op, inter), c in coeffs.items():
            rapidd.set_any_Ecoeff(c, op, interaction=inter)
        return difrate_electron(0.3, m, target="Si", screening=screening)['dRdE']

    for op in (4, 6, 9, 10):
        assert np.allclose(rate({(op, "contact"): 1.}), rate({(op, "contact"): 1.}, False), rtol=1e-12)
    s, u = rate({(1, "contact"): 1.}), rate({(1, "contact"): 1.}, False)
    assert np.sum(s) < 0.9 * np.sum(u)
    # contact c_11 (R_1 ~ q^2) comes mostly from large q, where the screening is weak
    s, u = rate({(11, "contact"): 1.}), rate({(11, "contact"): 1.}, False)
    assert np.sum(s) < np.sum(u) and not np.allclose(s, u, rtol=1e-3)
    # no density-spin interference: c_1 + c_4 = c_1 alone + c_4 alone, also with long-range couplings
    both = rate({(1, "contact"): 1., (4, "contact"): 2., (6, "long-range"): 0.5})
    sep = rate({(1, "contact"): 1.}) + rate({(4, "contact"): 2., (6, "long-range"): 0.5})
    assert np.allclose(both, sep, rtol=1e-10)
    # spin couplings: independent-particle response (no local field effects), different from the unscreened
    # density response with local field effects (arXiv:2402.06817; LFE reduce it by ~25-30%)
    r4, r1 = rate({(4, "contact"): 1.}, False), rate({(1, "contact"): 1.}, False)
    assert 1.2 < np.sum(r4) / (0.75 / 4. * np.sum(r1)) < 1.6
    rapidd.reset_coefficients()


@pytest.mark.parametrize("op", [3, 5, 7, 8, 12, 13, 14, 15])
@pytest.mark.parametrize("interaction", ["contact", "long-range"])
def test_unsupported_operators(qcdark2_halo, op, interaction):
    rapidd.reset_coefficients()
    rapidd.set_any_Ecoeff(1., 1)
    rapidd.set_any_Ecoeff(1., op, interaction=interaction)
    with pytest.raises(UnsupportedCrystalCouplingError, match="c_%d " % op):
        difrate_electron(0.3, 0.1, target="Si")
    rapidd.reset_coefficients()


def test_crystal_arguments(qcdark2_halo):
    rapidd.set_any_Ecoeff(1., 1)
    with pytest.raises(ValueError):
        difrate_electron(0.3, 0.1, target="Si", shells=["3p"])
    with pytest.raises(ValueError):
        difrate_electron(0.3, 0.1, target="Si", resp=2)
    with pytest.raises(ValueError):
        difrate_electron(0.3, 0.1, target="Ge")
    rapidd.reset_coefficients()


def test_truncation_warning(qcdark2_halo):
    _crystal_warned.clear()
    rapidd.set_any_Ecoeff(1., 6)
    with pytest.warns(CrystalTruncationWarning, match="q in"):
        r = difrate_electron(0.3, 1.0, target="Si")
    assert crystal_edge_fractions(r)[0] > 0.1
    rapidd.reset_coefficients()
    rapidd.set_any_Ecoeff(1., 1)
    with warnings.catch_warnings():
        warnings.simplefilter("error", CrystalTruncationWarning)
        difrate_electron(0.3, 0.01, target="Si")
    rapidd.reset_coefficients()


def test_damicm_limit_shape(qcdark2_halo):
    '''
    Sanity check only: the shape of 1/S(m), S = signal with 2-4 e- weighted by the DAMIC-M pattern efficiencies,
    follows the published DAMIC-M QCDark heavy-mediator limit within a factor 3 between 2 MeV and 1 GeV (QCDark1 with
    analytic screening there, no diffusion or pattern details here; the largest deviations are near threshold).
    '''
    lim = np.loadtxt(os.path.join(rapidd.core.experiments_dir, "DAMIC-M", "DAMIC-M_2025_DMe_heavy_QCDark.dat"))
    masses = np.array([m for m in lim[:, 0] if 2e-3 <= m <= 1.])
    shape = []
    for m in masses:
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(coeff_from_sigma_e(1e-37, m), 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", CrystalTruncationWarning)
            dq = difrate_dQ(difrate_electron(0.3, m, target="Si"), Q_max=4)
        S = sum(eff * dq['dRdQ'][Q - 1] for Q, eff in DAMICM_EFF.items())
        shape.append(1e-37 / S)
    rapidd.reset_coefficients()
    ratio = np.interp(masses, lim[:, 0], lim[:, 1]) / np.array(shape)
    ratio /= np.exp(np.mean(np.log(ratio)))
    assert np.all((ratio > 1. / 3.) & (ratio < 3.)), ratio
