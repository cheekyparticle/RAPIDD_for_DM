'''DM-electron ionization rate for the NREFT operators: independent evaluation of the velocity
integral, R_j W_j decomposition, input checks, model helpers and limits.'''
import os

import numpy as np
import pytest
from scipy.special import erf

import rapidd

from test_electron_responses import R_paper, kinematics, couplings

M_E = rapidd.M_ELECTRON
HBARC = 1.973269804e-14       # GeV cm
HBAR = 6.582119569e-25        # GeV s
C_KMS = 299792.458
M_XE_KG = 131.293 * 1.66053906660e-27
GEV_TO_KG = 1.78266192e-27
V0, VESC, VE = 220., 544., 244.   # halo of arXiv:1912.08204
XE_5P = 9                         # index of 5p in rapidd.XE_SHELLS
EB_5P = 12.44330433672413e-9      # GeV


def shm_speed_distribution(v):
    '''Lab-frame speed distribution f(v) (km/s)^-1 of the SHM of App. E of arXiv:1912.08204,
    with the angular integral done numerically (independent of the RAPIDD halo code)'''
    nesc = erf(VESC / V0) - 2 * VESC / V0 * np.exp(-(VESC / V0)**2) / np.sqrt(np.pi)
    cos = np.linspace(-1., 1., 4001)
    w2 = v[:, None]**2 + VE**2 + 2 * v[:, None] * VE * cos[None, :]
    f = np.exp(-w2 / V0**2) * (w2 < VESC**2)
    ang = np.sum((f[:, 1:] + f[:, :-1]) / 2 * np.diff(cos), axis=1)
    return 2 * np.pi * v**2 * ang / (nesc * np.pi**1.5 * V0**3)


def independent_dRdlnE(coeffs, mchi, rho, ik_list, drop_vperp_const=False, resp=(1, 2, 3, 4)):
    '''dR/dlnE [events/kg/day] of the shell 5p from Eqs. (ionization spectrum) and (RW):
    n_chi/(128 pi mchi^2 m_e^2) int dq q int d^3v f(v)/v Theta(v - vmin) sum_j R_j(v) W_j,
    with R_j evaluated at each v with the full (v_el^perp)^2 (or with (v_el^perp)^2 = v^2 if
    drop_vperp_const, to check that the constant part matters), summing the terms j in resp.
    Same grids and q quadrature (left sum) as RAPIDD.'''
    base = os.path.join(rapidd.core.atomic_resp_dir, "atomic_response_%i", "Xe_5p.txt")
    W = [np.loadtxt(base % j, comments="//") for j in (1, 2, 3, 4)]
    k = np.exp(np.log(1e-7) + np.log(1e3) * np.arange(101) / 99)
    E = np.sqrt(M_E**2 + k**2) - M_E
    q = np.exp(np.log(1e-6) + np.log(1e3) * np.arange(201) / 199)
    v = np.linspace(0., VESC + VE, 6001)[1:]
    fv = shm_speed_distribution(v)
    n_chi = rho / mchi * HBARC**3  # GeV^3
    out = []
    for ik in ik_list:
        dE = E[ik] + EB_5P
        integ = 0.
        for iq in range(200):
            vmin = (dE / q[iq] + q[iq] / (2 * mchi)) * C_KMS
            sel = v >= vmin
            if sel.sum() < 2:
                continue
            vn = v[sel] / C_KMS  # units of c
            c = couplings(coeffs, q[iq])
            vperp2, qv = kinematics(q[iq], dE, mchi, vn)
            if drop_vperp_const:
                vperp2 = vn**2
            R = R_paper(c, 0.5, q[iq] / M_E, vperp2, qv * np.ones_like(vn))
            RW = sum(R[j - 1] * W[j - 1][ik, iq] for j in resp)
            g = fv[sel] * C_KMS / vn * RW  # f(v) dv/v: f in (km/s)^-1, v in units of c
            vel_int = np.sum((g[1:] + g[:-1]) / 2 * np.diff(vn))
            integ += (q[iq + 1] - q[iq]) * q[iq] * vel_int
        rate = n_chi / (128 * np.pi * mchi**2 * M_E**2) * integ  # GeV, per atom
        out.append(rate / HBAR * 86400 / M_XE_KG)
    return np.array(out)


@pytest.mark.parametrize("coeffs,mchi", [({(7, "contact"): 1e-3}, 0.01),
                                         ({(7, "contact"): 1e-3}, 0.1),
                                         ({(13, "long-range"): 1e-3}, 0.1),
                                         ({(3, "contact"): 1e-3, (5, "contact"): 2e-3}, 1.0)])
def test_velocity_integral(paper_halo, coeffs, mchi):
    '''Both the v^2 and the constant part of (v_el^perp)^2 are in the rate (e.g. O_7, 10 MeV)'''
    for (op, inter), val in coeffs.items():
        rapidd.set_any_Ecoeff(val, op, interaction=inter)
    r = rapidd.difrate_electron(0.4, mchi, 0.5, shells=["5p"])
    ik = [i for i in (5, 25, 45, 60) if r['dRdlnE'][i] > 1e-3 * r['dRdlnE'].max()]
    assert len(ik) >= 2
    ref = independent_dRdlnE(coeffs, mchi, 0.4, ik)
    np.testing.assert_allclose(r['dRdlnE'][ik], ref, rtol=5e-3)
    if mchi == 0.01:
        # R_1 W_1 alone; dropping q^2/(4 mu^2)(mchi - m_e)/(mchi + m_e) - dE/mu (i.e. keeping only
        # the v^2 integrated with halo(vmin, 2)) changes it completely for light DM
        r1 = rapidd.difrate_electron(0.4, mchi, 0.5, shells=["5p"], resp=1)['dRdlnE'][ik]
        ref1 = independent_dRdlnE(coeffs, mchi, 0.4, ik, resp=(1,))
        np.testing.assert_allclose(r1, ref1, rtol=5e-3)
        wrong = independent_dRdlnE(coeffs, mchi, 0.4, ik, drop_vperp_const=True, resp=(1,))
        assert np.all(np.abs(wrong / ref1 - 1) > 0.5)


def test_response_decomposition(paper_halo):
    '''sum_j (R_j W_j only) = total; R_4 W_4 of O_15 is negative'''
    rapidd.set_any_Ecoeff(1e-3, 15, interaction="long-range")
    tot = rapidd.difrate_electron(0.4, 0.1)
    parts = [rapidd.difrate_electron(0.4, 0.1, resp=j)['dRdE'] for j in (1, 2, 3, 4)]
    np.testing.assert_allclose(sum(parts), tot['dRdE'], rtol=1e-12, atol=1e-30)
    assert np.all(parts[1] == 0.)              # O_15 has no R_2
    assert np.all(parts[3] <= 0.) and parts[3].min() < 0.
    np.testing.assert_allclose(rapidd.difrate_electron(0.4, 0.1, resp=(1, 3))['dRdE'], parts[0] + parts[2], rtol=1e-12)


def test_shell_selection(paper_halo):
    rapidd.set_any_Ecoeff(1e-3, 8)
    r_all = rapidd.difrate_electron(0.4, 0.1)
    r_5 = rapidd.difrate_electron(0.4, 0.1, shells=rapidd.XE_SHELLS_PAPER)
    np.testing.assert_allclose(r_5['dRdE'], sum(r_all['dRdE_shell'][s] for s in rapidd.XE_SHELLS_PAPER), rtol=1e-12)
    assert set(r_5['dRdE_shell']) == set(rapidd.XE_SHELLS_PAPER)
    np.testing.assert_allclose(r_5['dRdlnE'], r_5['E'] * r_5['dRdE'])
    tot = rapidd.total_rate_electron(0.4, 0.1, shells=rapidd.XE_SHELLS_PAPER)
    np.testing.assert_allclose(tot, np.trapezoid(r_5['dRdlnE'], np.log(r_5['E'])) if hasattr(np, "trapezoid")
                               else np.trapz(r_5['dRdlnE'], np.log(r_5['E'])), rtol=1e-12)
    with pytest.raises(ValueError):
        rapidd.difrate_electron(0.4, 0.1, shells=["6s"])
    with pytest.raises(ValueError):
        rapidd.difrate_electron(0.4, 0.1, resp=5)


def test_missing_atomic_response(paper_halo):
    '''W_2 is needed by O_7 but not by O_1'''
    rapidd.core._read_atomic_responses(b"/nonexistent", 2)  # marks W_2 as not read
    try:
        rapidd.set_any_Ecoeff(1e-3, 1)
        rapidd.difrate_electron(0.4, 0.1)
        rapidd.set_any_Ecoeff(1e-3, 7)
        with pytest.raises(RuntimeError, match="W_j"):
            rapidd.difrate_electron(0.4, 0.1)
        rapidd.difrate_electron(0.4, 0.1, resp=(1, 3, 4))  # W_2 not used
    finally:
        rapidd.read_atomic_responses()


def test_missing_halo_column(paper_halo, tmp_path):
    '''halo(vmin, 2) is needed by the operators depending on v_el^perp'''
    path = str(tmp_path / "halo_i0.dat")
    rapidd.calc_halo(path, profile="SHM", vesc=544., v0=220., ve=244., i=0)
    rapidd.read_halo(path)
    rapidd.set_any_Ecoeff(1e-3, 4)
    rapidd.difrate_electron(0.4, 0.1)
    rapidd.set_any_Ecoeff(1e-3, 8)
    with pytest.raises(RuntimeError, match="halo"):
        rapidd.difrate_electron(0.4, 0.1)
    rapidd.difrate_electron(0.4, 0.1, resp=(2, 3))  # the v^2 part of R_1 is not needed


def test_reset_and_coefficients():
    rapidd.set_any_Ecoeff(2., 5)
    rapidd.set_any_Ecoeff(3., 5, interaction="long-range")
    assert rapidd.Ce_val(5) == 2. and rapidd.Ce_LR_val(5) == 3.
    rapidd.reset_coefficients()
    assert rapidd.Ce_val(5) == 0. and rapidd.Ce_LR_val(5) == 0.
    with pytest.raises(ValueError):
        rapidd.set_any_Ecoeff(1., 5, interaction="short")


def test_reference_cross_section():
    m = 0.05
    c = rapidd.coeff_from_sigma_e(1e-38, m)
    np.testing.assert_allclose(rapidd.sigma_e_from_coeff(c, m), 1e-38, rtol=1e-12)
    # sigma_e = mu^2 c^2/(16 pi m^2 m_e^2) = c^2/(16 pi (m + m_e)^2)
    np.testing.assert_allclose(rapidd.sigma_e_from_coeff(c, m), c**2 / (16 * np.pi * (m + M_E)**2) * 3.893793721e-28,
                               rtol=1e-12)
    rapidd.reset_coefficients()
    rapidd.set_dark_photon_coeff(m, 1e-38, interaction="long-range")
    assert rapidd.Ce_val(1) == 0. and np.isclose(rapidd.Ce_LR_val(1), c)
    rapidd.reset_coefficients()


def test_model_coefficients():
    '''Eqs. (anapole/magnetic dipole effective couplings) of arXiv:1912.08204, g_e = 2'''
    e = np.sqrt(4 * np.pi * 7.2973525643e-3)
    qr2 = (7.2973525643e-3 * M_E)**2
    m, g = 0.1, 1e-3
    an = rapidd.electron_model_coeffs("anapole", m, g)
    assert np.isclose(an[(8, "contact")], 8 * e * M_E * m * g) and np.isclose(an[(9, "contact")], -8 * e * M_E * m * g)
    md = rapidd.electron_model_coeffs("magnetic_dipole", m, g)
    assert np.isclose(md[(1, "contact")], 4 * e * M_E * g) and np.isclose(md[(4, "contact")], 16 * e * m * g)
    assert np.isclose(md[(5, "long-range")], 16 * e * M_E**2 * m / qr2 * g)
    assert np.isclose(md[(6, "long-range")], -16 * e * M_E**2 * m / qr2 * g)
    ed = rapidd.electron_model_coeffs("electric_dipole", m, g)
    assert list(ed) == [(11, "long-range")] and np.isclose(ed[(11, "long-range")], 16 * e * m * M_E**2 / qr2 * g)
    with pytest.raises(ValueError):
        rapidd.electron_model_coeffs("charge_radius", m, g)
    rapidd.reset_coefficients()
    set_ = rapidd.set_any_E_EFT_coeff("magnetic_dipole", m, g, operators=(4,))
    assert list(set_) == [(4, "contact")] and rapidd.Ce_val(1) == 0. and rapidd.Ce_val(4) == md[(4, "contact")]
    rapidd.reset_coefficients()


def test_limits_scaling(paper_halo):
    '''The limit does not depend on the coupling used to compute the signal (rate ~ g^2)'''
    for exp in ("Xenon10", "Xenon1T"):
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(1.0, 8)
        l1 = rapidd.electron_limit(rapidd.electron_signal(exp, 0.4, 0.1), coupling=1.0)
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(1e-4, 8)
        l2 = rapidd.electron_limit(rapidd.electron_signal(exp, 0.4, 0.1), coupling=1e-4)
        np.testing.assert_allclose(l1, l2, rtol=1e-8)
    sig = {'signal': np.array([10., 0., 1.]), 'n_obs': np.array([0, 0, 5])}
    # strongest bin: N_90(0) = 2.3026 events against 10 signal events
    np.testing.assert_allclose(rapidd.electron_limit(sig, coupling=2.), 2. * np.sqrt(2.302585092994046 / 10.))
    assert np.isnan(rapidd.electron_limit({'signal': -1., 'n_obs': -1, 'n_bkg': -1.}))
    assert np.isinf(rapidd.electron_limit({'signal': np.zeros(7), 'n_obs': np.zeros(7, dtype=int)}))
