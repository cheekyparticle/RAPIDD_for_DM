"""Chiral EFT structure factors (rapidd.chiral_eft): SI of Hoferichter, Klos, Menendez, Schwenk 2019 (HKMS), SD of
Hoferichter, Menendez, Schwenk 2020 (HMS) and Hu et al. 2022 (Hu2022), and the equivalent cross sections of the
LZ limits."""
import numpy as np
import pytest

import rapidd
from rapidd import chiral_eft as ce
from rapidd.core import vev

CARD = {'rhoDM': 0.3, 'vmp': 238., 'vescape': 544., 'vSun_r': 11.1, 'vSun_phi': 12.2, 'vSun_theta': 7.3,
        'vEarth_mod': 29.8,
        'SPu': 0.0153, 'SPd': 0.0191, 'SPs': 0.0447, 'SPg': 1 - 0.0153 - 0.0191 - 0.0447,
        'SNu': 0.0110, 'SNd': 0.0273, 'SNs': 0.0447, 'SNg': 1 - 0.0110 - 0.0273 - 0.0447,
        'VPu': 2.0, 'VPd': 1.0, 'VNu': 1.0, 'VNd': 2.0,
        'AVPu': 0.842, 'AVPd': -0.427, 'AVPs': -0.085, 'AVNu': -0.427, 'AVNd': 0.842, 'AVNs': -0.085,
        'SigPu': 0.84, 'SigPd': -0.23, 'SigPs': -0.046, 'SigNu': -0.23, 'SigNd': 0.84, 'SigNs': -0.046}


@pytest.fixture(autouse=True)
def default_halo():
    rapidd.read_halo(rapidd.halo_path)
    rapidd.reset_coefficients()
    yield
    rapidd.reset_coefficients()


def test_oscillator_length_matches_tables():
    """b of the fits (hbar omega = 45 A^-1/3 - 25 A^-2/3 MeV) = tabulated b of HMS Tables 6-7 and 1605.08043 Table 4"""
    for A, b_fm in ((19, 1.7623), (40, 1.9399), (73, 2.1076), (129, 2.2873), (131, 2.2925), (136, 2.3051)):
        assert ce.oscillator_length(A) * ce.HBARC == pytest.approx(b_fm, abs=1e-4)
    for k, (A, Z, J) in ce.sd_isotopes().items():
        assert ce._sd_entry(A, Z)["b_fm"] == pytest.approx(ce.oscillator_length(A) * ce.HBARC, abs=1e-4)


def test_si_normalization_at_q0():
    """F_+^M(0) = A, F_-^M(0) = Z - N for all the isotopes"""
    for target, (Z, As, _) in ce.TARGETS.items():
        for A in As:
            R = ce.si_responses(target, A, 0.)
            assert R["m+"] == pytest.approx(A)
            assert R["m-"] == pytest.approx(2 * Z - A)
            if "L2" in R:
                assert R["L2"]["m+"] == 0. and R["L2"]["m-"] == 0.


def test_xenon_fits_agree_with_1605_08043():
    """The ChiralEFT4DM xenon M and Phi'' fits are those of Table 4 of arXiv:1605.08043 (spot checks)"""
    R = ce._si_fits("Xe")
    assert R["m+"]["131"][1] == -131.26 and R["m-"]["129"][2] == -12.3687
    assert R["phi+"]["132"][0] == -28.7972 and R["phi-"]["136"][4] == 0.00559688
    assert R["pi"]["128"][0] == -2.42605
    assert R["m+"]["131L2"][1] == 2.17516


def test_sd_normalization_reproduces_spin_expectation_values():
    """S_N(0) = (2J+1)(J+1)/(pi J) <S_N>^2 (HMS Eq. 69 and Table 4): within the rounding of <S_N> (3 digits)"""
    for name, (A, Z, J) in ce.sd_isotopes().items():
        S_p, S_n = ce.spin_expectation(A, Z)
        for aA, S in (((1., 0.), S_p), ((0., 1.), S_n)):
            val = ce.sd_structure_function(A, Z, 0., aA=aA, two_body="none")
            ref = (2 * J + 1) * (J + 1) / (np.pi * J) * S**2
            assert val == pytest.approx(ref, rel=2. * 0.0005 / abs(S) + 1e-3), (name, aA)


def test_c_D_band_of_HMS_table_5():
    assert ce.c_D_from_reduction(0.09, -0.30) == pytest.approx(-6.08, abs=0.01)
    assert ce.c_D_from_reduction(0.11, -0.20) == pytest.approx(0.30, abs=0.01)


def test_fermi_integrals_against_analytic_expressions():
    """I_1^sigma, I_2^sigma by quadrature vs the analytic expressions of Klos et al. (arXiv:1304.7684, App. B)"""
    rho = 0.10
    m = ce.M_PI
    k = (1.5 * np.pi**2 * rho * ce.HBARC**3)**(1. / 3.)
    I1, I2, IP, Ic6 = ce.fermi_integrals(rho, 0.)
    assert I1[0] == pytest.approx(1 - 3 * m**2 / k**2 + 3 * m**3 / k**3 * np.arctan(k / m), rel=1e-10)
    assert I2[0] == pytest.approx(I1[0], rel=1e-10)
    assert abs(IP[0]) < 1e-12 and abs(Ic6[0]) < 1e-12
    for Q in (0.05, 0.2, 0.5):
        at = np.arctan2(2 * m * k, m**2 + Q**2 / 4 - k**2)
        lg = np.log((m**2 + (k - Q / 2)**2) / (m**2 + (k + Q / 2)**2))
        I1a = (8 * k * Q * (48 * (k**2 + m**2)**2 + 32 * (k**2 - 3 * m**2) * Q**2 - 3 * Q**4) + 768 * m**3 * Q**3 * at
               + 3 * (16 * (k**2 + m**2)**2 - 8 * (k**2 - 5 * m**2) * Q**2 + Q**4) * (4 * (k**2 + m**2) - Q**2) * lg) \
            / (512 * k**3 * Q**3)
        I2a = (8 * k * (2 * k**2 - 3 * m**2) * Q + 24 * m**3 * Q * at + 3 * m**2 * (4 * k**2 - Q**2 + 4 * m**2) * lg) \
            / (16 * k**3 * Q)
        I1, I2, _, _ = ce.fermi_integrals(rho, Q)
        assert I1[0] == pytest.approx(I1a, rel=1e-10)
        assert I2[0] == pytest.approx(I2a, rel=1e-10)


def test_deltas_at_q0():
    """delta'(0) = delta''(0) = delta a(0), the reduction of the isovector coupling fixed by c_D"""
    for name, (rho, da0) in ce.TWO_BODY_SETTINGS.items():
        d1, d2 = ce.sd_deltas(np.array([0.]), name)
        assert d1[0] == pytest.approx(da0, abs=1e-10) and d2[0] == pytest.approx(da0, abs=1e-10)
    d1, d2 = ce.sd_deltas(np.array([0., 0.1]), "none")
    assert not d1.any() and not d2.any()
    # interpolation of the cached table vs direct evaluation
    q = np.array([0.0173, 0.1234, 0.3])
    rho, da0 = ce.TWO_BODY_SETTINGS["central"]
    da, daP = ce.delta_a(q, rho, ce.c_D_from_reduction(rho, da0))
    d1, d2 = ce.sd_deltas(q, "central")
    assert np.allclose(d1, -q**2 * ce.R2_A / 6 + da, atol=1e-5)
    assert np.allclose(d2, -ce.G_PI_NN * ce.F_PION / (ce.G_A * ce.M_N) * q**2 / (q**2 + ce.M_PI**2) + da + daP, atol=1e-5)


@pytest.mark.parametrize("target,basis,tol,tol0", [("Xe", "iso_GCN5082", 0.02, 1e-3), ("Ge", "pn_Fitz", 0.025, 0.02),
                                                    ("Ar", "pn_Fitz", 0.005, 1e-3), ("F", "pn_Fitz", 0.005, 1e-3)])
def test_standard_SI_matches_one_body_NREFT(target, basis, tol, tol0):
    """Standard SI (c^M only): the chiral EFT and the NREFT shell-model O_1 rates of RAPIDD agree within 2%
    (shell-model shapes); at q = 0 both are normalized to A^2, except the Ge responses of pn_Fitz, whose F_M(0) are
    the values of arXiv:1203.3542 rounded to 2 digits (e.g. 1000 instead of Z^2 = 1024 for 70Ge p p)"""
    E = np.array([0.5, 5., 20., 60.])
    c = 1e-8
    rapidd.set_any_Ncoeff(c * vev**2, 1, 'p')
    rapidd.set_any_Ncoeff(c * vev**2, 1, 'n')
    r_nreft = rapidd.difrate_dER(0.3, 100., E, target=target, basis=basis)
    r_chiral = ce.difrate_chiral(0.3, 100., E, ce.couplings_standard_SI(c), target=target)
    assert np.allclose(r_chiral / r_nreft, 1., atol=tol)
    assert r_chiral[0] / r_nreft[0] == pytest.approx(1., abs=tol0)


def test_axial_and_tensor():
    """TT = AA without two-body currents; the two-body currents only act on the isovector coupling"""
    E = np.array([1., 10., 40.])
    c = 1e-6
    for a in ((c, 0.), (0., c), (c, -0.5 * c)):
        aa_1b = ce.difrate_chiral(0.3, 50., E, ce.ChiralCouplings(c4A=a), "Xe", two_body="none")
        tt = ce.difrate_chiral(0.3, 50., E, ce.ChiralCouplings(c4T=a), "Xe")
        assert np.allclose(aa_1b, tt, rtol=1e-12)
    iso = ce.ChiralCouplings(c4A=(c, c))
    assert np.allclose(ce.difrate_chiral(0.3, 50., E, iso, "Xe"), ce.difrate_chiral(0.3, 50., E, iso, "Xe", "none"))
    # neutron-only coupling in xenon: the two-body currents reduce the rate (HMS Fig. 3)
    n = ce.ChiralCouplings(c4A=(0., c))
    ratio = ce.difrate_chiral(0.3, 50., E, n, "Xe") / ce.difrate_chiral(0.3, 50., E, n, "Xe", "none")
    assert np.all((ratio > 0.6) & (ratio < 0.8))
    # proton-only coupling in xenon: dominated by the two-body currents
    p = ce.ChiralCouplings(c4A=(c, 0.))
    assert np.all(ce.difrate_chiral(0.3, 50., E, p, "Xe") > 10 * ce.difrate_chiral(0.3, 50., E, p, "Xe", "none"))
    # no SD rate in argon (J = 0)
    assert not ce.difrate_chiral(0.3, 50., E, n, "Ar").any()


def test_mdm_couplings_reproduce_c1_c4():
    from rapidd.alpha_map import chiral_couplings_mdm, c1_dirac_mdm, c4_dirac_mdm
    alphaq = {'alpha_u': (1e-9, 2e-9, 1e-6, -3e-7), 'alpha_d': (2e-9, -1e-9, 5e-7, 2e-6), 'alpha_s': (3e-9, 0, 1e-7, 0),
              'alpha_c': (1e-8, 0, 0, 0), 'alpha_b': (-1, -1, -1, -1)}
    c = chiral_couplings_mdm(CARD, alphaq, 30.)
    c1 = c1_dirac_mdm(CARD, alphaq)
    c4 = c4_dirac_mdm(CARD, alphaq)
    assert c.cM[0] + c.cM[1] == pytest.approx(c1[0]) and c.cM[0] - c.cM[1] == pytest.approx(c1[1])
    assert c.c4A[0] + c.c4T[0] == pytest.approx(c4[0]) and c.c4A[1] + c.c4T[1] == pytest.approx(c4[1])
    assert c.cpi != 0. and c.cb != 0. and c.cPhi[0] != 0. and c.cMdot[0] != 0.


def test_helm_form_factor():
    """Helm form factor of Lewin & Smith (1996): F^2(0) = 1, first zero at q r_n = 4.493 (j_1(x) = 0)"""
    A = 131
    c, a, s = 1.23 * A**(1. / 3.) - 0.60, 0.52, 0.9
    rn = np.sqrt(c**2 + 7. / 3. * np.pi**2 * a**2 - 5. * s**2) / ce.HBARC
    assert ce.helm_form_factor_sq(0., A) == pytest.approx(1.)
    assert ce.helm_form_factor_sq(4.4934094579 / rn, A) < 1e-15
    # small q: F^2 ~ 1 - q^2 (r_n^2/5 + s^2) (Gaussian smearing)
    q = 1e-3
    assert ce.helm_form_factor_sq(q, A) == pytest.approx(1. - q**2 * (rn**2 / 5. + 2 * (0.9 / ce.HBARC)**2 / 2.), rel=1e-6)


def test_fermi_integrals_against_hu2022_closed_forms():
    """I^P and I_c6 by quadrature vs the closed forms of the Hu2022 notebook (DM_func.delta_2BC); I_1, I_2 are
    checked in test_fermi_integrals_against_analytic_expressions"""
    rho, m = 0.10, ce.M_PI
    k = (1.5 * np.pi**2 * rho * ce.HBARC**3)**(1. / 3.)
    for Q in (0.02, 0.1, 0.3, 0.55):
        lg = np.log((4 * m**2 + (2 * k - Q)**2) / (4 * m**2 + (2 * k + Q)**2))
        IP = -3. / (512 * k**3 * Q**3) * (8 * k * Q * (48 * (k**2 + m**2)**2 - 32 * k**2 * Q**2 - 3 * Q**4)
                                         + 3 * (4 * (k**2 + m**2) - Q**2) * (4 * m**2 + (2 * k - Q)**2)
                                         * (4 * m**2 + (2 * k + Q)**2) * lg)
        Ic6 = -9. / (128 * k**3 * Q) * (32 * k**3 * Q + 32 * k * m**2 * Q + 8 * k * Q**3
                                        + (16 * (k**2 + m**2)**2 + 8 * (m**2 - k**2) * Q**2 + Q**4) * lg)
        _, _, IPq, Ic6q = ce.fermi_integrals(rho, Q)
        assert IPq[0] == pytest.approx(IP, rel=1e-8)
        assert Ic6q[0] == pytest.approx(Ic6, rel=1e-8)


def test_hu2022_deltas_reproduce_their_notebook():
    """delta', delta'' of Hu2022 Eq. 3 vs DM_func.delta_2BC + dSigmap_1B / dSigmapp_1B of the AbInitio4DarkMatter
    notebook (commit e10cccb), evaluated with its constants and LECs: (interaction, rho [fm^-3], q [MeV]) ->
    (delta', delta'')"""
    ref = {("EM1.8_2.0", 0.08, 0.): (-0.168014, -0.168014), ("EM1.8_2.0", 0.08, 50.): (-0.176290, -0.239818),
           ("N4LO_LNL2", 0.12, 50.): (-0.381233, -0.370618), ("N4LO_LNL2", 0.12, 300.): (-0.652698, -0.782926),
           ("DNNLOgo", 0.12, 150.): (-0.268858, -0.516348), ("DNNLOgo", 0.08, 0.): (-0.119977, -0.119977)}
    for (inte, rho, q), (r1, r2) in ref.items():
        d1, d2 = ce.hu_deltas(np.array([q * 1e-3]), inte, rho)
        assert d1[0] == pytest.approx(r1, abs=2e-6) and d2[0] == pytest.approx(r2, abs=2e-6), (inte, rho, q)


def test_hu2022_normalization_at_q0():
    """One-body Hu2022 S_N(0) = (2J+1)(J+1)/(pi J) <S_N>^2 with the spin expectation values of their notebook (the
    same normalization as HMS): within 0.5% (the fits are not constrained at q = 0) or 2e-5 absolute"""
    for name, (A, Z, J) in ce.sd_isotopes("Hu2022").items():
        for inte in ce.HU2022_INTERACTIONS:
            S_p, S_n = ce.spin_expectation(A, Z, "Hu2022", inte)
            for aA, S in (((1., 0.), S_p), ((0., 1.), S_n)):
                val = ce.sd_structure_function(A, Z, 0., aA=aA, two_body=(inte, None), sd_model="Hu2022")
                ref = (2 * J + 1) * (J + 1) / (np.pi * J) * S**2
                assert val == pytest.approx(ref, rel=5e-3, abs=2e-5), (name, inte, aA)


def test_bands_contain_central():
    E = np.linspace(0.5, 60., 25)
    for c in (ce.ChiralCouplings(c4A=(0., 1e-6)), ce.ChiralCouplings(c4A=(1e-6, 0.)),
              ce.ChiralCouplings(c4A=(1e-6, -1e-6)) + ce.couplings_standard_SI(1e-9)):
        for model in ce.SD_MODELS:
            central, lo, hi = ce.difrate_chiral_band(0.3, 100., E, c, "Xe", sd_model=model)
            assert np.all(lo <= central * (1 + 1e-12)) and np.all(central <= hi * (1 + 1e-12)), model
            assert np.all(hi > lo)
    # the SI part has no band
    si = ce.couplings_standard_SI(1e-9)
    central, lo, hi = ce.difrate_chiral_band(0.3, 100., E, si, "Xe", sd_model="Hu2022")
    assert np.array_equal(lo, central) and np.array_equal(hi, central)


def test_hu2022_tensor_not_allowed():
    with pytest.raises(ValueError, match="axial-vector"):
        ce.difrate_chiral(0.3, 100., 10., ce.ChiralCouplings(c4T=(0., 1e-6)), sd_model="Hu2022")
    with pytest.raises(ValueError, match="axial-vector"):
        ce.sd_structure_function(129, 54, 0.05, aT=(0., 1.), two_body=("EM1.8_2.0", 0.1), sd_model="Hu2022")
    with pytest.raises(ValueError):
        ce.difrate_chiral(0.3, 100., 10., ce.ChiralCouplings(c4A=(0., 1e-6)), sd_model="Klos2013")


def test_hms_central_within_hu2022_band():
    """S_n(q) of 129Xe with two-body currents: the Hu2022 central value lies in the Hu2022 band (three interactions,
    rho = 0.08 ... 0.12 fm^-3) for q = 0 ... 200 MeV, the HMS central value only for q <= 55 MeV: above, HMS is
    below the Hu2022 band (0.94 times its lower edge at 60 MeV, 0.60 at 120 MeV, back to ~0.95 at 170-200 MeV).
    For 131Xe the HMS S_n is above the Hu2022 band at all q (1.08-1.32 times its upper edge), already at q = 0
    (<S_n> = -0.272 in HMS, -0.194 ... -0.239 in Hu2022)."""
    q = np.linspace(0., 0.2, 21)
    hu = np.array([ce.sd_structure_function(129, 54, q, aA=(0., 1.), two_body=(i, r), sd_model="Hu2022")
                   for i in ce.HU2022_INTERACTIONS for r in ce.HU2022_RHO_BAND])
    hu_central = np.mean([ce.sd_structure_function(129, 54, q, aA=(0., 1.), two_body=(i, ce.HU2022_RHO_CENTRAL),
                                                   sd_model="Hu2022") for i in ce.HU2022_INTERACTIONS], axis=0)
    hms = ce.sd_structure_function(129, 54, q, aA=(0., 1.), two_body="central")
    low = q <= 0.055
    assert np.all((hms[low] >= hu.min(axis=0)[low]) & (hms[low] <= hu.max(axis=0)[low]))
    assert np.all(hms[~low & (q <= 0.15)] < hu.min(axis=0)[~low & (q <= 0.15)])
    assert np.all((hu_central >= hu.min(axis=0)) & (hu_central <= hu.max(axis=0)))


def test_equivalent_sigma_reference_models():
    """Pure SD-n (SD-p) axial-vector model: sigma_eq SD-n (SD-p) = sigma_SD_n (sigma_SD_p) of MadDM; standard SI
    isoscalar model: sigma_eq SI / sigma_SI_p = N(HKMS)/N(Helm), a few per cent"""
    from rapidd.alpha_map import sigmaSD_nucleon_mdm, sigmaSI_nucleon_mdm
    m = 100.
    for alphaq, key, idx in (({'alpha_u': (0., 0., 1e-6, 0.), 'alpha_d': (0., 0., 1e-6 * 0.427 / 0.842, 0.)}, 'SD_p', 0),
                             ({'alpha_u': (0., 0., 1e-6 * 0.427 / 0.842, 0.), 'alpha_d': (0., 0., 1e-6, 0.)}, 'SD_n', 1)):
        card = dict(CARD, AVPs=0., AVNs=0.)
        c = ce.ChiralCouplings(c4A=(1e-6, 0.)) if idx == 0 else ce.ChiralCouplings(c4A=(0., 1e-6))
        s = ce.equivalent_sigma(c, m)
        ref = 3 * (m * (ce.mproton, ce.mneutron)[idx] / (m + (ce.mproton, ce.mneutron)[idx]))**2 * 1e-12 \
            / (16 * np.pi) * ce.GEV2_TO_CM2
        assert s[key] == pytest.approx(ref, rel=1e-6)
    # with the MadDM alphas (AA only, the proton-only or neutron-only combination of u and d)
    for alphaq, key, idx in (({'alpha_u': (0., 0., 0.842e-6, 0.), 'alpha_d': (0., 0., 0.427e-6, 0.)}, 'SD_p', 0),
                             ({'alpha_u': (0., 0., 0.427e-6, 0.), 'alpha_d': (0., 0., 0.842e-6, 0.)}, 'SD_n', 1)):
        from rapidd.calc_dRdE import equivalent_sigma_mdm
        card = dict(CARD, AVPs=0., AVNs=0.)
        res = equivalent_sigma_mdm(card, alphaq, m)
        sig = sigmaSD_nucleon_mdm(card, alphaq, m)
        assert sig[1 - idx] < 1e-6 * sig[idx]
        assert res['sigmaN_%s_eq' % key] == pytest.approx(sig[idx], rel=1e-6)
    c = 1e-8
    s = ce.equivalent_sigma(ce.couplings_standard_SI(c), m)
    sig_p = (m * ce.mproton / (m + ce.mproton))**2 * c**2 / np.pi * ce.GEV2_TO_CM2
    assert 0.93 < s["SI"] / sig_p < 0.99
    assert s["SD_n"] > 0 and s["SD_p"] > s["SD_n"]   # the SI signal seen through the SD references


def test_equivalent_sigma_isospin_cancellation():
    """SI with f_n/f_p = -Z/(A-Z) of 132Xe: the xenon signal (nearly) cancels, sigma_eq << sigma_SI_p"""
    c = 1e-8
    fn = -54. / 78. * c
    s = ce.equivalent_sigma(ce.couplings_standard_SI(c, fn), 100.)
    sig_p = (100. * ce.mproton / (100. + ce.mproton))**2 * c**2 / np.pi * ce.GEV2_TO_CM2
    assert s["SI"] < 1e-2 * sig_p
    # proton-only SI: sigma_eq ~ (Z/A)^2 sigma_p
    s = ce.equivalent_sigma(ce.couplings_standard_SI(c, 0.), 100.)
    assert 0.8 * (54. / 131.)**2 < s["SI"] / sig_p < 1.1 * (54. / 131.)**2


def test_equivalent_sigma_independent_of_reference_and_density():
    from rapidd.calc_dRdE import equivalent_sigma_mdm
    c = ce.ChiralCouplings(c4A=(1e-6, -2e-6), c4T=(3e-7, 0.)) + ce.couplings_standard_SI(1e-9, 2e-9)
    a = ce.equivalent_sigma(c, 50.)
    b = ce.equivalent_sigma(c, 50., sigma_ref=1e-46)
    for k in ce.EQUIVALENT_SIGMA_KEYS:
        assert a[k] == pytest.approx(b[k], rel=1e-10)
    alphaq = {'alpha_u': (1e-9, 0., 1e-6, 0.), 'alpha_d': (1e-9, 0., 0., 2e-7)}
    r1 = equivalent_sigma_mdm(dict(CARD, rhoDM=0.3), alphaq, 50.)
    r2 = equivalent_sigma_mdm(dict(CARD, rhoDM=0.55), alphaq, 50.)
    assert set(r1) == {'sigmaN_SI_eq', 'sigmaN_SD_n_eq', 'sigmaN_SD_p_eq'}
    for k in r1:
        assert r1[k] == pytest.approx(r2[k], rel=1e-12)


def test_equivalent_sigma_no_signal_and_threshold():
    s = ce.equivalent_sigma(ce.ChiralCouplings(), 100.)
    assert s == {"SI": 0., "SD_n": 0., "SD_p": 0.}
    # below the LZ threshold the reference models have no signal: no verdict
    s = ce.equivalent_sigma(ce.couplings_standard_SI(1e-8), 2.)
    assert all(np.isnan(v) for v in s.values())


def test_DDrate_save_chiral(tmp_path):
    """DDrate_save with the chiral EFT structure factors (default): total = SI + SD for Xe, Ar, Ge"""
    import os
    from rapidd.calc_dRdE import DDrate_save
    from rapidd.alpha_map import chiral_couplings_mdm
    alphaq = {'alpha_u': (1e-9, 0., 1e-6, 1e-6), 'alpha_d': (1e-9, 0., 0., 0.)}
    assert DDrate_save(CARD, alphaq, 50., str(tmp_path), targets=["xe", "ar", "ge"], save_npz=True) is None
    d = np.load(os.path.join(tmp_path, 'dRdE_nucleon.npz'))
    for t in ('Xe', 'Ar', 'Ge'):
        np.testing.assert_allclose(d['drde_' + t], d['drde_%s_SI' % t] + d['drde_%s_SD' % t], rtol=1e-12)
    assert not d['drde_Ar_SD'].any() and d['drde_Xe_SD'].all()
    E = d['energies'][[10, 300, 600]]
    c = chiral_couplings_mdm(CARD, alphaq, 50.)
    np.testing.assert_allclose(d['drde_Xe'][[10, 300, 600]], ce.difrate_chiral(0.3, 50., E, c, 'Xe'), rtol=1e-12)
    txt = np.loadtxt(os.path.join(tmp_path, 'dRdE_nucleon.dat'))
    np.testing.assert_allclose(txt[:, 1], d['drde_Xe'], rtol=1e-6)
    # one-body NREFT rates
    os.makedirs(tmp_path / 'nreft')
    DDrate_save(CARD, alphaq, 50., str(tmp_path / 'nreft'), nuclear_model='nreft', save_npz=True)
    assert 'drde_Xe_SI' not in np.load(os.path.join(tmp_path, 'nreft', 'dRdE_nucleon.npz'))
    with pytest.raises(ValueError):
        DDrate_save(CARD, alphaq, 50., str(tmp_path), nuclear_model='klos')
