'''Comparison with the figures of arXiv:1912.08204 (data extracted with
tests/data/extract_paper_figures.py into tests/data/paper_1912.08204).

Settings of the paper: SHM with v0 = 220, vesc = 544, vE = 244 km/s, rho = 0.4 GeV/cm^3, jchi = 1/2,
the 5 outer Xenon shells (4s, 4p, 4d, 5s, 5p), dR/dlnE in kg^-1 yr^-1, g_e = 2.

Findings (see the README / report of the electron_nreft branch):
- The XENON10 limits (Figs. 8, 10) agree within 1% at low masses and within 25% everywhere, and
  the model spectra (Fig. 9) within ~20%.
- The spectra of Figs. 5 and 7 and the integrated rates of Fig. 6 are ~1.1-1.45 times smaller in
  the paper than Eq. (ionization spectrum) with the paper's own atomic responses, halo and
  couplings (RAPIDD agrees with an independent evaluation of that equation to 1e-4, and its W_1 of
  5p with Fig. 2 of the paper). The offset depends on the kinematics (~1.45 when the rate comes
  from the tail of the velocity distribution, ~1.1 in the bulk), so it cannot be a convention
  factor: these comparisons are expected failures (xfail) at the 20% level, and are checked
  against the observed band instead.
- Fig. 6, O_3: for jchi = 1/2 the DM responses of the paper give exactly the same rate for O_3
  and O_5, while the paper's O_3 bars are 5-50 times larger than O_5 (and than RAPIDD).
- All the comparisons use the sign of W_2 of the paper (set_w2_sign("Catena")); arXiv:2405.04855
  corrects it (RAPIDD default) and also finds Fig. 7 to be ~3/4 of the formula.
- XENON1T: the paper uses 4 S2 bins with g2 = 33, RAPIDD the S2-only data release response and
  its mass dependent ROI: only the order of magnitude is compared.
'''
import os

import numpy as np
import pytest

import rapidd

from conftest import load_csv, PAPER_DIR, PAPER_HALO

YR = 365.25


@pytest.fixture(autouse=True, scope="module")
def catena_w2_sign():
    '''The figures of arXiv:1912.08204 use its sign of W_2'''
    rapidd.set_w2_sign("Catena")
    yield
    rapidd.set_w2_sign("Liang")
RHO = 0.4
SHELLS = rapidd.XE_SHELLS_PAPER
OFFSET_REASON = ("Figs. 5-7 of arXiv:1912.08204 are 1.1-1.45 times below Eq. (ionization spectrum) "
                 "evaluated with the paper's inputs (kinematics dependent offset, see module docstring)")


def spectrum(m, coeffs, resp=None):
    rapidd.reset_coefficients()
    for (op, inter), val in coeffs.items():
        rapidd.set_any_Ecoeff(val, op, interaction=inter)
    r = rapidd.difrate_electron(RHO, m, 0.5, shells=SHELLS, resp=resp)
    return r['E'] * 1e3, r['dRdlnE'] * YR   # eV, kg^-1 yr^-1


def ratio(name, E, y, threshold=1e-3):
    '''RAPIDD / paper at the paper's points where the paper's rate is > threshold * peak'''
    x, yp = load_csv(name).T
    ok = y != 0
    yi = np.exp(np.interp(np.log(x), np.log(E[ok]), np.log(np.abs(y[ok]))))
    sel = yp > threshold * yp.max()
    return yi[sel] / yp[sel]


# (csv name, mDM [GeV], coefficients, resp)
C13, L13 = {(13, "contact"): 1e-3}, {(13, "long-range"): 1e-3}
C7, L15 = {(7, "contact"): 1e-3}, {(15, "long-range"): 1e-3}
SPECTRA = [("Xe_coup_13_C__m%s" % lab, m, C13, None) for lab, m in (("1000MeV", 1.), ("100MeV", .1), ("10MeV", .01))] \
    + [("Xe_coup_13_L__m%s" % lab, m, L13, None) for lab, m in (("1000MeV", 1.), ("100MeV", .1), ("10MeV", .01))] \
    + [("Xe_coup_7_C_mDM_10__%s" % lab, .01, C7, resp) for lab, resp in (("total", None), ("res1", 1), ("res2", 2), ("res3", 3))] \
    + [("Xe_coup_15_L_mDM_100__%s" % lab, .1, L15, resp) for lab, resp in (("total", None), ("res1", 1), ("res3", 3), ("res4_negative", 4))]


@pytest.mark.xfail(strict=True, reason=OFFSET_REASON)
@pytest.mark.parametrize("name,m,coeffs,resp", SPECTRA, ids=[s[0] for s in SPECTRA])
def test_spectra_figs_5_7(paper_halo, name, m, coeffs, resp):
    E, y = spectrum(m, coeffs, resp)
    r = ratio(name, E, y)
    assert np.all(np.abs(r - 1) < 0.2)


@pytest.mark.parametrize("name,m,coeffs,resp", SPECTRA, ids=[s[0] for s in SPECTRA])
def test_spectra_figs_5_7_offset_band(paper_halo, name, m, coeffs, resp):
    '''Guard: shape and sign as in the paper, normalization within the observed offset band'''
    E, y = spectrum(m, coeffs, resp)
    if name.endswith("res4_negative"):
        assert np.all(y <= 0)   # dagger in Fig. 7: negative contribution
    r = ratio(name, E, y)
    # the total of O_15 is a cancellation between R_3 W_3 and R_4 W_4: it amplifies the offset
    lo = 0.55 if name == "Xe_coup_15_L_mDM_100__total" else 0.95
    assert np.all((r > lo) & (r < 1.4)), r


MODELS = [("Spectrum_Anapole", "anapole", 1e-2, (8, 9)),
          ("Spectrum_Magnetic_Dipole", "magnetic_dipole", 1e-5, (1, 4, 5, 6)),
          ("Spectrum_Electric_Dipole", "electric_dipole", 1e-6, (11,))]
MODEL_CURVES = [(fig, model, g, ops) for fig, model, g, all_ops in MODELS for ops in (None,) + tuple((o,) for o in all_ops)]


@pytest.mark.parametrize("fig,model,g,ops", MODEL_CURVES,
                         ids=["%s_%s" % (c[0], "total" if c[3] is None else "O%d" % c[3][0]) for c in MODEL_CURVES])
def test_model_spectra_fig_9(paper_halo, fig, model, g, ops):
    '''Fig. 9 (mDM = 100 MeV): median within 25%, every point within 50% (the largest deviations
    are close to the end point, e.g. O_5 of the magnetic dipole)'''
    rapidd.reset_coefficients()
    rapidd.set_any_E_EFT_coeff(model, 0.1, g, operators=ops)
    r = rapidd.difrate_electron(RHO, 0.1, 0.5, shells=SHELLS)
    name = "%s__%s" % (fig, "total" if ops is None else "O%d" % ops[0])
    rr = ratio(name, r['E'] * 1e3, r['dRdlnE'] * YR)
    assert abs(np.median(rr) - 1) < 0.25, np.median(rr)
    assert np.all(np.abs(rr - 1) < 0.5), rr


def total_rate_rows():
    '''Xenon bars of Fig. 6: (interaction, operator, mDM [MeV], R/c^2 [kg^-1 yr^-1])'''
    data = np.genfromtxt(os.path.join(PAPER_DIR, "Total_Rate_Comparison.csv"), delimiter=",", dtype=None,
                         encoding=None, comments="#")
    return [(str(r[0]), int(r[1]), float(r[3]), float(r[4])) for r in data if str(r[2]) == "Xe"]


@pytest.fixture(scope="module")
def total_rate_ratios():
    rapidd.read_halo(PAPER_HALO)
    rapidd.read_atomic_responses()
    out = {}
    for inter, op, m, val in total_rate_rows():
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(1.0, op, interaction=inter)
        out[(inter, op, m)] = rapidd.total_rate_electron(RHO, m / 1e3, 0.5, shells=SHELLS) * YR / val
    rapidd.reset_coefficients()
    rapidd.read_halo(rapidd.halo_path)
    return out


@pytest.mark.xfail(strict=True, reason=OFFSET_REASON)
def test_total_rates_fig_6(total_rate_ratios):
    r = np.array([v for (inter, op, m), v in total_rate_ratios.items() if op != 3])
    assert np.all(np.abs(r - 1) < 0.2)


def test_total_rates_fig_6_offset_band(total_rate_ratios):
    '''Guard: all the operators but O_3, O_5, O_15 within the observed offset band'''
    for (inter, op, m), v in total_rate_ratios.items():
        if op in (3, 5, 15):
            continue
        assert 0.95 < v < 1.4, (inter, op, m, v)


@pytest.mark.xfail(strict=True, reason="Fig. 6: the paper's O_3 bars are 5-50 times above O_5, while the DM "
                                       "responses of the paper give the same rate for O_3 and O_5 (jchi = 1/2)")
def test_total_rate_O3_fig_6(total_rate_ratios):
    for inter in ("contact", "long-range"):
        for m in (10., 100., 1000.):
            assert abs(total_rate_ratios[(inter, 3, m)] - 1) < 0.2


def test_O3_equals_O5(paper_halo):
    '''For jchi = 1/2 the DM responses of O_3 and O_5 coincide (Eqs. (DM response 1-4))'''
    for inter in ("contact", "long-range"):
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(1.0, 3, interaction=inter)
        r3 = rapidd.difrate_electron(RHO, 0.1)['dRdE']
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(1.0, 5, interaction=inter)
        r5 = rapidd.difrate_electron(RHO, 0.1)['dRdE']
        np.testing.assert_allclose(r3, r5, rtol=1e-12)


SETTERS = {"c1": lambda inter: (lambda m: rapidd.set_any_Ecoeff(1., 1, interaction=inter)),
           "c8": lambda inter: (lambda m: rapidd.set_any_Ecoeff(1., 8, interaction=inter)),
           "c13": lambda inter: (lambda m: rapidd.set_any_Ecoeff(1., 13, interaction=inter))}
LIMITS = [("Constraint_Coupling_%s__%s_%s" % (tag, exp, c), exp, SETTERS[c](inter))
          for tag, inter in (("C", "contact"), ("L", "long-range")) for exp in ("XENON10", "XENON1T")
          for c in ("c1", "c8", "c13")] \
    + [("Constraint_%s__%s" % (fig, exp), exp, (lambda model: (lambda m: rapidd.set_any_E_EFT_coeff(model, m, 1.)))(model))
       for fig, model in (("Anapole", "anapole"), ("Magnetic_Dipole", "magnetic_dipole"),
                          ("Electric_Dipole", "electric_dipole")) for exp in ("XENON10", "XENON1T")]


def limit_ratios(name, exp, setter):
    x, y = load_csv(name).T
    sel = slice(None, None, 2)
    lim = rapidd.electron_limits_vs_mass(exp, x[sel] / 1e3, setter, rhoDM=RHO, shells=SHELLS)
    return x[sel], lim / y[sel]


@pytest.mark.parametrize("name,exp,setter", [l for l in LIMITS if l[1] == "XENON10"],
                         ids=[l[0] for l in LIMITS if l[1] == "XENON10"])
def test_xenon10_limits(paper_halo, name, exp, setter):
    '''Figs. 8 and 10, XENON10: within 30%'''
    m, r = limit_ratios(name, exp, setter)
    assert np.all(np.isfinite(r)), (m, r)
    assert np.all(np.abs(r - 1) < 0.3), (m, r)


@pytest.mark.parametrize("name,exp,setter", [l for l in LIMITS if l[1] == "XENON1T"],
                         ids=[l[0] for l in LIMITS if l[1] == "XENON1T"])
def test_xenon1t_limits_order_of_magnitude(paper_halo, name, exp, setter):
    '''Figs. 8 and 10, XENON1T: different analysis (S2-only release + ROI vs 4 bins, g2 = 33),
    only the order of magnitude is compared, for m > 25 MeV (the ROI table starts at 20 MeV, where
    the curves are steep and the two analyses differ by up to a factor ~12)'''
    m, r = limit_ratios(name, exp, setter)
    r = r[np.isfinite(r) & (m > 25.)]
    assert len(r) > 3
    assert np.all((r > 0.1) & (r < 10.)), r
