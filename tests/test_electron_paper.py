'''Comparison with the figures of arXiv:1912.08204 (data extracted with
tests/data/extract_paper_figures.py into tests/data/paper_1912.08204).

Settings of the paper: SHM with v0 = 220, vesc = 544, vE = 244 km/s, rho = 0.4 GeV/cm^3, jchi = 1/2,
the 5 outer Xenon shells (4s, 4p, 4d, 5s, 5p) and all the Argon shells (1s, 2s, 2p, 3s, 3p), dR/dlnE
in kg^-1 yr^-1, g_e = 2.

Atomic responses: the paper used DarkARC (Python, v1.0 of 2019, no longer maintained: q in [1, 1000] keV
on a 100 x 100 grid, l' <= 7 fixed, radial integrals with tolerances 1e-3/1e-2). The RAPIDD tables come
from DarkART 594ffb5 (C++, its maintained successor; lib/atomic_responses/README.md), converged in
r_points (5000 vs 10000: identical) and in l' (l' <= 40 fixed: W within 0.5%, transverse part
(q/m_e)^2 W_3 - W_4 within 0.25%). Where the two differ, the DarkART/RAPIDD results are the more reliable.

Findings (see examples/checks_ER_rates.ipynb):
- Figs. 5-7 use rho = 0.3 GeV/cm^3, not the 0.4 stated in the text: the spectra of Figs. 5 and 7 and
  the integrated rates of Fig. 6 are 0.752 times RAPIDD at rho = 0.4 (0.75 = 0.3/0.4), flat in E_e
  and the same for all the operators, masses and responses, while Figs. 8-10 (rho = 0.4) agree with
  RAPIDD at rho = 0.4. arXiv:2405.04855 also finds Fig. 7 to be ~3/4 of the formula. The values of
  Figs. 5-7 are therefore multiplied by RHO_FIX = 0.4/0.3 before the comparison: they then agree
  within 2% (0.5% for most curves and bars).
- The XENON10 limits (Figs. 8, 10), with the paper's sign of W_2, agree within 4% at all the
  extracted points (4-850 MeV): RAPIDD/paper is 0.99-1.04. Contact couplings at m > 100 MeV agree
  within 1%; at m < 20 MeV and for the long-range couplings, where the lowest S2 bin dominates, RAPIDD
  is 3.7% weaker. With the legacy left sum in E_e of the XENON10 response (dR/dE(E_i) (E_{i+1} - E_i)
  on the log grid of E_e, ~7.5% more rate, replaced by the integral of the interpolated spectrum in
  each bin of the number of quanta) RAPIDD/paper was 0.97-1.00, with the opposite pattern. With the
  sign of W_2 of arXiv:2405.04855 the limits on c_8, c_13 and the anapole are 1.4-3.8 times weaker.
- The model spectra (Fig. 9, rho = 0.4) agree within 1% (median), 7% at the end points.
- The previous excess of the paper at large q (contact interactions at 100 MeV - 1 GeV, amplified
  for O_5, O_15 by the cancellation between R_3 W_3 and R_4 W_4) came from the RAPIDD tables of Xe
  4d, produced with DarkART before its fix of the 4d binding energy (Dec. 2022); with the regenerated
  tables it is gone.
- Contact interactions at 1 GeV (Xenon): with the 5 shells and q <= 1 MeV the bars of O_5, O_6, O_13,
  O_14 (and O_3 with +c_3^2/4) of the paper (times 0.4/0.3) are 6-9% above RAPIDD and O_15 38% above
  (LARGE_Q_1GEV; the other bars at 1 GeV agree within 1.5%). These rates come from large q (the
  integrand is still large at q = 1 MeV), and are explained by two features of the paper's
  computation (test_1gev_contact_explained_by_inner_shells):
  (1) the paper also includes the shells 3s, 3p, 3d (its text says 5 shells), which contribute at
      large q: with 8 shells and q <= 1 MeV the ratios become O_6 0.970, O_13/O_14 0.971, O_3+ 0.972,
      O_5 0.976, O_15 0.865, while O_1, O_4, O_7-O_12 change by < 1%;
  (2) the paper's q integral extends effectively up to q ~ 1.05 MeV (evaluated with DarkART tables
      extended beyond 1 MeV, generated outside the repository and not shipped, see
      lib/atomic_responses/README.md; these numbers are not recomputed by the tests): with 8 shells and q_max = 1.05 MeV the O_6, O_13/O_14, O_3+ bars and the
      O_13 spectrum at 1 GeV agree within 0.3%.
  The remaining O_5 (+2%) and O_15 (+11%) are the difference of R_3 W_3 and R_4 W_4, up to ~400 times
  larger than the result in Xenon (~1100 for Argon 2s at 1 MeV): most likely a numerical inaccuracy of
  DarkARC (for Argon at 1 GeV the paper's O_5 and O_15 bars are 2.1 and 13 times RAPIDD,
  test_total_rates_fig_6_argon). Excluded as causes: an extension of q to the kinematic end point
  (200 times above the paper), linear or nearest-node interpolation of the tables, the halo, the
  truncation l' <= 7 of DarkARC (which however explains 0.4% of the long-range O_7, O_8, O_12, O_13,
  O_14 bars). The interpolation of the tables is not the cause either (DarkART tables with 399 q
  points change these rates by < 0.04%).
- q_max = 1 MeV ~ 2 m_e stays the RAPIDD default (as in DarkARC, DarkART, obscura): above it the NREFT
  and the non-relativistic wave functions are not valid. For mDM > 0.3 GeV and contact O_3, O_5, O_6,
  O_13-O_15 RAPIDD issues a QMaxTruncationWarning (tests/test_electron_qmax.py).
- Argon (all the 5 shells, test_*_argon): the bars of Fig. 6 and the spectra of Figs. 5 and 7 agree
  within 1% at 10 and 100 MeV and for the long-range interactions (within 2% tested). At 1 GeV contact
  RAPIDD/paper is O_6 0.982, O_13/O_14 0.973, O_3+ 0.964 (0.98-1.00 with q_max = 1.045 MeV), and
  O_5 0.472, O_15 0.076 (the R_3 W_3 - R_4 W_4 cancellation, not changed by q_max): pinned in
  LARGE_Q_1GEV_AR.
- Fig. 6, O_3: for jchi = 1/2 the DM responses of the paper give exactly the same rate for O_3 and
  O_5 (both in arXiv:1912.08204 and in arXiv:2405.04855, where the c_3^2 and c_5^2 terms of a_1,
  a_2 coincide), while the paper's O_3 bars are 3-60 times larger than O_5. They are reproduced
  (as all the other bars) with the opposite sign of the c_3^2 term of
  R_4, i.e. R_3 W_3 + R_4 W_4 -> c_3^2/4 (q^2/m_e^2 W_3 + W_4) instead of
  c_3^2/4 (q^2/m_e^2 W_3 - W_4) = c_3^2/4 sum |q/m_e x f|^2 >= 0: the figure, not Eq. (DM response 4),
  has the wrong sign. Only the bars with the opposite sign are tested.
- All the comparisons use the sign of W_2 of the paper (set_w2_sign("Catena")); arXiv:2405.04855
  corrects it (RAPIDD default) and also finds Fig. 7 to be ~3/4 of the formula.
- XENON1T: the paper uses 4 S2 bins with g2 = 33, RAPIDD the S2-only data release response and
  its mass dependent ROI: only the order of magnitude is compared (RAPIDD is validated against the
  data release in test_electron_xenon1t.py).
'''
import os

import numpy as np
import pytest

import rapidd

# expected for the contact O_3, O_5, O_6, O_13-O_15 above 0.3 GeV (tested in test_electron_qmax.py)
pytestmark = pytest.mark.filterwarnings("ignore::rapidd.electron.QMaxTruncationWarning")

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
RHO_FIX = 0.4 / 0.3   # Figs. 5-7 use rho = 0.3: their values are multiplied by 0.4/0.3 (module docstring)
# Contact at 1 GeV, 5 shells, q <= 1 MeV: RAPIDD/paper (rescaled). The paper is 6-38% above because it
# also includes 3s, 3p, 3d and integrates effectively up to q ~ 1.05 MeV; O_5 and O_15 keep a residual
# from the R_3 W_3 - R_4 W_4 cancellation (module docstring, test_1gev_contact_explained_by_inner_shells).
LARGE_Q_1GEV = {"3+": 0.931, 5: 0.944, 6: 0.913, 13: 0.924, 14: 0.924, 15: 0.724}
# Same with the 8 shells 3s ... 5p (RAPIDD default q <= 1 MeV)
SHELLS_8 = rapidd.XE_SHELLS_PAPER + ("3s", "3p", "3d")
LARGE_Q_1GEV_8_SHELLS = {"3+": 0.972, 5: 0.976, 6: 0.970, 13: 0.971, 14: 0.971, 15: 0.865}
# Argon, all the shells, contact at 1 GeV: q_max (O_6, O_13, O_14, O_3+) and R_3 W_3 - R_4 W_4
# cancellation in the paper (O_5, O_15), module docstring
LARGE_Q_1GEV_AR = {"3+": 0.964, 5: 0.472, 6: 0.982, 13: 0.973, 14: 0.973, 15: 0.076}


def spectrum(m, coeffs, resp=None, shells=SHELLS, target="Xe"):
    rapidd.reset_coefficients()
    for (op, inter), val in coeffs.items():
        rapidd.set_any_Ecoeff(val, op, interaction=inter)
    r = rapidd.difrate_electron(RHO, m, 0.5, shells=shells, resp=resp, target=target)
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


@pytest.mark.parametrize("name,m,coeffs,resp", SPECTRA, ids=[s[0] for s in SPECTRA])
def test_spectra_figs_5_7(paper_halo, name, m, coeffs, resp):
    '''Shape and sign as in the paper, with the paper's curves rescaled by RHO_FIX = 0.4/0.3: within 2%
    (5% close to the end points of the 1 GeV spectra). O_13 contact at 1 GeV: RAPIDD/paper 0.92-0.93
    (pinned as the bar LARGE_Q_1GEV[13]), down to 0.59 at the end point: the paper also includes 3s,
    3p, 3d (median 0.980 with 8 shells, test_1gev_contact_explained_by_inner_shells) and integrates up
    to q ~ 1.05 MeV (with both, within 0.3% up to 1 keV, module docstring).'''
    E, y = spectrum(m, coeffs, resp)
    if name.endswith("res4_negative"):
        assert np.all(y <= 0)   # dagger in Fig. 7: negative contribution
    r = ratio(name, E, y) / RHO_FIX
    if name == "Xe_coup_13_C__m1000MeV":  # as the bar of O_13 contact at 1 GeV, down to 0.59 at the end point
        assert abs(np.median(r) - LARGE_Q_1GEV[13]) < 0.01 and np.all((r > 0.55) & (r < 0.95)), r
    else:
        assert np.all((r > 0.95) & (r < 1.02)), r


MODELS = [("Spectrum_Anapole", "anapole", 1e-2, (8, 9)),
          ("Spectrum_Magnetic_Dipole", "magnetic_dipole", 1e-5, (1, 4, 5, 6)),
          ("Spectrum_Electric_Dipole", "electric_dipole", 1e-6, (11,))]
MODEL_CURVES = [(fig, model, g, ops) for fig, model, g, all_ops in MODELS for ops in (None,) + tuple((o,) for o in all_ops)]


@pytest.mark.parametrize("fig,model,g,ops", MODEL_CURVES,
                         ids=["%s_%s" % (c[0], "total" if c[3] is None else "O%d" % c[3][0]) for c in MODEL_CURVES])
def test_model_spectra_fig_9(paper_halo, fig, model, g, ops):
    '''Fig. 9 (mDM = 100 MeV, rho = 0.4): median within 2%, every point within 10% (the largest
    deviations are close to the end point)'''
    rapidd.reset_coefficients()
    rapidd.set_any_E_EFT_coeff(model, 0.1, g, operators=ops)
    r = rapidd.difrate_electron(RHO, 0.1, 0.5, shells=SHELLS)
    name = "%s__%s" % (fig, "total" if ops is None else "O%d" % ops[0])
    rr = ratio(name, r['E'] * 1e3, r['dRdlnE'] * YR)
    assert abs(np.median(rr) - 1) < 0.02, np.median(rr)
    assert np.all(np.abs(rr - 1) < 0.1), rr


def total_rate_rows(target="Xe"):
    '''Xenon (or Argon) bars of Fig. 6: (interaction, operator, mDM [MeV], R/c^2 [kg^-1 yr^-1])'''
    data = np.genfromtxt(os.path.join(PAPER_DIR, "Total_Rate_Comparison.csv"), delimiter=",", dtype=None,
                         encoding=None, comments="#")
    return [(str(r[0]), int(r[1]), float(r[3]), float(r[4])) for r in data if str(r[2]) == target]


@pytest.fixture(scope="module")
def total_rate_ratios():
    '''RAPIDD/paper for the Xenon bars of Fig. 6, keyed by (interaction, operator, mDM [MeV]); for O_3
    also with the opposite sign of R_4 W_4 (key (interaction, "3+", mDM))'''
    rapidd.read_halo(PAPER_HALO)
    rapidd.read_atomic_responses()
    out = {}
    for inter, op, m, val in total_rate_rows():
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(1.0, op, interaction=inter)
        out[(inter, op, m)] = rapidd.total_rate_electron(RHO, m / 1e3, 0.5, shells=SHELLS) * YR / val
        if op == 3:
            r4 = rapidd.total_rate_electron(RHO, m / 1e3, 0.5, shells=SHELLS, resp=4) * YR / val
            out[(inter, "3+", m)] = out[(inter, op, m)] - 2 * r4
    rapidd.reset_coefficients()
    rapidd.read_halo(rapidd.halo_path)
    return out


def test_total_rates_fig_6(total_rate_ratios):
    '''All the operators but O_3, paper's bars rescaled by RHO_FIX = 0.4/0.3: within 2%, except the
    contact O_5, O_6, O_13, O_14, O_15 at 1 GeV, checked against their current values (LARGE_Q_1GEV)
    within 1%: the paper also includes 3s, 3p, 3d and integrates up to q ~ 1.05 MeV, plus the
    numerical R_3 W_3 - R_4 W_4 cancellation for O_5 and O_15 (module docstring and
    test_1gev_contact_explained_by_inner_shells)'''
    for (inter, op, m), v in total_rate_ratios.items():
        if op in (3, "3+"):
            continue
        r = v / RHO_FIX
        if (inter, m) == ("contact", 1000.) and op in LARGE_Q_1GEV:
            assert abs(r - LARGE_Q_1GEV[op]) < 0.01, (inter, op, m, r)
        else:
            assert 0.98 < r < 1.02, (inter, op, m, r)


def test_total_rate_O3_fig_6_opposite_R4_sign(total_rate_ratios):
    '''The paper's O_3 bars (3-60 times above O_5, while Eq. (DM response 4) gives O_3 = O_5 for
    jchi = 1/2) are reproduced, as the other bars, with the opposite sign of R_4 W_4, i.e. with
    +c_3^2/4 in R_4 (see module docstring). The comparison with the sign of the paper's equation
    (that of RAPIDD) is not tested.'''
    for inter in ("contact", "long-range"):
        for m in (10., 100., 1000.):
            r = total_rate_ratios[(inter, "3+", m)] / RHO_FIX
            if (inter, m) == ("contact", 1000.):
                assert abs(r - LARGE_Q_1GEV["3+"]) < 0.01, r
            else:
                assert 0.98 < r < 1.02, (inter, m, r)


def test_1gev_contact_explained_by_inner_shells(paper_halo):
    '''Contact at 1 GeV: with the 8 shells 3s ... 5p (the paper's text says 5) and q <= 1 MeV the
    ratios to the paper (rescaled) become LARGE_Q_1GEV_8_SHELLS (O_6 0.970, O_13/O_14 0.971, O_3+ 0.972,
    O_5 0.976, O_15 0.865) within 0.005, while O_1, O_4, O_7-O_12 change by < 1%; the O_13 spectrum
    goes from a median 0.928 to 0.980. The rest (q up to ~1.05 MeV in the paper, cancellation for
    O_5, O_15) is discussed in the module docstring.'''
    for inter, op, m, val in total_rate_rows():
        if (inter, m) != ("contact", 1000.):
            continue
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(1.0, op, interaction=inter)
        r5, r8 = (rapidd.total_rate_electron(RHO, 1., 0.5, shells=sh) * YR / val / RHO_FIX for sh in (SHELLS, SHELLS_8))
        if op == 3:
            r5 -= 2 * rapidd.total_rate_electron(RHO, 1., 0.5, shells=SHELLS, resp=4) * YR / val / RHO_FIX
            r8 -= 2 * rapidd.total_rate_electron(RHO, 1., 0.5, shells=SHELLS_8, resp=4) * YR / val / RHO_FIX
            op = "3+"
        if op in LARGE_Q_1GEV_8_SHELLS:
            assert abs(r5 - LARGE_Q_1GEV[op]) < 0.01, (op, r5)
            assert abs(r8 - LARGE_Q_1GEV_8_SHELLS[op]) < 0.005, (op, r8)
            assert r8 > r5 + 0.03
        else:
            assert abs(r8 / r5 - 1) < 0.01, (op, r5, r8)
            assert 0.98 < r8 < 1.02, (op, r8)
    meds = []
    for sh in (SHELLS, SHELLS_8):
        E, y = spectrum(1., C13, shells=sh)
        meds.append(np.median(ratio("Xe_coup_13_C__m1000MeV", E, y) / RHO_FIX))
    assert abs(meds[0] - 0.928) < 0.005 and abs(meds[1] - 0.980) < 0.005, meds


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


########################################################################################
# Argon (all the 5 shells, as in the paper)
########################################################################################

@pytest.fixture
def paper_halo_argon(paper_halo):
    '''Halo of arXiv:1912.08204 and the Argon atomic responses'''
    rapidd.read_atomic_responses(target="Ar")
    yield


@pytest.fixture(scope="module")
def total_rate_ratios_argon():
    '''RAPIDD/paper (rescaled by RHO_FIX) for the Argon bars of Fig. 6, keyed by (interaction, operator,
    mDM [MeV]); O_3 also with the opposite sign of R_4 W_4 (key (interaction, "3+", mDM))'''
    rapidd.read_halo(PAPER_HALO)
    rapidd.read_atomic_responses(target="Ar")
    out = {}
    for inter, op, m, val in total_rate_rows("Ar"):
        rapidd.reset_coefficients()
        rapidd.set_any_Ecoeff(1.0, op, interaction=inter)
        out[(inter, op, m)] = rapidd.total_rate_electron(RHO, m / 1e3, 0.5, target="Ar") * YR / val / RHO_FIX
        if op == 3:
            r4 = rapidd.total_rate_electron(RHO, m / 1e3, 0.5, target="Ar", resp=4) * YR / val / RHO_FIX
            out[(inter, "3+", m)] = out[(inter, op, m)] - 2 * r4
    rapidd.reset_coefficients()
    rapidd.read_halo(rapidd.halo_path)
    return out


def test_total_rates_fig_6_argon(total_rate_ratios_argon):
    '''Argon bars of Fig. 6 (times 0.4/0.3), O_3 with +c_3^2/4 as for Xenon: within 2% (within 1% in
    practice) at 10 and 100 MeV and for all the long-range bars. Contact at 1 GeV pinned within 1%
    (LARGE_Q_1GEV_AR): O_6 0.982, O_13/O_14 0.973, O_3+ 0.964 (q_max = 1 MeV, the paper integrates up
    to ~1.05 MeV), O_5 0.472 and O_15 0.076 (the paper's O_5, O_15 are 2.1 and 13 times RAPIDD: the
    R_3 W_3 - R_4 W_4 cancellation, ~1100 for 2s at 1 MeV, in DarkARC); O_1, O_4, O_7-O_12 within 2%.'''
    assert len(total_rate_ratios_argon) == 2 * 3 * (len(rapidd.ELECTRON_OPERATORS) + 1)
    for (inter, op, m), r in total_rate_ratios_argon.items():
        if op == 3:
            continue
        if (inter, m) == ("contact", 1000.) and op in LARGE_Q_1GEV_AR:
            assert abs(r / LARGE_Q_1GEV_AR[op] - 1) < 0.01, (inter, op, m, r)
        else:
            assert 0.98 < r < 1.02, (inter, op, m, r)


C6, L6 = {(6, "contact"): 1e-3}, {(6, "long-range"): 1e-3}
SPECTRA_AR = [("Ar_coup_6_C__m%s" % lab, m, C6, None) for lab, m in (("1000MeV", 1.), ("100MeV", .1), ("10MeV", .01))] \
    + [("Ar_coup_6_L__m%s" % lab, m, L6, None) for lab, m in (("1000MeV", 1.), ("100MeV", .1), ("10MeV", .01))] \
    + [("Ar_coup_7_C_mDM_10__%s" % lab, .01, C7, resp) for lab, resp in (("total", None), ("res1", 1), ("res2", 2), ("res3", 3))] \
    + [("Ar_coup_15_L_mDM_100__%s" % lab, .1, L15, resp) for lab, resp in (("total", None), ("res1", 1), ("res3", 3), ("res4_negative", 4))]


@pytest.mark.parametrize("name,m,coeffs,resp", SPECTRA_AR, ids=[s[0] for s in SPECTRA_AR])
def test_spectra_figs_5_7_argon(paper_halo_argon, name, m, coeffs, resp):
    '''Argon spectra of Figs. 5 (O_6) and 7 (O_7 contact at 10 MeV, O_15 long-range at 100 MeV, per
    term R_j W_j), rescaled by RHO_FIX: within 2% (within 1.1% in practice). O_6 contact at 1 GeV:
    median 0.993, within 5% below 500 eV and down to 0.65 at the end point (q_max = 1 MeV).'''
    E, y = spectrum(m, coeffs, resp, shells=None, target="Ar")
    if name.endswith("res4_negative"):
        assert np.all(y <= 0)
    r = ratio(name, E, y) / RHO_FIX
    if name == "Ar_coup_6_C__m1000MeV":
        x, yp = load_csv(name).T
        low = x[yp > 1e-3 * yp.max()] < 500.
        assert abs(np.median(r) - 0.993) < 0.005, np.median(r)
        assert np.all((r[low] > 0.95) & (r[low] < 1.0)) and np.all((r > 0.6) & (r < 1.0)), r
    else:
        assert np.all((r > 0.98) & (r < 1.02)), r


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
    '''Figs. 8 and 10, XENON10: within 4%'''
    m, r = limit_ratios(name, exp, setter)
    assert np.all(np.isfinite(r)), (m, r)
    assert np.all(np.abs(r - 1) < 0.04), (m, r)


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
