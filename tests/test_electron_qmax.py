'''QMaxTruncationWarning: the q integral of the DM-electron rate stops at q = 1 MeV ~ 2 m_e (end of the
atomic response tables, as in DarkARC/DarkART and arXiv:1912.08204). For mDM > 0.3 GeV the contact
rates of O_3, O_5, O_6, O_13, O_14, O_15 get a sizeable contribution from larger q, where the NREFT and
the non-relativistic wave functions are not valid: RAPIDD warns once per combination (target,
operators, mass bin). O_1 and O_4 (the only operators used by MadDM) never trigger it.'''
import warnings

import numpy as np
import pytest

import rapidd
import rapidd.electron as electron

from test_electron_regression import CARD, ALPHAS, M_E

RHO = 0.4
OUTER_SHELL = {"Xe": "5p", "Ar": "3p"}


@pytest.fixture(autouse=True)
def fresh_warnings(default_halo):
    '''Each test starts with no combination already warned'''
    electron._qmax_warned.clear()
    yield
    electron._qmax_warned.clear()


def qmax_warnings(record):
    return [w for w in record if issubclass(w.category, rapidd.QMaxTruncationWarning)]


def run(m, coeffs, n_calls=1, target="Xe"):
    rapidd.reset_coefficients()
    for (op, inter), val in coeffs.items():
        rapidd.set_any_Ecoeff(val, op, interaction=inter)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        for _ in range(n_calls):
            rapidd.total_rate_electron(RHO, m, 0.5, shells=[OUTER_SHELL[target]], target=target)
    return qmax_warnings(rec)


@pytest.mark.parametrize("op", rapidd.QMAX_OPERATORS)
def test_warning_contact(op):
    '''Contact O_3, O_5, O_6, O_13, O_14, O_15 above 0.3 GeV: one warning with the lost fraction'''
    w = run(0.5, {(op, "contact"): 1.}, n_calls=3)
    assert len(w) == 1
    msg = str(w[0].message)
    assert "q = 1 MeV ~ 2 m_e" in msg and "NREFT" in msg and "O_%d ~" % op in msg
    f, lower = rapidd.qmax_lost_fraction(op, 0.5)
    assert not lower and ("%.2g%%" % (100 * f)) in msg


def test_warning_once_per_combination():
    '''Repeated calls do not warn again; a new mass bin, set of operators or target does'''
    rapidd.read_atomic_responses(target="Ar")
    assert len(run(0.6, {(6, "contact"): 1.})) == 1
    assert len(run(0.65, {(6, "contact"): 1.})) == 0        # same bin (0.5, 0.7]
    assert len(run(0.9, {(6, "contact"): 1.})) == 1         # bin (0.7, 1]
    assert len(run(0.9, {(6, "contact"): 1., (15, "contact"): 2.})) == 1
    assert len(run(0.9, {(6, "contact"): 1., (1, "contact"): 2.})) == 0  # O_1 does not change the set
    assert len(run(0.9, {(6, "contact"): 1.}, target="Ar")) == 1
    w = run(3., {(15, "contact"): 1.})
    assert len(w) == 1 and "O_15 > 38%" in str(w[0].message)  # above 1 GeV: lower bound


@pytest.mark.parametrize("m,coeffs", [
    (0.3, {(15, "contact"): 1.}),                           # threshold: mDM > 0.3 GeV
    (0.1, {(6, "contact"): 1., (13, "contact"): 1.}),
    (2.0, {(1, "contact"): 1., (4, "contact"): 1.}),        # MadDM operators
    (2.0, {(op, "long-range"): 1. for op in rapidd.QMAX_OPERATORS}),
    (2.0, {(op, "contact"): 1. for op in (7, 8, 9, 10, 11, 12)}),
])
def test_no_warning(m, coeffs):
    assert run(m, coeffs, n_calls=2) == []


@pytest.mark.parametrize("model", rapidd.ELECTRON_MODELS)
def test_models_no_warning(model):
    '''The models of arXiv:1912.08204 (set_any_E_EFT_coeff) have no contact O_3, O_5, O_6, O_13-O_15'''
    rapidd.reset_coefficients()
    rapidd.set_any_E_EFT_coeff(model, 1., 1e-3)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        rapidd.difrate_electron(RHO, 1., 0.5, shells=["5p"])
    assert qmax_warnings(rec) == []


def test_set_any_E_EFT_coeff_contact_operator(monkeypatch):
    '''Any contact coefficient set through set_any_E_EFT_coeff is checked as well'''
    monkeypatch.setattr(electron, "electron_model_coeffs",
                        lambda model, mDM, coupling, **kw: {(13, "contact"): coupling})
    rapidd.reset_coefficients()
    electron.set_any_E_EFT_coeff("anapole", 1., 1e-3)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        rapidd.difrate_electron(RHO, 1., 0.5, shells=["5p"])
    assert len(qmax_warnings(rec)) == 1


@pytest.mark.parametrize("dm_spin", (1, 2))
def test_maddm_entry_point_no_warning(tmp_path, dm_spin):
    '''DDrate_electron_save (MadDM: O_1, O_4) never warns, even at large mass'''
    from rapidd.calc_dRdE import DDrate_electron_save
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        DDrate_electron_save(CARD, ALPHAS, dm_spin, M_E, 3.0, str(tmp_path), halo_path=rapidd.halo_path)
    assert qmax_warnings(rec) == []


def test_lost_fraction_table():
    '''Interpolation of QMAX_LOST_FRACTION: increasing with the mass, exact at the nodes'''
    for op in rapidd.QMAX_OPERATORS:
        ms = np.linspace(0.31, 1.0, 50)
        f = [rapidd.qmax_lost_fraction(op, m)[0] for m in ms]
        assert np.all(np.diff(f) > 0)
        for m, row in rapidd.QMAX_LOST_FRACTION.items():
            assert rapidd.qmax_lost_fraction(op, m) == (row[op], False)
    assert rapidd.qmax_lost_fraction(1, 2.) == (0., False)


def test_q_max_is_end_of_tables():
    '''Q_MAX_GEV is the last point of the q grid of the atomic response tables (q_max = 1 MeV)'''
    assert rapidd.Q_MAX_GEV == 1e-3
    with open(rapidd.core.os.path.join(rapidd.core.atomic_resp_dir, "README.md")) as f:
        assert "q_max = 1000.0" in f.read()
