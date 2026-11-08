'''Annual modulation (rapidd.modulation), shared by DM-nucleus and DM-electron scattering:
time bins, Earth velocity, halo tables (vEarth_mod), common interface and physics checks.'''
import inspect
import os
import sys

import numpy as np
import pytest

import rapidd
from rapidd import modulation
from rapidd.calc_dRdE import DDrate_save, DDrate_electron_save
from rapidd.modulation import modulation_times, modulation_halo_from_card, trapezoid_weights

from conftest import DATA_DIR

sys.path.insert(0, os.path.join(DATA_DIR, "regression_nuclear_modulation"))
from generate_reference import CARD, ALPHAQ  # noqa: E402
sys.path.pop(0)

M_E = 0.000510998950
ALPHAS_E = {"alpha_e": (1.0e-4, 2.0e-4, 1.5e-4, 0.5e-4)}
MODULATION_KWARGS = ("halo_path", "modulation", "modulation_halo_dir", "day_bins", "verbose", "save_rates_per_bin")


def test_modulation_times():
    t, dt = modulation_times(4)
    assert dt == pytest.approx(365.24 / 4, rel=1e-15)
    np.testing.assert_allclose(t, (np.arange(4) + 0.5) * dt, rtol=1e-15)
    assert rapidd.YEAR_DAYS == 365.24
    with pytest.raises(ValueError):
        modulation_times(0)


def test_earth_velocity_period_and_phase():
    '''omega = 2 pi/365.24: |v_lab| has period 365.24 days and is maximum at the beginning of June
    (~72 days after March 22); the annual-average speed (March 9 vector) is between min and max'''
    np.testing.assert_allclose(rapidd.lab_frame_speed([0., 100.3]), rapidd.lab_frame_speed([365.24, 465.54]),
                               rtol=1e-13)
    t = np.linspace(0., 365.24, 36525)
    v = rapidd.lab_frame_speed(t)
    assert 70. < t[v.argmax()] < 76.
    assert v.min() < rapidd.lab_frame_speed_annual_avg() < v.max()
    # v_earth scales the Earth velocity only
    v0 = rapidd.lab_frame_speed(t, v_earth=0.)
    np.testing.assert_allclose(v0, v0[0], rtol=1e-14)


def test_same_interface():
    '''DM-nucleus and DM-electron: same modulation arguments, defaults and order'''
    def tail(f):
        p = list(inspect.signature(f).parameters.values())
        i = [q.name for q in p].index("halo_path")
        return [(q.name, q.default) for q in p[i:]]
    assert tail(DDrate_save) == tail(DDrate_electron_save)
    assert [n for n, _ in tail(DDrate_save)] == list(MODULATION_KWARGS)
    assert dict(tail(DDrate_save))["day_bins"] == 12
    assert dict(tail(DDrate_save))["verbose"] is False
    assert dict(tail(DDrate_save))["save_rates_per_bin"] is True


def _count_calc_halo(monkeypatch):
    calls = []
    calc = rapidd.core.calc_halo

    def counting(*args, **kwargs):
        calls.append(kwargs)
        return calc(*args, **kwargs)
    monkeypatch.setattr(rapidd.core, "calc_halo", counting)
    monkeypatch.setattr(modulation, "calc_halo", counting)
    return calls


def test_halo_path_from_card_ignores_vEarth_mod(tmp_path, monkeypatch):
    '''The annual-average table does not depend on vEarth_mod (Eq. 12): no table is computed
    for the default card nor for a different vEarth_mod; it is computed if vescape changes'''
    calls = _count_calc_halo(monkeypatch)
    card = dict(CARD)
    assert rapidd.halo_path_from_card(card, str(tmp_path)) == rapidd.halo_path
    card['vEarth_mod'] = 25.
    assert rapidd.halo_path_from_card(card, str(tmp_path)) == rapidd.halo_path
    assert calls == [] and os.listdir(tmp_path) == []
    card['vescape'] = 600.
    assert rapidd.halo_path_from_card(card, str(tmp_path)) == str(tmp_path / "SHM.dat")
    assert len(calls) == 1 and os.path.isfile(tmp_path / "SHM.dat")


def test_modulation_tables_vEarth_mod(tmp_path, monkeypatch):
    '''The modulation tables use vEarth_mod (default 29.8 if missing from the card); they are
    reused if already computed for the same parameters and recomputed if vEarth_mod changes'''
    calls = _count_calc_halo(monkeypatch)
    d = str(tmp_path / "mod")
    card = {k: v for k, v in CARD.items() if k != 'vEarth_mod'}
    paths = modulation_halo_from_card(card, d, day_bins=3)
    assert len(calls) == 1 and calls[0]['v_earth'] == rapidd.default_halo_v_earth
    t_default = [open(p).read() for p in paths]
    # same parameters (vEarth_mod = default given explicitly): not recomputed
    modulation_halo_from_card(dict(card, vEarth_mod=29.8), d, day_bins=3)
    assert len(calls) == 1
    # different vEarth_mod: recomputed, different tables
    paths = modulation_halo_from_card(dict(card, vEarth_mod=15.), d, day_bins=3)
    assert len(calls) == 2 and calls[1]['v_earth'] == 15.
    assert all(open(p).read() != t for p, t in zip(paths, t_default))
    # different day_bins: recomputed
    assert len(modulation_halo_from_card(dict(card, vEarth_mod=15.), d, day_bins=2)) == 2
    assert len(calls) == 3


def test_verbose(tmp_path, capfd):
    '''verbose=False silences the C loop over the time bins'''
    t, _ = modulation_times(2)
    rapidd.calc_halo(None, times=t, table_dir=str(tmp_path / "a"), verbose=False)
    assert capfd.readouterr().out == ""
    rapidd.calc_halo(None, times=t, table_dir=str(tmp_path / "b"), verbose=True)
    out = capfd.readouterr().out
    assert "Calculating for t=91.31 days" in out and "Done!" in out


def _check_common_outputs(npz, names, day_bins, E, weights, static):
    t, dt = modulation_times(day_bins)
    np.testing.assert_allclose(npz['t_bin'], t, rtol=1e-15)
    assert float(npz['dt_bin']) == pytest.approx(dt, rel=1e-15)
    for n in names:
        b = npz['drde_%s_bin' % n]
        assert b.shape == (day_bins, len(E))
        np.testing.assert_allclose(npz['rate_%s_bin' % n], dt * (b @ weights), rtol=1e-12)
        np.testing.assert_allclose(npz['drde_%s_year' % n], b.sum(axis=0) * dt, rtol=1e-12)
        np.testing.assert_allclose(npz['drde_%s_avg' % n], npz['drde_%s_year' % n] / 365.24, rtol=1e-12)
        # time average close to the spectrum with the annual-average halo
        tot_avg, tot_static = npz['drde_%s_avg' % n] @ weights, npz['drde_%s' % n] @ weights
        assert abs(tot_avg / tot_static - 1) < 0.01


def test_nuclear_modulation_physics(tmp_path):
    '''Heavy elastic DM in Xe: above the phase-reversal energy (~16 keV for 100 GeV) the rate is
    maximum at the beginning of June with an amplitude of ~1%; the rate integrated from 0.01 keV
    is dominated by low energies (small vmin) and has the opposite phase'''
    day_bins = 12
    DDrate_save(CARD, ALPHAQ, 100., str(tmp_path), modulation=True, day_bins=day_bins)
    npz = np.load(tmp_path / "DDrates.npz")
    E = npz['energies']
    _check_common_outputs(npz, ("Xe", "Ar", "Ge"), day_bins, E, trapezoid_weights(E), None)
    t, dt = modulation_times(day_bins)
    w = 2 * np.pi / 365.24
    M = np.c_[np.ones_like(t), np.cos(w * t), np.sin(w * t)]

    def fit(r):
        c = np.linalg.lstsq(M, r, rcond=None)[0]
        return np.hypot(c[1], c[2]) / c[0], (np.arctan2(c[2], c[1]) / w) % 365.24
    sel = E >= 25.
    amp, phase = fit(npz['drde_Xe_bin'][:, sel] @ trapezoid_weights(E[sel]))
    assert 70. < phase < 76. and 0.002 < amp < 0.05
    amp, phase = fit(npz['rate_Xe_bin'])
    assert 70. + 365.24 / 2 < phase < 76. + 365.24 / 2 and amp < 0.05
    # files
    rt = np.loadtxt(tmp_path / "rate_vs_time.dat")
    np.testing.assert_allclose(rt[:, 0], t, rtol=1e-10)
    np.testing.assert_allclose(rt[:, 2:], np.c_[npz['rate_Xe_bin'], npz['rate_Ar_bin'], npz['rate_Ge_bin']],
                               rtol=1e-10)
    sp = np.loadtxt(tmp_path / "dRdE_modulation.dat")
    np.testing.assert_allclose(sp[:, 1], npz['drde_Xe_avg'], rtol=1e-10)


def test_electron_modulation(tmp_path):
    day_bins = 6
    DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path), modulation=True, day_bins=day_bins)
    npz = np.load(tmp_path / "DDrates_e.npz")
    E = npz['E_keV']
    k, Eg, _ = rapidd.difrate_electron_dEr_all(0.3, 0.1, 0.5, M_E)
    _check_common_outputs(npz, ("tot",), day_bins, E, np.diff(Eg) * 1e6, None)
    rt = np.loadtxt(tmp_path / "rate_vs_time_e_recoil.dat")
    np.testing.assert_allclose(rt[:, 2], npz['rate_tot_bin'], rtol=1e-10)
    assert os.path.isfile(tmp_path / "dRdE_modulation_e_recoil.dat")


@pytest.mark.parametrize("target", ("nucleus", "electron"))
def test_save_rates_per_bin_false(tmp_path, target):
    if target == "nucleus":
        DDrate_save(CARD, ALPHAQ, 100., str(tmp_path), modulation=True, day_bins=2, save_rates_per_bin=False)
        npz, names = np.load(tmp_path / "DDrates.npz"), ("Xe", "Ar", "Ge")
    else:
        DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path), modulation=True, day_bins=2,
                             save_rates_per_bin=False)
        npz, names = np.load(tmp_path / "DDrates_e.npz"), ("tot",)
    for n in names:
        assert 'drde_%s_bin' % n not in npz.files
        for key in ('rate_%s_bin', 'drde_%s_year', 'drde_%s_avg'):
            assert key % n in npz.files
    assert 't_bin' in npz.files and 'dt_bin' in npz.files
