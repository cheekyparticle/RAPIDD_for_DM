'''Targets and files of the MadDM entry points DDrate_save and DDrate_electron_save, and the
DM-electron cross section without rates (electron_cross_section, MadDM direct_mode = limit).'''
import os
import sys

import numpy as np
import pytest

import rapidd
from rapidd.calc_dRdE import (DDrate_save, DDrate_electron_save, electron_cross_section, select_targets,
                              NUCLEAR_TARGET_NAMES)

from conftest import DATA_DIR

sys.path.insert(0, os.path.join(DATA_DIR, "regression_nuclear_modulation"))
from generate_reference import CARD, ALPHAQ  # noqa: E402
sys.path.pop(0)

M_E = 0.000510998950
ALPHAS_E = {"alpha_e": (1.0e-4, 2.0e-4, 1.5e-4, 0.5e-4)}


def test_select_targets():
    assert select_targets(["ge", "XE", "xe"], NUCLEAR_TARGET_NAMES) == ("Xe", "Ge")
    assert select_targets("ar", ("Xe", "Ar")) == ("Ar",)
    with pytest.raises(ValueError):
        select_targets(["ge"], ("Xe", "Ar"))
    with pytest.raises(ValueError):
        select_targets([], NUCLEAR_TARGET_NAMES)


def test_nuclear_targets(tmp_path):
    '''one column per target, in the order Xe, Ar, Ge, F, independent of the other targets'''
    (tmp_path / "all").mkdir()
    (tmp_path / "ge").mkdir()
    DDrate_save(CARD, ALPHAQ, 100., str(tmp_path / "all"), halo_path=rapidd.halo_path)
    DDrate_save(CARD, ALPHAQ, 100., str(tmp_path / "ge"), targets=["ge", "xe"], halo_path=rapidd.halo_path)
    assert sorted(os.listdir(tmp_path / "ge")) == ["dRdE_nucleon.dat"]
    full = np.loadtxt(tmp_path / "all" / "dRdE_nucleon.dat")
    sel = np.loadtxt(tmp_path / "ge" / "dRdE_nucleon.dat")
    assert full.shape == (1000, 5) and sel.shape == (1000, 3)
    np.testing.assert_array_equal(sel, full[:, [0, 1, 3]])
    with open(tmp_path / "ge" / "dRdE_nucleon.dat") as f:
        assert f.readline().split()[1:6:2] == ["E_R", "drde_Xe", "drde_Ge"]
    with open(tmp_path / "all" / "dRdE_nucleon.dat") as f:
        assert f.readline().split()[1:10:2] == ["E_R", "drde_Xe", "drde_Ar", "drde_Ge", "drde_F"]
    assert np.all(full[:, 4] > 0)


def test_electron_targets(tmp_path):
    '''Argon only: rates but no XENON10 / XENON1T signals; the Xenon columns do not depend on Argon'''
    (tmp_path / "xe").mkdir()
    (tmp_path / "ar").mkdir()
    (tmp_path / "xear").mkdir()
    res_xe = DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path / "xe"), halo_path=rapidd.halo_path)
    res_ar = DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path / "ar"), targets=["ar"],
                                  halo_path=rapidd.halo_path)
    DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path / "xear"), targets=["ar", "xe"],
                         halo_path=rapidd.halo_path)
    assert sorted(os.listdir(tmp_path / "xe")) == ["dRdE_electron.dat", "dRdS2_electron_XENON10.dat",
                                                   "dRdS2_electron_XENON1T.dat", "signal_e_recoil.dat"]
    assert sorted(os.listdir(tmp_path / "ar")) == ["dRdE_electron.dat"]
    assert "Xenon10_signal" in res_xe and "Xenon10_signal" not in res_ar and "Xenon1T_signal" not in res_ar
    assert res_ar["sigma_e"] == res_xe["sigma_e"]
    xe = np.loadtxt(tmp_path / "xe" / "dRdE_electron.dat")
    ar = np.loadtxt(tmp_path / "ar" / "dRdE_electron.dat")
    both = np.loadtxt(tmp_path / "xear" / "dRdE_electron.dat")
    n_xe, n_ar = len(rapidd.XE_SHELLS), len(rapidd.AR_SHELLS)
    assert xe.shape[1] == 3 + n_xe and ar.shape[1] == 3 + n_ar and both.shape[1] == 4 + n_xe + n_ar
    np.testing.assert_array_equal(both[:, [0, 1, 2]], xe[:, [0, 1, 2]])
    np.testing.assert_array_equal(both[:, 3], ar[:, 2])
    np.testing.assert_array_equal(both[:, 4:4 + n_xe], xe[:, 3:])
    np.testing.assert_array_equal(both[:, 4 + n_xe:], ar[:, 3:])
    # total = sum of the shells (up to the 7 digits of the text file)
    np.testing.assert_allclose(both[:, 2], both[:, 4:4 + n_xe].sum(axis=1), rtol=1e-5)
    np.testing.assert_allclose(both[:, 3], both[:, 4 + n_xe:].sum(axis=1), rtol=1e-5)


def test_electron_silicon(tmp_path, monkeypatch):
    '''Silicon: dRdE_electron_Si.dat and dRdQ_electron_Si.dat on their own grids; the Xenon files do not change.
    With a coupling that needs W_2 ... W_5 Si is skipped with a warning, the other targets are computed.'''
    from rapidd import calc_dRdE
    from rapidd.electron import difrate_electron, difrate_dQ
    for d in ("xe", "si", "xesi", "skip"):
        (tmp_path / d).mkdir()
    DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path / "xe"), halo_path=rapidd.halo_path)
    res = DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path / "si"), targets=["si"],
                               halo_path=rapidd.halo_path)
    DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path / "xesi"), targets=["SI", "xe"],
                         halo_path=rapidd.halo_path)
    assert sorted(os.listdir(tmp_path / "si")) == ["dRdE_electron_Si.dat", "dRdQ_electron_Si.dat"]
    assert sorted(os.listdir(tmp_path / "xesi")) == sorted(os.listdir(tmp_path / "xe") +
                                                           ["dRdE_electron_Si.dat", "dRdQ_electron_Si.dat"])
    np.testing.assert_array_equal(np.loadtxt(tmp_path / "xesi" / "dRdE_electron.dat"),
                                  np.loadtxt(tmp_path / "xe" / "dRdE_electron.dat"))
    E, drde = np.loadtxt(tmp_path / "si" / "dRdE_electron_Si.dat", unpack=True)
    Q, drdq = np.loadtxt(tmp_path / "si" / "dRdQ_electron_Si.dat", unpack=True)
    assert E.size == rapidd.CRYSTAL_NE - 1 and E[0] == pytest.approx(1e-4) and E[-1] == pytest.approx(0.05)
    assert list(Q) == list(range(1, calc_dRdE.CRYSTAL_Q_MAX + 1)) and np.all(drdq[:3] > 0)
    np.testing.assert_allclose(drdq, res["dRdQ_Si"], rtol=1e-6)
    # same rate as rapidd.electron with the coefficients of the MadDM alphas
    calc_dRdE.electron_cross_section(ALPHAS_E, 2, M_E, 0.1)
    r = difrate_electron(CARD['rhoDM'], 0.1, 0.5, m_e=M_E, target="Si")
    np.testing.assert_allclose(drde, r['dRdE'], rtol=1e-6)
    np.testing.assert_allclose(drdq, difrate_dQ(r)['dRdQ'], rtol=1e-6)
    assert res["skipped_targets"] == {}
    # unsupported coupling (not generated by the MadDM alphas, c_1 and c_4 only): Si skipped, Xenon computed
    orig = calc_dRdE.electron_cross_section

    def with_c3(*args):
        out = orig(*args)
        rapidd.set_any_Ecoeff(1e-3, 3)
        return out
    monkeypatch.setattr(calc_dRdE, "electron_cross_section", with_c3)
    with pytest.warns(calc_dRdE.CrystalTargetSkippedWarning, match="c_3"):
        res = DDrate_electron_save(CARD, ALPHAS_E, 2, M_E, 0.1, str(tmp_path / "skip"), targets=["xe", "si"],
                                   halo_path=rapidd.halo_path)
    assert list(res["skipped_targets"]) == ["Si"] and "dRdQ_Si" not in res
    assert "dRdE_electron_Si.dat" not in os.listdir(tmp_path / "skip")
    assert "dRdE_electron.dat" in os.listdir(tmp_path / "skip")
    rapidd.reset_coefficients()


@pytest.mark.parametrize("dm_spin,m", ((2, 0.1), (1, 0.05), (2, 3.0)))
def test_electron_cross_section(tmp_path, dm_spin, m):
    '''same DM response and sigma_e as DDrate_electron_save, without the rates'''
    res = DDrate_electron_save(CARD, ALPHAS_E, dm_spin, M_E, m, str(tmp_path), halo_path=rapidd.halo_path)
    xs = electron_cross_section(ALPHAS_E, dm_spin, M_E, m)
    assert xs == {"DM_response": res["DM_response"], "sigma_e": res["sigma_e"]}
    assert electron_cross_section(ALPHAS_E, 3, M_E, m) == {"DM_response": -1}
