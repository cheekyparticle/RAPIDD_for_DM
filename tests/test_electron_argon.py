'''DM-electron rate of Argon (target="Ar"): shells, binding energies and mass of the target,
independent evaluation of the rate of single shells, separation from the Xenon tables (the Xenon
default is unchanged). The comparison with the Argon figures of arXiv:1912.08204 is in
test_electron_paper.py (test_*_argon).'''
import contextlib

import numpy as np
import pytest

import rapidd

from test_electron_rate import independent_dRdlnE

HARTREE_EV = 27.211386245988
# -orbital energies [Hartree] of the DarkART data files data/Ar_<shell>.txt (Bunge et al. 1993)
AR_EB_HARTREE = {"1s": 118.610349, "2p": 9.571464, "2s": 12.322152, "3p": 0.591016, "3s": 1.277352}
M_AR_KG = 39.948 * 1.66053906660e-27


@pytest.fixture
def argon(paper_halo):
    rapidd.read_atomic_responses(target="Ar")
    yield


def test_target_data():
    assert rapidd.ELECTRON_TARGETS == ("Xe", "Ar")
    assert rapidd.AR_SHELLS == ("1s", "2p", "2s", "3p", "3s") == rapidd.AR_SHELLS_PAPER
    for s, shell in enumerate(rapidd.AR_SHELLS):
        np.testing.assert_allclose(rapidd.AR_BINDING_ENERGIES[s], AR_EB_HARTREE[shell] * HARTREE_EV * 1e-9, rtol=1e-14)
    # the Xenon binding energies of RAPIDD follow the same convention (1s: 1224.397767 Hartree)
    np.testing.assert_allclose(rapidd.XE_BINDING_ENERGIES[0], 1224.397767 * HARTREE_EV * 1e-9, rtol=1e-14)
    assert rapidd.ELECTRON_TARGET_MASS_AU == {"Xe": 131.293, "Ar": 39.948}
    assert rapidd.ELECTRON_TARGET_SHELLS["Xe"] == rapidd.XE_SHELLS


def test_tables_shape():
    '''Same grid as Xenon: 100 k' x 200 q after 3 comment lines'''
    import os
    for j in (1, 2, 3, 4):
        for shell in rapidd.AR_SHELLS:
            W = np.loadtxt(os.path.join(rapidd.core.atomic_resp_dir, "atomic_response_%d" % j, "Ar_%s.txt" % shell),
                           comments="//")
            assert W.shape == (rapidd.GRIDSIZE_K, rapidd.GRIDSIZE_Q) and np.all(np.isfinite(W))
            if j == 1:
                assert np.all(W >= 0)


@pytest.mark.parametrize("coeffs,mchi,shell", [({(7, "contact"): 1e-3}, 0.01, "3p"),
                                               ({(13, "long-range"): 1e-3}, 0.1, "3p"),
                                               ({(1, "contact"): 1e-3}, 0.1, "3s"),
                                               ({(15, "long-range"): 1e-3}, 0.1, "2p"),
                                               ({(6, "contact"): 1e-3}, 1.0, "2s")])
def test_independent_rate(argon, coeffs, mchi, shell):
    '''RAPIDD Argon shells vs an independent evaluation of Eqs. (ionization spectrum), (RW)'''
    for (op, inter), val in coeffs.items():
        rapidd.set_any_Ecoeff(val, op, interaction=inter)
    rapidd.electron._qmax_warned.clear()
    warns = mchi > 0.3 and any(op in rapidd.QMAX_OPERATORS and inter == "contact" for op, inter in coeffs)
    with pytest.warns(rapidd.QMaxTruncationWarning) if warns else contextlib.nullcontext():
        r = rapidd.difrate_electron(0.4, mchi, 0.5, shells=[shell], target="Ar")
    ik = [i for i in (5, 25, 45, 60, 75) if r['dRdlnE'][i] > 1e-3 * r['dRdlnE'].max()]
    assert len(ik) >= 2
    eb = rapidd.AR_BINDING_ENERGIES[rapidd.AR_SHELLS.index(shell)]
    ref = independent_dRdlnE(coeffs, mchi, 0.4, ik, table="Ar_" + shell, eb=eb, m_target_kg=M_AR_KG)
    np.testing.assert_allclose(r['dRdlnE'][ik], ref, rtol=5e-3)
    assert list(r['dRdE_shell']) == [shell] and r['target'] == "Ar"


def test_xenon_unchanged_by_argon(paper_halo):
    '''Reading the Argon tables does not change the Xenon rates (separate tables), target="Xe" is the
    default'''
    rapidd.set_any_Ecoeff(1e-3, 7)
    xe = rapidd.difrate_electron(0.4, 0.05)['dRdE']
    rapidd.read_atomic_responses(target="Ar")
    ar = rapidd.difrate_electron(0.4, 0.05, target="Ar")
    np.testing.assert_array_equal(rapidd.difrate_electron(0.4, 0.05, target="Xe")['dRdE'], xe)
    assert ar['rate_shell'].shape == (5, rapidd.GRIDSIZE_K)
    assert np.all(ar['dRdE'] >= 0) and ar['dRdE'].max() > 0
    np.testing.assert_allclose(sum(ar['dRdE_shell'].values()), ar['dRdE'], rtol=1e-12)
    # the totals of the two targets are comparable (Fig. 6 of arXiv:1912.08204)
    ratio = rapidd.total_rate_electron(0.4, 0.05, target="Ar") / rapidd.total_rate_electron(0.4, 0.05)
    assert 0.1 < ratio < 10


def test_argon_input_checks(paper_halo):
    rapidd.set_any_Ecoeff(1e-3, 1)
    with pytest.raises(ValueError):
        rapidd.difrate_electron(0.4, 0.1, target="Ge")
    with pytest.raises(ValueError):
        rapidd.read_atomic_responses(target="Ge")
    rapidd.read_atomic_responses(target="Ar")
    with pytest.raises(ValueError):
        rapidd.difrate_electron(0.4, 0.1, shells=["5p"], target="Ar")
    with pytest.raises(ValueError):
        rapidd.difrate_electron(0.4, 0.1, shells=[5], target="Ar")
    rate = rapidd.difrate_electron(0.4, 0.1, target="Ar")
    with pytest.raises(ValueError):
        rapidd.electron_signal("Xenon10", mDM=0.1, rate=rate)
