'''XENON1T S2-only DM-electron signal against the published limits of the data release
(tests/data/xenon1t_s2only/5d_results_dmelectron.csv, from
https://github.com/XENON1T/s2only_data_release/blob/master/limits/5d_results_dmelectron.csv,
arXiv:1907.11485): limit_cm2 is the excluded reference cross section of the heavy mediator dark photon
(F_DM = 1) and n_expected_signal_mean the expected signal in the ROI at that cross section.

XENON1T computes the signal following arXiv:1703.00910 with the ionization form factors of Essig/Yu
(ddlDM), RAPIDD uses the W_1 of arXiv:1912.08204. The S2 response is evaluated at E_e plus the
energy of the de-excitation photon of the shell (63.1, 87.9, 50.3 eV for 4d, 4p, 4s, 0 for the
others: "electron" convention, default), the conservative signal of XENON1T. With the same halo
(SHM, rho = 0.3 GeV/cm^3, v0 = 220, vesc = 544, vE = 232 km/s, assumed: the standard XENON1T halo):
- the spectra of wimprates (Essig form factors) folded in the same way give 0.91-1.01 times
  n_expected_signal_mean for 0.1 GeV <= m <= 1 GeV (0.76 at 50 MeV, 0.25 at 30 MeV, where the signal
  is close to the threshold of the response and the conversion of XENON1T, in produced electrons,
  is not reproduced by a shift of the energy; see examples/checks_ER_rates.ipynb): this is the
  XENON1T conversion above ~0.1 GeV;
- the RAPIDD signal at the published limit is ~0.38-0.51 of n_expected_signal_mean for
  50 MeV <= m <= 5 GeV (0.20 at 30 MeV): the ratio of the RAPIDD and Essig form factors.
With the total deposited energy E_e + E_b ("deposited" convention) the RAPIDD signal is
0.55-0.73 of n_expected_signal_mean (wimprates: 1.28-1.40, more than XENON1T).

Before the energy of the de-excitation was added to E_e, the ratio was ~0.07-0.25; before the
regeneration of the Xe 4d atomic responses and the Simpson q quadrature (deposited energy),
~0.45-0.56.
'''
import os

import numpy as np
import pytest

import rapidd

from conftest import DATA_DIR

X1T_RESULTS = os.path.join(DATA_DIR, "xenon1t_s2only", "5d_results_dmelectron.csv")
RHO = 0.3


@pytest.fixture(scope="module")
def xenon1t_halo(tmp_path_factory):
    path = str(tmp_path_factory.mktemp("halo") / "halo_table_SHM_XENON1T.dat")
    rapidd.calc_halo(path, vesc=544., v0=220., ve=232., i=0)
    rapidd.read_halo(path)
    rapidd.read_atomic_responses()
    yield
    rapidd.reset_coefficients()
    rapidd.read_halo(rapidd.halo_path)


def published():
    d = np.genfromtxt(X1T_RESULTS, delimiter=",", names=True)
    return d[(d['mass_gev'] >= 0.03) & (d['mass_gev'] <= 5.)]


@pytest.mark.parametrize("row", published(), ids=lambda r: "m%g" % r['mass_gev'])
def test_signal_at_published_limit(xenon1t_halo, row):
    m = float(row['mass_gev'])
    rapidd.reset_coefficients()
    rapidd.set_dark_photon_coeff(m, float(row['limit_cm2']))
    sig = rapidd.electron_signal("Xenon1T", rhoDM=RHO, mDM=m)
    assert sig['n_obs'] == int(row['n_observed'])
    np.testing.assert_allclose(sig['n_bkg'], row['n_expected_background_mean'], rtol=1e-6)
    ratio = sig['signal'] / row['n_expected_signal_mean']
    if m >= 0.05:
        assert 0.35 < ratio < 0.55, ratio
    else:
        assert 0.15 < ratio < 0.25, ratio


@pytest.mark.parametrize("row", published(), ids=lambda r: "m%g" % r['mass_gev'])
def test_signal_at_published_limit_deposited(xenon1t_halo, row):
    m = float(row['mass_gev'])
    rapidd.reset_coefficients()
    rapidd.set_dark_photon_coeff(m, float(row['limit_cm2']))
    rapidd.set_er_energy_convention("deposited")
    try:
        sig = rapidd.electron_signal("Xenon1T", rhoDM=RHO, mDM=m)
    finally:
        rapidd.set_er_energy_convention("electron")
    ratio = sig['signal'] / row['n_expected_signal_mean']
    assert 0.5 < ratio < 0.8, ratio
