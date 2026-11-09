'''XENON1T S2-only DM-electron signal against the published limits of the data release
(tests/data/xenon1t_s2only/5d_results_dmelectron.csv, from
https://github.com/XENON1T/s2only_data_release/blob/master/limits/5d_results_dmelectron.csv,
arXiv:1907.11485): limit_cm2 is the excluded reference cross section of the heavy mediator dark photon
(F_DM = 1) and n_expected_signal_mean the expected signal in the ROI at that cross section.

XENON1T computes the signal following arXiv:1703.00910 with the ionization form factors of Essig/Yu
(ddlDM), RAPIDD uses the W_1 of arXiv:1912.08204. With the same halo (SHM, rho = 0.3 GeV/cm^3,
v0 = 220, vesc = 544, vE = 232 km/s, assumed: the standard XENON1T halo) the RAPIDD signal at the
published limit is ~0.55-0.73 of n_expected_signal_mean for 30 MeV <= m <= 5 GeV (0.24 at 20 MeV,
where the signal comes from the tail of the halo). The ratio is explained by the form factors:
folding the spectra of wimprates (Essig form factors) with the same halo and the same response gives
1.8-2.2 times the RAPIDD signal (see examples/checks_ER_rates.ipynb).

Before the deposited energy E_e + E_b was used in the response, the ratio was ~0.07-0.25; before the
regeneration of the Xe 4d atomic responses and the Simpson q quadrature, ~0.45-0.56.
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
    assert 0.5 < ratio < 0.8, ratio
