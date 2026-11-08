"""Reference outputs of DDrate_save (DM-nucleus spectra in Xe, Ar, Ge for the MadDM couplings
c1, c4), without and with the annual modulation (rapidd.modulation, DAY_BINS time bins), used by
tests/test_nuclear_regression.py. Same structure as regression_direct_electron_port.

Run from the repository root:
    python tests/data/regression_nuclear_modulation/generate_reference.py
It writes the DDrate_save outputs (ddrate_<case>/) next to this script.
"""
import os
import shutil

import rapidd
from rapidd.calc_dRdE import DDrate_save

HERE = os.path.dirname(os.path.abspath(__file__))
# maddm_card: halo and the default quark -> nucleon form factors of MadDM
CARD = {'rhoDM': 0.3, 'vmp': 238., 'vescape': 544., 'vSun_r': 11.1, 'vSun_phi': 12.2, 'vSun_theta': 7.3,
        'vEarth_mod': 29.8,
        'SPu': 0.0153, 'SPd': 0.0191, 'SPs': 0.0447, 'SPg': 1 - 0.0153 - 0.0191 - 0.0447,
        'SNu': 0.0110, 'SNd': 0.0273, 'SNs': 0.0447, 'SNg': 1 - 0.0110 - 0.0273 - 0.0447,
        'VPu': 2.0, 'VPd': 1.0, 'VNu': 1.0, 'VNd': 2.0,
        'AVPu': 0.842, 'AVPd': -0.427, 'AVPs': -0.085, 'AVNu': -0.427, 'AVNd': 0.842, 'AVNs': -0.085,
        'SigPu': 0.84, 'SigPd': -0.23, 'SigPs': -0.046, 'SigNu': -0.23, 'SigNd': 0.84, 'SigNs': -0.046}
# DM-quark alphas (SIe, SIo, SDe, SDo) [GeV^-2]
ALPHAQ = {'alpha_u': (1.0e-6, 0.5e-6, 1.0e-6, 0.), 'alpha_d': (1.0e-6, 0., 2.0e-6, 0.)}
DD_CASES = ((100., False), (100., True), (10., True))  # (mDM, modulation)
DAY_BINS = 4  # time bins of the modulation cases


def case_name(m, modulation):
    return "ddrate_m%g%s" % (m, "_mod" if modulation else "")


def main():
    for m, modulation in DD_CASES:
        out = os.path.join(HERE, case_name(m, modulation))
        shutil.rmtree(out, ignore_errors=True)
        os.makedirs(out)
        DDrate_save(CARD, ALPHAQ, m, out, halo_path=rapidd.halo_path, modulation=modulation, day_bins=DAY_BINS)
        shutil.rmtree(os.path.join(out, "halo_table_modulation"), ignore_errors=True)


if __name__ == "__main__":
    main()
