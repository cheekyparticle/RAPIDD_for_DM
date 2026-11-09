"""Reference outputs of the DM-electron scattering (SI/SD only: Ce(1), Ce(4), atomic response W_1),
used by tests/test_electron_regression.py. They were first produced with the branch
direct_electron_port and regenerated on 2026-11-09 with the current code, after the fixes listed in
tests/test_electron_regression.py (XENON1T deposited energy, q quadrature, Xe_4d tables, XENON10
E integral).

Run from the repository root with the library built:
    python tests/data/regression_direct_electron_port/generate_reference.py
It writes reference.npz and the DDrate_electron_save outputs (ddrate_<case>/) next to this
script. Regenerate it only after a documented change of the results. Options:
    --ddrate-only      only the DDrate_electron_save outputs
    --modulation-only  only the annual-modulation case ddrate_spin2_m0.1_mod (day_bins time bins
                       centered at t_i = (i + 0.5) 365.24/day_bins, omega = 2 pi/365.24, vEarth_mod
                       from the card)
"""
import sys
import os
import shutil
import numpy as np
import rapidd
from rapidd.calc_dRdE import DDrate_electron_save

HERE = os.path.dirname(os.path.abspath(__file__))
M_E = 0.000510998950  # GeV
RHO = 0.3             # GeV/cm^3
MASSES = (0.05, 0.1, 1.0)
JCHIS = (0.5, 0.0)
# (label, Ce(1), Ce(4))
COEFFS = (("c1", 1.0e-3, 0.0), ("c4", 0.0, 1.0e-3), ("c1c4", 1.0e-3, 2.0e-3))

# card and alphas for DDrate_electron_save (alpha_e = (aSIe, aSIo, aSDe, aSDo) in GeV^-2)
CARD = {'rhoDM': RHO, 'vmp': 238., 'vescape': 544., 'vSun_r': 11.1, 'vSun_phi': 12.2, 'vSun_theta': 7.3}
ALPHAS = {"alpha_e": (1.0e-4, 2.0e-4, 1.5e-4, 0.5e-4)}
DD_CASES = ((2, 0.1, False), (2, 1.0, False), (1, 0.05, False), (2, 0.1, True))  # (MadDM dm_spin, mDM, modulation)
DAY_BINS = 4  # time bins of the modulation case


def case_key(label, m, jchi):
    return "%s_m%g_j%g" % (label, m, jchi)


def main(modulation_only=False, ddrate_only=False):
    if modulation_only or ddrate_only:
        write_ddrate_cases(modulation_only=modulation_only)
        return
    ref = {}
    rapidd.read_halo(rapidd.halo_path)
    rapidd.read_atomic_responses(response_index=1)
    for label, c1, c4 in COEFFS:
        for m in MASSES:
            for jchi in JCHIS:
                rapidd.reset_coefficients()
                rapidd.set_any_Ecoeff(c1, 1)
                rapidd.set_any_Ecoeff(c4, 4)
                key = case_key(label, m, jchi)
                k, E, rate_shell = rapidd.difrate_electron_dEr_all(RHO, m, jchi, M_E)
                x10 = rapidd.difrate_dS2("Xenon10", E, rate_shell, m)
                x1t = rapidd.difrate_dS2("Xenon1T", E, rate_shell, m)
                ref[key + "/k"] = k
                ref[key + "/E"] = E
                ref[key + "/rate_shell"] = rate_shell
                ref[key + "/dm_response"] = np.array(rapidd.dm_response_electron(m, jchi))
                ref[key + "/X10_dRdS2"] = x10['dRdS2']
                ref[key + "/X10_dRdS2_shell"] = x10['dRdS2_shell']
                ref[key + "/X10_signal"] = x10['signal']
                ref[key + "/X1T_dRdS2"] = x1t['dRdS2']
                ref[key + "/X1T_signal"] = np.array(x1t['signal'])
                ref[key + "/X1T_bkg"] = np.array(x1t['n_bkg'])
                ref[key + "/X1T_obs"] = np.array(x1t['n_obs'])
    np.savez(os.path.join(HERE, "reference.npz"), **ref)
    write_ddrate_cases()


def write_ddrate_cases(modulation_only=False):
    for dm_spin, m, modulation in DD_CASES:
        if modulation_only and not modulation:
            continue
        out = os.path.join(HERE, "ddrate_spin%i_m%g%s" % (dm_spin, m, "_mod" if modulation else ""))
        shutil.rmtree(out, ignore_errors=True)
        os.makedirs(out)
        res = DDrate_electron_save(CARD, ALPHAS, dm_spin, M_E, m, out, halo_path=rapidd.halo_path,
                                   modulation=modulation, day_bins=DAY_BINS)
        shutil.rmtree(os.path.join(out, "halo_table_modulation"), ignore_errors=True)
        with open(os.path.join(out, "result.txt"), "w") as f:
            for key in sorted(res):
                f.write("%s %r\n" % (key, res[key]))


if __name__ == "__main__":
    main(modulation_only="--modulation-only" in sys.argv[1:], ddrate_only="--ddrate-only" in sys.argv[1:])
