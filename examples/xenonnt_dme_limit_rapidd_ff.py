"""XENONnT S2-only DM-electron limit (arXiv:2601.11296, heavy mediator, F_DM = 1) rescaled to the
RAPIDD atomic responses (W_1 of arXiv:1912.08204, DarkART).

The published limit (black curves of Fig. 4(d), lib/experiments/XENONnT/XENONnT_DMe_SI.dat) is
computed with wimprates and the ionization form factors of Essig, Volansky, Yu (arXiv:1703.00910,
ddlDM tables). The data release has no response matrix for electronic recoils, so the limit cannot
be recomputed with the RAPIDD spectrum. It is rescaled by the ratio of the signals,

    sigma_e^RAPIDD(m) = sigma_e^XENONnT(m) * S_wimprates(m) / S_RAPIDD(m),

with S the rate integrated over the electronic recoil energy E_ER = E_e + E_sec in the region of
interest [0.04, 0.7] keV of the paper, E_sec the energy of the de-excitation photon of the shell
("electron" convention of XENONnT and XENON1T, rapidd.XE_DEEXCITATION_ENERGIES), for the 5 outer
shells and the same halo for the two spectra (the RAPIDD default table, SHM of arXiv:2105.00599, as
XENONnT: wimprates is given the same eta(v_min)). The rescaling assumes the same spectral shape
within the region of interest: the efficiency and the charge yield of XENONnT, which are not
released, are not included (with [0.1, 0.7] keV the ratio decreases by up to 20%, at 50-100 MeV).

Requires wimprates (pip install wimprates); run from the repository root with the library built:
    python examples/xenonnt_dme_limit_rapidd_ff.py
It writes lib/experiments/XENONnT/XENONnT_DMe_SI_RAPIDD_FF.dat.
"""
import os
import sys
import warnings

import numpy as np
from scipy.interpolate import interp1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import rapidd  # noqa: E402

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import wimprates as wr
    import numericalunits as nu

E_ROI = (0.04, 0.7)  # keV, arXiv:2601.11296
SHELLS = ('4s', '4p', '4d', '5s', '5p')
LIMIT_IN = os.path.join(rapidd.core.experiments_dir, 'XENONnT', 'XENONnT_DMe_SI.dat')
LIMIT_OUT = os.path.join(rapidd.core.experiments_dir, 'XENONnT', 'XENONnT_DMe_SI_RAPIDD_FF.dat')

# with recent scipy dme_ionization_ff returns an array of length 1 to quad: return a scalar
_ff = wr.electron.dme_ionization_ff
wr.electron.dme_ionization_ff = lambda shell, e, q: float(np.squeeze(_ff(shell, e, q)))
WR_UNIT = 1. / (nu.kg * nu.day * nu.keV)


def use_halo(path):
    '''same eta(v_min) in RAPIDD and wimprates'''
    rapidd.read_halo(path)
    h = np.loadtxt(path, skiprows=1)  # v_min [km/s], eta [(km/s)^-1]
    wr.electron.inverse_mean_speed_kms = interp1d(h[:, 0], h[:, 1], bounds_error=False, fill_value=0.)


def spectra(m, sigma_e=1e-40, rho=0.3, E_wr=np.geomspace(1e-3, 1., 80)):
    '''per-shell dR/dE_e [events/kg/day/keV] of RAPIDD and wimprates (dark photon, F_DM = 1)'''
    rapidd.reset_coefficients()
    rapidd.set_dark_photon_coeff(m, sigma_e)
    r = rapidd.difrate_electron(rho, m)
    E = r['E_GeV'][:-1] * 1e6
    out = {}
    for sh in SHELLS:
        y_wr = np.array([wr.rate_dme(e * nu.keV, int(sh[0]), wr.l_to_number[sh[1]], m * nu.GeV / nu.c0**2,
                                     sigma_e * nu.cm**2) for e in E_wr]) / WR_UNIT * rho / 0.3
        out[sh] = (E, r['rate_shell'][rapidd.XE_SHELLS.index(sh)] * 1e-6, E_wr, y_wr)
    return out


def roi_rate(spec, which, E_lo=E_ROI[0], E_hi=E_ROI[1], convention='electron'):
    '''rate [events/kg/day] in [E_lo, E_hi] keV of E_ER = E_e + E_sec, summed over the shells'''
    E_ER = np.geomspace(E_lo, E_hi, 400)
    tot = np.zeros_like(E_ER)
    for sh, (E, y, E_wr, y_wr) in spec.items():
        s = rapidd.XE_SHELLS.index(sh)
        E_sec = (rapidd.XE_BINDING_ENERGIES if convention == 'deposited' else rapidd.XE_DEEXCITATION_ENERGIES)[s] * 1e6
        x, yy = (E, y) if which == 'rapidd' else (E_wr, y_wr)
        tot += np.interp(E_ER, x + E_sec, yy, left=0., right=0.)
    return np.trapezoid(tot, E_ER)


def ratios(masses):
    '''S_RAPIDD / S_wimprates in the region of interest (and with E_ER > 0.1 keV)'''
    use_halo(rapidd.halo_path)
    out = []
    for m in masses:
        spec = spectra(m)
        out.append((roi_rate(spec, 'rapidd') / roi_rate(spec, 'wimprates'),
                    roi_rate(spec, 'rapidd', E_lo=0.1) / roi_rate(spec, 'wimprates', E_lo=0.1)))
    return np.array(out)


def main():
    rapidd.read_atomic_responses()
    rapidd.set_er_energy_convention('electron')
    m, lim = np.loadtxt(LIMIT_IN, unpack=True)
    r = ratios(m)
    with open(LIMIT_OUT, 'w') as f:
        f.write('# XENONnT S2-only DM-electron limit (arXiv:2601.11296, 90% CL after the power constraint, heavy mediator\n'
                '# F_DM = 1) rescaled to the RAPIDD atomic responses: sigma_e = XENONnT_DMe_SI.dat x S_wimprates/S_RAPIDD,\n'
                '# with S the rate in E_ER = E_e + E_deexcitation in [0.04, 0.7] keV (5 outer shells, same halo).\n'
                '# Approximate: the XENONnT efficiency and charge yield are not included. Written by\n'
                '# examples/xenonnt_dme_limit_rapidd_ff.py.\n'
                '# DM mass [GeV]   sigma_e [cm^2]   S_RAPIDD/S_wimprates [0.04, 0.7] keV   (same, [0.1, 0.7] keV)\n')
        for mi, li, (ri, ri01) in zip(m, lim, r):
            f.write('%.6e %.6e %.4f %.4f\n' % (mi, li / ri, ri, ri01))
            print('m = %6.3f GeV: S_RAPIDD/S_wimprates = %.3f ([0.1, 0.7] keV: %.3f), limit x %.2f'
                  % (mi, ri, ri01, 1 / ri))
    rapidd.reset_coefficients()


if __name__ == '__main__':
    main()
