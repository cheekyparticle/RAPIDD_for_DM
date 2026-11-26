"""DAMIC-M DM-electron limit (arXiv:2503.14617, heavy mediator, F_DM = 1) rescaled to the RAPIDD silicon response
(QCDark2 RPA dielectric function with local field effects, lib/crystal_responses).

The published limit (lib/experiments/DAMIC-M/DAMIC-M_2025_DMe_heavy_QCDark.dat) is computed with the QCDark crystal
form factor (QCDark1, C. E. Dreyer, R. Essig, M. Fernandez-Serra, A. Singal, C. Zhen, arXiv:2306.14944) and the
analytic "modified Thomas-Fermi" screening of G. Cappellini, R. Del Sole, L. Reining, F. Bechstedt (PRB 47, 9892
(1993): eps_0 = 11.3, tau = 1.563, omega_p = 16.6 eV, q_TF = 4.13 keV, |f_crystal|^2 / eps_MTF^2), the
Ramanathan-Kurinsky charge yield (arXiv:2004.10709, 100 K) and the standard halo of arXiv:2105.00599 (Baxter et
al.). DAMIC-M Collaboration, arXiv:2503.14617, PRL 135, 071002 (2025); limits from
github.com/DAMIC-M/LBC_2025_HSBenchmark_Limits (CC BY 4.0). The pattern
response of DAMIC-M is not released, so the limit cannot be recomputed. It is rescaled by the ratio of the signals,

    sigma_e^RAPIDD(m) = sigma_e^DAMIC-M(m) * S_DAMIC-M(m) / S_RAPIDD(m),

with S = sum_{Q=2,3,4} eps_Q R_Q, R_Q the rate with Q electron-hole pairs and eps_Q = 0.38, 0.65, 0.79 the
efficiencies of a Q e- deposit to be found as a pattern (arXiv:2503.14617; the candidates have 2-4 e-). S_DAMIC-M
is the QCDark1 rate with MTF screening, reimplemented here from the QCDark1 formula (Eq. 13 of arXiv:2306.14944,
d_rate and d_rate_RamanathanQ of github.com/asingal14/QCDark) and its published Si form factor; S_RAPIDD is the
RAPIDD rate (rapidd.difrate_electron(target="Si"), RPA screening). Both use the RAPIDD default halo (SHM of
arXiv:2105.00599, the halo of DAMIC-M) and the same charge yield. Not modelled: charge diffusion, pattern
identification details beyond eps_Q, the profile likelihood with the observed {31} candidate and the backgrounds
(the rescaling assumes they act in the same way on the two signals), Earth shielding (verne; it matters only at the
low-mass end, where the limit rises steeply).

Masses whose kinematic end point (1/2) m (v_esc + v_E)^2 is below the onset of Q = 4 in the charge yield table
(6.75 eV, m < 1.9 MeV) are not written: there the signal comes from the end point next to the Q = 2, 3 onsets and the
ratio depends on grid and halo-tail details (it varies between 0.2 and 200 at 0.75-1.6 MeV).

Diagnostic column: S_QCDark2-MTF / S_DAMIC-M, with the QCDark2 Im eps (no LFE screening) and the MTF screening,
separates the change of the form factor (QCDark2 vs QCDark1) from the change of the screening (RPA + LFE vs MTF).

Requires the QCDark1 silicon form factor (results/Si_final.hdf5 of github.com/asingal14/QCDark, commit efdb7d1)
and h5py; run from the repository root with the library built:
    python examples/damicm_dme_limit_rapidd.py <QCDark>/results/Si_final.hdf5
It writes lib/experiments/DAMIC-M/DAMIC-M_2025_DMe_heavy_QCDark_RAPIDD.dat.
"""
import os
import sys
import warnings

import h5py
import numpy as np
from scipy.integrate import simpson

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import rapidd  # noqa: E402
from rapidd.electron import (coeff_from_sigma_e, difrate_electron, difrate_dQ, read_ionization_table,  # noqa: E402
                             CrystalTruncationWarning)

ALPHA = 1. / 137.03599908     # QCDark constants
M_E_EV = 5.1099894e5
C_KMS = 299792.458
KG_EV = 5.609588e35
CM2SEC = 1. / C_KMS * 1e-5
SEC2YR = 1. / (60. * 60. * 24 * 365.25)
MTF_SI = {'eps0': 11.3, 'tau': 1.563, 'omegaP': 16.6, 'qTF': 4.13e3}   # Cappellini et al., Si (eV)
EFF = {2: 0.38, 3: 0.65, 4: 0.79}                                     # arXiv:2503.14617
SIGMA_E = 1e-37
LIMIT_IN = os.path.join(rapidd.core.experiments_dir, 'DAMIC-M', 'DAMIC-M_2025_DMe_heavy_QCDark.dat')
LIMIT_OUT = os.path.join(rapidd.core.experiments_dir, 'DAMIC-M', 'DAMIC-M_2025_DMe_heavy_QCDark_RAPIDD.dat')


def eps_mtf(q, E, p=MTF_SI):
    '''modified Thomas-Fermi dielectric function (real), q and E in eV'''
    val = (1. / (p['eps0'] - 1.) + p['tau'] * (q / p['qTF'])**2 + q**4 / (4. * M_E_EV**2 * p['omegaP']**2)
           - (E / p['omegaP'])**2)
    return 1. + 1. / val


class Halo:
    '''eta(v_min) [1/c] from a RAPIDD halo table (v_min [km/s], eta [(km/s)^-1])'''
    def __init__(self, path=rapidd.halo_path):
        h = np.loadtxt(path, skiprows=1)
        self.v, self.eta = h[:, 0], h[:, 1]

    def __call__(self, vmin_c):
        return np.interp(vmin_c * C_KMS, self.v, self.eta, left=self.eta[0], right=0.) * C_KMS


def qcdark1_dRdE(path, m_X, halo, rho=0.3e9, sigma=SIGMA_E):
    '''
    QCDark1 heavy-mediator dR/dE [events/kg/year/eV] at the E bin centres with MTF screening:
    dR/dlnE = n_X (kg/M_cell) sigma alpha (m_e/mu)^2 int dq E/q^2 eta |f_crystal|^2 / eps_MTF^2 (Simpson on the q bin
    centres), as d_rate of QCDark1.
    '''
    with h5py.File(path, 'r') as f:
        dq = f['results'].attrs['dq'] * ALPHA * M_E_EV
        dE = f['results'].attrs['dE']
        M_cell = f['results'].attrs['mCell']
        f2 = f['results/f2'][...]
    q = dq * np.arange(f2.shape[0]) + dq / 2.
    E = dE * np.arange(f2.shape[1]) + dE / 2.
    mu = m_X * M_E_EV / (m_X + M_E_EV)
    vmin = q[:, None] / (2. * m_X) + E[None, :] / q[:, None]
    with np.errstate(divide='ignore'):
        scr = 1. / eps_mtf(q[:, None], E[None, :])**2
    integrand = E[None, :] / q[:, None]**2 * halo(vmin) * f2 * scr
    pref = (rho / m_X) * (KG_EV / M_cell) * sigma * ALPHA * (M_E_EV / mu)**2
    return E, pref * simpson(integrand, x=q, axis=0) / CM2SEC / SEC2YR / E


def qcdark2_mtf_dRdE(m_X, halo, rho=0.3e9, sigma=SIGMA_E):
    '''Same rate with W_1 (Im eps) of the RAPIDD tables (QCDark2) and the MTF screening, per kg of 2 x 28.0855 u'''
    W = np.loadtxt(os.path.join(rapidd.crystal_resp_dir, 'Si_W1.txt'), comments='//').T    # (q, E)
    q = (0.01 + 0.02 * np.arange(W.shape[0])) * ALPHA * M_E_EV
    E = 0.1 * np.arange(W.shape[1])
    mu = m_X * M_E_EV / (m_X + M_E_EV)
    vmin = q[:, None] / (2. * m_X) + E[None, :] / np.maximum(q[:, None], 1e-30)
    with np.errstate(divide='ignore'):
        scr = 1. / eps_mtf(q[:, None], E[None, :])**2
    I = simpson(q[:, None] * halo(vmin) * W * scr, x=q, axis=0)
    M_cell = rapidd.CRYSTAL_CELL_MASS_AU['Si'] * rapidd.amu * 1e9
    rho_n, sig_n = rho / 50677.3093773**3, sigma * 50677.3093773**2
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.where(E > 0, (rho_n / m_X) * sig_n / (8. * mu**2) * I / M_cell / E, 0.)
    return E, r * 50677.3093773 * C_KMS * 1e5 / SEC2YR * KG_EV


def rate_per_Q(E, dRdE, dE, n_max=4):
    '''R_Q [events/kg/year] (rectangle rule on the E grid, p_Q interpolated, 1.1 <= E <= 50 eV), as QCDark'''
    Et, p = read_ionization_table()
    E = np.round(E, 5)
    sel = (E >= Et.min()) & (E <= Et.max())
    return np.array([np.sum(dRdE[sel] * dE * np.interp(E[sel], Et, p[n])) for n in range(n_max)])


def signal(RQ):
    return sum(eff * RQ[Q - 1] for Q, eff in EFF.items())


def rapidd_RQ(m):
    rapidd.reset_coefficients()
    rapidd.set_any_Ecoeff(coeff_from_sigma_e(SIGMA_E, m), 1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', CrystalTruncationWarning)
        r = difrate_electron(0.3, m, target='Si')
    rapidd.reset_coefficients()
    return difrate_dQ(r, Q_max=4)['dRdQ'] * 365.25


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    qcdark1 = sys.argv[1]
    rapidd.read_halo(rapidd.halo_path)
    rapidd.read_crystal_responses()
    halo = Halo()
    m, lim = np.loadtxt(LIMIT_IN, unpack=True)
    Et, p = read_ionization_table()
    E_onset = Et[np.argmax(p[max(EFF) - 1] > 0)]            # onset of the largest Q of the window [eV]
    v_max = (halo.v[halo.eta > 0].max()) / C_KMS             # end of the halo table (v_esc + v_E) [c]
    rows = []
    for mi, li in zip(m, lim):
        if 0.5 * mi * 1e9 * v_max**2 < E_onset:
            print('m = %8.5f GeV: skipped (kinematic end point below the Q = %d onset, %.2f eV)' % (mi, max(EFF), E_onset))
            continue
        m_ev = mi * 1e9
        E1, r1 = qcdark1_dRdE(qcdark1, m_ev, halo)
        S_dm = signal(rate_per_Q(E1, r1, E1[1] - E1[0]))
        E2, r2 = qcdark2_mtf_dRdE(m_ev, halo)
        S_mtf = signal(rate_per_Q(E2, r2, 0.1))
        S_rp = signal(rapidd_RQ(mi))
        ratio = S_rp / S_dm if S_dm > 0 else np.nan
        rows.append((mi, li / ratio if ratio > 0 else np.nan, ratio, S_mtf / S_dm if S_dm > 0 else np.nan))
        print('m = %8.5f GeV: S_RAPIDD/S_DAMIC-M = %.3f (QCDark2 Im eps + MTF / DAMIC-M: %.3f), limit x %.2f'
              % (mi, ratio, rows[-1][3], 1. / ratio))
    with open(LIMIT_OUT, 'w') as f:
        f.write('# DAMIC-M DM-electron limit (arXiv:2503.14617, 90%% CL, heavy mediator F_DM = 1, QCDark1 + analytic MTF\n'
                '# screening) rescaled to the RAPIDD silicon response (QCDark2 RPA dielectric function with local field\n'
                '# effects, arXiv:2603.12326): sigma_e = DAMIC-M_2025_DMe_heavy_QCDark.dat x S_DAMIC-M/S_RAPIDD, with\n'
                '# S = sum_{Q=2..4} eps_Q R_Q (eps_Q = 0.38, 0.65, 0.79 pattern efficiencies), Ramanathan-Kurinsky 100 K\n'
                '# charge yield and the same halo (SHM of arXiv:2105.00599) for both. S_DAMIC-M reimplements the QCDark1\n'
                '# rate with its published Si form factor (github.com/asingal14/QCDark efdb7d1, results/Si_final.hdf5).\n'
                '# Approximate: charge diffusion, pattern details, backgrounds and Earth shielding are not modelled.\n'
                '# Masses with (1/2) m (v_esc + v_E)^2 below the Q = 4 onset (%.2f eV, m < 1.9 MeV) are omitted: there the\n'
                '# signal comes from the kinematic end point and the ratio is not reliable.\n'
                '# Written by examples/damicm_dme_limit_rapidd.py.\n'
                '# DM mass [GeV]   sigma_e [cm^2]   S_RAPIDD/S_DAMIC-M   S_QCDark2-MTF/S_DAMIC-M (diagnostic)\n'
                % E_onset)
        for row in rows:
            f.write('%.6e %.6e %.4f %.4f\n' % row)


if __name__ == '__main__':
    main()
