"""Crystal response W_1 of silicon for the RAPIDD DM-electron rate, from the QCDark2 RPA dielectric function.

Usage:
    python qcdark2_to_rapidd.py <QCDark2>/dielectric_functions/composite/Si_comp.h5
        <QCDark2>/dielectric_functions/nolfe/Si_nolfe.h5 [--commit HASH]
        [--out-dir ../] [--reference ../../../tests/data/si_qcdark2_reference.txt]

Reads only the HDF5 datasets q [alpha m_e], E [eV], epsilon (q x E, complex) and the attributes V_cell [bohr^3],
M_cell [eV] (requires h5py; the QCDark2 Python code is not used). Writes, on the (Delta E, q) grid of the files,

    Si_W1.txt          W_1(q, Delta E) = Delta E q^2 V_cell Im eps(q, Delta E)       / (pi^2 alpha)   (unscreened)
    Si_W1_screened.txt W_1(q, Delta E) = Delta E q^2 V_cell Im[-1/eps(q, Delta E)]   / (pi^2 alpha)   (RPA screened)
    Si_W1_nolfe.txt    W_1(q, Delta E) = Delta E q^2 V_cell Im eps_noLFE(q, Delta E) / (pi^2 alpha)

with eps the composite dielectric function (Si_comp.h5, local field effects, LFE, for q <= 8 alpha m_e) and eps_noLFE
the one without LFE (Si_nolfe.h5, only the diagonal G = G' elements, epsilon^noLFE of arXiv:2603.12326, Sec. 2.4).
Im eps_noLFE = -U Im chi^0 is the independent-particle response, used for the couplings to the electron spin: in a
non-spin-polarised material they are not screened and their response is the single-particle one, as in
arXiv:2402.06817 (local-field corrections neglected, Sec. 3.1; unscreened responses from single-particle states,
Sec. 4.1).

in natural units (dimensionless), i.e. the first crystal response of arXiv:2105.02233 (Eq. R_crystal_2D), obtained
from |f_crystal|^2 = q^3 W_1/(8 Delta E alpha m_e^2) (arXiv:2105.02233, App. B) and
|f_crystal|^2 = q^5 V_cell Im eps/(8 pi^2 (alpha m_e)^2) (arXiv:2104.12786, Eq. 16).

References: QCDark2, C. Dreyer, R. Essig, M. Fernandez-Serra, M. Hott, A. Singal, arXiv:2603.12326, PRD 114, 035003
(2026), github.com/meganhott/QCDark2 (CC0); R. Catena, T. Emken, M. Matas, N. A. Spaldin, E. Urdshals,
arXiv:2105.02233, Phys. Rev. Research 3, 033149 (2021); R. Catena, N. A. Spaldin, arXiv:2402.06817, Phys. Rev.
Research 6, 033230 (2024); DarkELF, S. Knapen, J. Kozaczuk, T. Lin, arXiv:2104.12786, PRD 105, 015014 (2022). Im eps and Im[-1/eps] smaller
than CLIP_REL times the maximum of their q row (numerical noise below the gap, |Im eps| < 1e-16) are set to zero.

With --reference, also writes the heavy- and light-mediator dR/dE [events/kg/year/eV] of the QCDark2 rate formula
(arXiv:2603.12326, Eqs. 6 and 10: Simpson's rule on the q grid of the file, halo and sigma_e of REF_ASTRO, the
QCDark2 M_cell) reimplemented here from the original eps (screened and unscreened from Si_comp.h5, unscreened from
Si_nolfe.h5), used by the RAPIDD regression tests.
"""
import argparse
import os

import h5py
import numpy as np
from scipy.integrate import simpson
from scipy.special import erf

ALPHA = 1.0 / 137.03599908      # as in QCDark2
M_E_EV = 5.1099894e5            # eV, as in QCDark2
C_KMS = 299792.458
KG_EV = 5.609588e35             # eV per kg
SEC_PER_YEAR = 60. * 60. * 24 * 365.25
CLIP_REL = 1e-10
FMT = "%.6e"
REF_ASTRO = {'v0': 238., 'vEarth': 250.2, 'vEscape': 544., 'rhoX': 0.3e9, 'sigma_e': 1e-37}  # QCDark2 default_astro
REF_MASSES_EV = (1e6, 3e6, 1e7, 1e8, 1e9)


def read_eps(path):
    with h5py.File(path, 'r') as f:
        q = np.asarray(f['q'][()], dtype=float)
        E = np.asarray(f['E'][()], dtype=float)
        eps = np.asarray(f['epsilon'][()], dtype=complex)
        V_cell = float(f.attrs['V_cell'])
        M_cell = float(f.attrs['M_cell'])
    if eps.shape != (q.size, E.size):
        raise ValueError("epsilon has shape %s, expected (len(q), len(E)) = %s" % (eps.shape, (q.size, E.size)))
    return q, E, eps, V_cell, M_cell


def clip_noise(x):
    '''Zero for negative values and for values below CLIP_REL times the maximum of the q row'''
    rowmax = np.max(x, axis=1, keepdims=True)
    return np.where(x > CLIP_REL * rowmax, x, 0.)


def w1_tables(q, E, eps, V_cell):
    '''W_1 unscreened and RPA screened (dimensionless) on the (q, E) grid: q in alpha m_e, E in eV, V_cell in bohr^3'''
    qe = q * ALPHA * M_E_EV                  # eV
    Vc = V_cell / (ALPHA * M_E_EV)**3        # eV^-3
    pref = E[None, :] * qe[:, None]**2 * Vc / (np.pi**2 * ALPHA)
    return pref * clip_noise(np.imag(eps)), pref * clip_noise(np.imag(-1. / eps))


def eta_mb(q_ev, E, m_X, astro):
    '''Inverse mean speed [1/c] of the truncated Maxwell-Boltzmann distribution (as get_eta_MB of QCDark2)'''
    vesc, ve, v0 = astro['vEscape'] / C_KMS, astro['vEarth'] / C_KMS, astro['v0'] / C_KMS
    vmin = q_ev[:, None] / (2. * m_X) + E[None, :] / q_ev[:, None]
    a = -4. * ve * np.exp(-(vesc / v0)**2) + np.sqrt(np.pi) * v0 * (erf((vmin + ve) / v0) - erf((vmin - ve) / v0))
    b = (-2. * (ve + vesc - vmin) * np.exp(-(vesc / v0)**2)
         + np.sqrt(np.pi) * v0 * (erf(vesc / v0) - erf((vmin - ve) / v0)))
    val = np.where(vmin < vesc - ve, a, np.where(vmin < vesc + ve, b, 0.))
    K = v0**3 * (-2. * np.pi * (vesc / v0) * np.exp(-(vesc / v0)**2) + np.pi**1.5 * erf(vesc / v0))
    return v0**2 * np.pi / (2. * ve * K) * val


def qcdark2_dRdE(q, E, eps, V_cell, M_cell, m_X, mediator, screening, astro=REF_ASTRO):
    '''
    dR/dE [events/kg/year/eV] of the QCDark2 formula: (V_cell/M_cell) n_X sigma_e/(4 pi mu^2) int dq q F_DM^2 S eta,
    S = q^2 Im eps/(2 pi alpha) * |eps/eps_screen|^2 (screening "RPA": eps_screen = eps; None: 1), F_DM = 1 (heavy)
    or (alpha m_e/q)^2 (light). Same grid, quadrature and units as get_dR_dE of QCDark2.
    '''
    q_ev = q * ALPHA * M_E_EV
    rho_T = M_cell / KG_EV / V_cell                                   # kg/bohr^3
    mu = m_X * M_E_EV / (m_X + M_E_EV)
    pref = (1. / rho_T) * (astro['rhoX'] / m_X) * (astro['sigma_e'] / mu**2) / (4. * np.pi)
    F = np.ones_like(q) if mediator == 'heavy' else 1. / q**2
    if screening == 'RPA':
        loss = np.imag(-1. / eps)
    elif screening is None:
        loss = np.imag(eps)
    else:
        raise ValueError("screening must be 'RPA' or None")
    S = loss * q[:, None]**2 / (2. * np.pi * ALPHA)
    integrand = q[:, None] * F[:, None]**2 * S * eta_mb(q_ev, E, m_X, astro)
    cm2sec = 1. / C_KMS * 1e-5
    return pref * simpson(integrand, x=q_ev, axis=0) / cm2sec * SEC_PER_YEAR


def write_table(path, W, q, E, title, commit, V_cell, M_cell, source):
    with open(path, 'w') as f:
        f.write("// %s of silicon (crystal response W_1 of arXiv:2105.02233, dimensionless, natural units)\n" % title)
        f.write("// Source: QCDark2, C. Dreyer, R. Essig, M. Fernandez-Serra, M. Hott, A. Singal, arXiv:2603.12326 "
                "(github.com/meganhott/QCDark2, CC0, commit %s), dielectric_functions/%s\n" % (commit, source))
        f.write("// V_cell = %.10g bohr^3, M_cell(QCDark2) = %.10g eV\n" % (V_cell, M_cell))
        f.write("// Rows: Delta E [eV] = %g + %g i, i = 0..%d\n" % (E[0], E[1] - E[0], E.size - 1))
        f.write("// Columns: q [alpha m_e] = %g + %g j, j = 0..%d\n" % (q[0], q[1] - q[0], q.size - 1))
        for iE in range(E.size):
            f.write("\t".join(FMT % w if w != 0. else "0" for w in W[:, iE]) + "\n")


def check_grid(q, E):
    dq, dE = np.diff(q), np.diff(E)
    if not (np.allclose(dq, dq[0], rtol=1e-6) and np.allclose(dE, dE[0], rtol=1e-6)):
        raise ValueError("the q and E grids must be linear")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('h5', help="QCDark2 composite dielectric function (Si_comp.h5)")
    p.add_argument('h5_nolfe', help="QCDark2 dielectric function without LFE (Si_nolfe.h5)")
    p.add_argument('--commit', default="unknown", help="QCDark2 commit hash, written in the headers")
    p.add_argument('--out-dir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    p.add_argument('--reference', default=None, help="also write the QCDark2-formula reference rates to this file")
    a = p.parse_args()
    q, E, eps, V_cell, M_cell = read_eps(a.h5)
    check_grid(q, E)
    q0, E0, eps0, V0, M0 = read_eps(a.h5_nolfe)
    if not (np.array_equal(q0, q) and np.array_equal(E0, E) and V0 == V_cell and M0 == M_cell):
        raise ValueError("the two dielectric functions must have the same grid, V_cell and M_cell")
    W1, W1s = w1_tables(q, E, eps, V_cell)
    W1n, _ = w1_tables(q, E, eps0, V_cell)
    src = "composite/Si_comp.h5 (RPA, LFE for q <= 8 alpha m_e)"
    write_table(os.path.join(a.out_dir, 'Si_W1.txt'), W1, q, E, "W_1 (unscreened, Im eps)", a.commit, V_cell, M_cell,
                src)
    write_table(os.path.join(a.out_dir, 'Si_W1_screened.txt'), W1s, q, E, "W_1 (RPA screened, Im[-1/eps])",
                a.commit, V_cell, M_cell, src)
    write_table(os.path.join(a.out_dir, 'Si_W1_nolfe.txt'), W1n, q, E,
                "W_1 (independent-particle response, Im eps without LFE)", a.commit, V_cell, M_cell,
                "nolfe/Si_nolfe.h5 (RPA, no LFE)")
    if a.reference:
        cols, names = [E], ["E_eV"]
        for m in REF_MASSES_EV:
            for med in ('heavy', 'light'):
                for scr in ('RPA', None):
                    cols.append(qcdark2_dRdE(q, E, eps, V_cell, M_cell, m, med, scr))
                    names.append("m%.0e_%s_%s" % (m, med, 'RPA' if scr else 'none'))
                cols.append(qcdark2_dRdE(q, E, eps0, V_cell, M_cell, m, med, None))
                names.append("m%.0e_%s_noLFE" % (m, med))
        hdr = ("dR/dE [events/kg/year/eV] of the QCDark2 rate formula (arXiv:2603.12326 Eqs. 6, 10, reimplemented in "
               "lib/crystal_responses/generate/qcdark2_to_rapidd.py) from Si_comp.h5 (RPA, none) and Si_nolfe.h5 (noLFE, "
               "unscreened) (QCDark2 commit %s)\n"
               "halo: v0 = %g, vEarth = %g, vesc = %g km/s, rhoX = %g eV/cm^3; sigma_e = %g cm^2; M_cell = %.10g eV; "
               "light mediator: F_DM = (alpha m_e/q)^2\n%s"
               % (a.commit, REF_ASTRO['v0'], REF_ASTRO['vEarth'], REF_ASTRO['vEscape'], REF_ASTRO['rhoX'],
                  REF_ASTRO['sigma_e'], M_cell, " ".join(names)))
        np.savetxt(a.reference, np.column_stack(cols), fmt="%.8e", header=hdr)


if __name__ == '__main__':
    main()
