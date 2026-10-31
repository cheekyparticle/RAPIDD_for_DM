"""DM-electron scattering in Xenon for the NREFT operators of arXiv:1912.08204
(Catena, Emken, Spaldin, Tarantino).

Conventions
-----------
- Amplitude: M = sum_i c_i(q) <O_i>, c_i(q) = c_i^s + (q_ref/q)^2 c_i^l, with q_ref = alpha m_e,
  dimensionless c_i^s (contact) and c_i^l (long-range). Operators O_1, O_3 ... O_15 of Tab. I
  (DM spin 1/2); O_2, O_16 and O_17 ... O_20 do not enter the DM responses of the paper.
- jchi is the DM spin (MadDM dm_spin = 2 jchi + 1): the operators with S_chi vanish for jchi = 0.
- The rate is dR/dlnE = n_chi/(128 pi mchi^2 m_e^2) int dq q sum_j R_j W_j (Eqs. (ionization
  spectrum) and (RW)), with the atomic responses W_1 ... W_4 of the 11 Xenon shells.
- Units: masses in GeV, energies in keV, rates in events/kg/day(/keV), rhoDM in GeV/cm^3.
- The halo table must contain halo(vmin, 2) (i >= 2) for the operators depending on v_el^perp.
"""
import numpy as np

from .core import (set_any_Ecoeff, read_atomic_responses, difrate_electron_dEr_all, difrate_dS2,
                   calc_halo, read_halo, M_ELECTRON, XE_SHELLS, N_SHELL_XE, N_ATOMIC_RESP,
                   ELECTRON_OPERATORS, ELECTRON_INTERACTIONS)
from .stats import poisson_upper_limit

_trapezoid = getattr(np, "trapezoid", None) or np.trapz  # numpy < 2 has only trapz

ALPHA_EM = 7.2973525643e-3                   # fine-structure constant (CODATA 2022)
E_CHARGE = np.sqrt(4. * np.pi * ALPHA_EM)    # electron charge (natural units)
GEV2_TO_CM2 = 3.893793721e+8 * 1.e-36        # (hbar c)^2: GeV^-2 -> cm^2, as in phys_consts.h
XE_SHELLS_PAPER = ("4s", "4p", "4d", "5s", "5p")  # Xenon shells used in arXiv:1912.08204
ELECTRON_MODELS = ("anapole", "magnetic_dipole", "electric_dipole")
ELECTRON_EXPERIMENTS = ("Xenon10", "Xenon1T")

# halo parameters of arXiv:1912.08204 (App. E)
PAPER_HALO = {'rhoDM': 0.4, 'v0': 220., 'vesc': 544., 've': 244.}


def q_ref(m_e=M_ELECTRON):
    '''Reference momentum q_ref = alpha m_e [GeV] of the long-range couplings'''
    return ALPHA_EM * m_e


def reduced_mass(m1, m2):
    return m1 * m2 / (m1 + m2)


def write_paper_halo(table_path, profile="SHM", i=2):
    '''
    Write the halo table of arXiv:1912.08204 (SHM with v0 = 220 km/s, vesc = 544 km/s and
    Earth speed v_E = 244 km/s) with the integrals up to halo(vmin, i).
    The local density of the paper is PAPER_HALO['rhoDM'] = 0.4 GeV/cm^3.
    '''
    calc_halo(table_path, profile=profile, vesc=PAPER_HALO['vesc'], v0=PAPER_HALO['v0'],
              ve=PAPER_HALO['ve'], i=i)


def _shell_indices(shells):
    if shells is None:
        return list(range(N_SHELL_XE))
    idx = []
    for s in shells:
        if isinstance(s, str):
            if s not in XE_SHELLS:
                raise ValueError(f"Unknown Xenon shell {s}: choose from {XE_SHELLS}")
            idx.append(XE_SHELLS.index(s))
        else:
            if not 0 <= int(s) < N_SHELL_XE:
                raise ValueError(f"Shell index {s} out of range [0, {N_SHELL_XE - 1}]")
            idx.append(int(s))
    return sorted(set(idx))


def resp_mask(resp=None):
    '''Bit mask of the terms R_j W_j: None -> all, j -> only R_j W_j, iterable -> the listed ones'''
    if resp is None:
        return 2**N_ATOMIC_RESP - 1
    js = [resp] if np.ndim(resp) == 0 else list(resp)
    mask = 0
    for j in js:
        if not 1 <= int(j) <= N_ATOMIC_RESP:
            raise ValueError(f"Response index {j} must be in [1, {N_ATOMIC_RESP}]")
        mask |= 1 << (int(j) - 1)
    return mask


def difrate_electron(rhoDM, mDM, jchi=0.5, shells=None, resp=None, m_e=M_ELECTRON):
    '''
    Differential ionization rate of Xenon for the DM-electron coefficients set with
    set_any_Ecoeff (contact and long-range). Requires read_halo and read_atomic_responses.

    Parameters
    ----------
    rhoDM : float
        Local DM density [GeV/cm^3].
    mDM : float
        DM mass [GeV].
    jchi : float
        DM spin (1/2 for the operators of arXiv:1912.08204, 0 for scalar DM).
    shells : iterable of str/int or None
        Xenon shells to include (names as in XE_SHELLS or indices). None: all the 11 shells
        (as in MadDM); XE_SHELLS_PAPER: the 5 outer shells used in arXiv:1912.08204.
    resp : int, iterable of int or None
        Keep only the terms R_j W_j with j in resp (None: all). Single terms can be negative.
    m_e : float
        Electron mass [GeV].

    Returns
    -------
    dict with
        k, E : outgoing electron momentum and energy [keV] (GRIDSIZE_K points)
        dRdE, dRdlnE : total over the selected shells [events/kg/day/keV], [events/kg/day]
        dRdE_shell, dRdlnE_shell : dict shell name -> array, for the selected shells
        E_GeV, rate_shell : grids [GeV] (GRIDSIZE_K+1) and per-shell dR/dE [events/kg/day/GeV]
            of all the 11 shells (zero for the shells not selected), input of difrate_dS2
    '''
    if mDM <= 0 or rhoDM < 0:
        raise ValueError("mDM must be > 0 and rhoDM >= 0")
    idx = _shell_indices(shells)
    k, E, rate_shell = difrate_electron_dEr_all(rhoDM, mDM, jchi, m_e, resp_mask=resp_mask(resp))
    keep = np.zeros(N_SHELL_XE, dtype=bool)
    keep[idx] = True
    rate_shell = np.where(keep[:, None], rate_shell, 0.)
    E_keV = E[:-1] * 1e6
    dRdE_shell = {XE_SHELLS[s]: rate_shell[s] / 1e6 for s in idx}
    dRdE = rate_shell.sum(axis=0) / 1e6
    return {'k': k[:-1] * 1e6, 'E': E_keV, 'dRdE': dRdE, 'dRdlnE': E_keV * dRdE,
            'dRdE_shell': dRdE_shell,
            'dRdlnE_shell': {s: E_keV * r for s, r in dRdE_shell.items()},
            'E_GeV': E, 'rate_shell': rate_shell}


def total_rate_electron(rhoDM, mDM, jchi=0.5, shells=None, resp=None, m_e=M_ELECTRON, rate=None):
    '''
    Total ionization rate [events/kg/day]: integral of dR/dlnE over ln E on the grid of
    difrate_electron (trapezoidal rule in ln E, from E = 0.01 eV to the kinematic end point).
    If rate (the output of difrate_electron) is given, the other arguments are ignored.
    Note: the XENON10 response uses the legacy left sum dR/dE(E_i) (E_{i+1} - E_i), which on
    this log grid is ~7% larger.
    '''
    if rate is None:
        rate = difrate_electron(rhoDM, mDM, jchi=jchi, shells=shells, resp=resp, m_e=m_e)
    return float(_trapezoid(rate['dRdlnE'], np.log(rate['E'])))


def electron_signal(experiment, rhoDM=None, mDM=None, jchi=0.5, shells=None, m_e=M_ELECTRON, rate=None):
    '''
    S2 signal of the DM-electron scattering for the coefficients set.

    Parameters
    ----------
    experiment : {"Xenon10", "Xenon1T"} (case insensitive)
        Xenon10: 7 S2 bins (15 kg day, arXiv:1703.00910 detector model).
        Xenon1T: S2-only data release response (s2_response_er.csv) and the mass
        dependent ROI of 5d_results_dmelectron.dat (signal = -1 outside 0.02-5 GeV).
    rhoDM, mDM, jchi, shells, m_e : see difrate_electron (mDM is always needed).
    rate : dict or None
        Output of difrate_electron (computed if None).

    Returns
    -------
    dict of difrate_dS2: S2, dRdS2, signal, n_obs (+ dRdS2_shell for Xenon10, n_bkg for
    Xenon1T), plus the rate used.
    '''
    names = {e.lower(): e for e in ELECTRON_EXPERIMENTS}
    if str(experiment).lower() not in names:
        raise ValueError(f"experiment must be one of {ELECTRON_EXPERIMENTS} (case insensitive)")
    experiment = names[str(experiment).lower()]
    if mDM is None:
        raise ValueError("mDM is needed")
    if rate is None:
        if rhoDM is None:
            raise ValueError("rhoDM is needed when rate is not given")
        rate = difrate_electron(rhoDM, mDM, jchi=jchi, shells=shells, m_e=m_e)
    res = difrate_dS2(experiment, rate['E_GeV'], rate['rate_shell'], mDM)
    res['rate'] = rate
    return res


########################################################################################
# Models of arXiv:1912.08204 (Sec. IV and App. D)
########################################################################################

def sigma_e_from_coeff(coeff, mDM, m_e=M_ELECTRON):
    '''Reference cross section sigma_e = mu^2 c^2/(16 pi mDM^2 m_e^2) [cm^2] (Eq. (reference
    cross section)) of a dimensionless coefficient c'''
    mu = reduced_mass(mDM, m_e)
    return mu**2 * coeff**2 / (16. * np.pi * mDM**2 * m_e**2) * GEV2_TO_CM2


def coeff_from_sigma_e(sigma_e, mDM, m_e=M_ELECTRON):
    '''Dimensionless coefficient c >= 0 with reference cross section sigma_e [cm^2]'''
    if sigma_e < 0:
        raise ValueError("sigma_e must be >= 0")
    mu = reduced_mass(mDM, m_e)
    return np.sqrt(16. * np.pi * mDM**2 * m_e**2 * sigma_e / GEV2_TO_CM2) / mu


def set_dark_photon_coeff(mDM, sigma_e, interaction="contact", m_e=M_ELECTRON):
    '''
    Dark photon model: only O_1, with c_1^s (heavy mediator, "contact") or c_1^l (light
    mediator, "long-range") fixed by the reference cross section sigma_e [cm^2]. In terms of
    the Lagrangian parameters c_1^s = 4 mDM m_e g_D e eps/m_A'^2, c_1^l = 4 mDM m_e g_D e eps/q_ref^2.
    Returns the coefficient set.
    '''
    if interaction not in ELECTRON_INTERACTIONS:
        raise ValueError(f"interaction must be one of {ELECTRON_INTERACTIONS}")
    if mDM <= 0:
        raise ValueError("mDM must be > 0")
    c1 = coeff_from_sigma_e(sigma_e, mDM, m_e)
    set_any_Ecoeff(c1, 1, interaction=interaction)
    return c1


def electron_model_coeffs(model, mDM, coupling, m_e=M_ELECTRON, g_e=2.0):
    '''
    NREFT coefficients of the anapole, magnetic dipole and electric dipole DM models of
    arXiv:1912.08204 (Sec. IV.B, App. D), with e = sqrt(4 pi alpha) and q_ref = alpha m_e.

    Parameters
    ----------
    model : {"anapole", "magnetic_dipole", "electric_dipole"}
        anapole (Majorana DM): c_8^s = 8 e m_e mDM g/L^2, c_9^s = -4 g_e e m_e mDM g/L^2;
        magnetic_dipole (Dirac DM): c_1^s = 4 e m_e g/L, c_4^s = 8 g_e e mDM g/L,
            c_5^l = 16 e m_e^2 mDM/q_ref^2 g/L, c_6^l = -8 g_e e m_e^2 mDM/q_ref^2 g/L;
        electric_dipole (Dirac DM): c_11^l = 16 e mDM m_e^2/q_ref^2 g/L.
    mDM : float
        DM mass [GeV].
    coupling : float
        g/Lambda^2 [GeV^-2] for the anapole, g/Lambda [GeV^-1] for the dipoles.
    m_e : float
        Electron mass [GeV].
    g_e : float
        Electron g-factor (2 in arXiv:1912.08204).

    Returns
    -------
    dict (operator, interaction) -> coefficient
    '''
    if mDM <= 0:
        raise ValueError("mDM must be > 0")
    e, qr2 = E_CHARGE, q_ref(m_e)**2
    if model == "anapole":
        return {(8, "contact"): 8. * e * m_e * mDM * coupling,
                (9, "contact"): -4. * g_e * e * m_e * mDM * coupling}
    if model == "magnetic_dipole":
        return {(1, "contact"): 4. * e * m_e * coupling,
                (4, "contact"): 8. * g_e * e * mDM * coupling,
                (5, "long-range"): 16. * e * m_e**2 * mDM / qr2 * coupling,
                (6, "long-range"): -8. * g_e * e * m_e**2 * mDM / qr2 * coupling}
    if model == "electric_dipole":
        return {(11, "long-range"): 16. * e * mDM * m_e**2 / qr2 * coupling}
    raise ValueError(f"Model {model} not implemented. Choose from {ELECTRON_MODELS}.")


def set_any_E_EFT_coeff(model, mDM, coupling, m_e=M_ELECTRON, g_e=2.0, operators=None):
    '''
    Set the DM-electron coefficients of a model of arXiv:1912.08204 (see electron_model_coeffs).
    The other coefficients are not changed (call reset_coefficients first).
    operators: if given, set only the coefficients of these operators (e.g. to compute the
    contribution of a single operator, as in Fig. 9 of the paper).
    Returns the dict (operator, interaction) -> coefficient set.
    '''
    coeffs = electron_model_coeffs(model, mDM, coupling, m_e=m_e, g_e=g_e)
    if operators is not None:
        coeffs = {key: val for key, val in coeffs.items() if key[0] in operators}
    for (op, interaction), val in coeffs.items():
        set_any_Ecoeff(val, op, interaction=interaction)
    return coeffs


########################################################################################
# Exclusion limits (90% CL), using the scaling rate ~ coupling^2
########################################################################################

def electron_limit(signal, coupling=1.0, cl=0.9):
    '''
    Upper limit on a coupling g, given the signal computed for the value coupling, assuming
    rate ~ g^2 (one coefficient, or a model with a single overall coupling).

    Xenon10 (7 bins): Poisson limit without background in each bin independently, keeping the
    strongest one (as in arXiv:1912.08204). Xenon1T (ROI of the S2-only data release):
    Poisson limit with the expected background subtracted (as rapidd.stats.poisson_likelihood_limit).
    Returns np.inf if there is no signal and np.nan if the mass is outside the Xenon1T table.

    Parameters
    ----------
    signal : dict
        Output of electron_signal (or difrate_dS2).
    coupling : float
        Coupling used to compute the signal.
    cl : float
        Confidence level.
    '''
    sig = np.atleast_1d(np.asarray(signal['signal'], dtype=float))
    obs = np.atleast_1d(np.asarray(signal['n_obs']))
    if 'n_bkg' in signal:  # Xenon1T ROI
        if sig[0] < 0:
            return np.nan
        if sig[0] <= 0:
            return np.inf
        n_up = poisson_upper_limit(obs[0], cl) - signal['n_bkg']
        if n_up <= 0:
            return 0.
        return abs(coupling) * np.sqrt(n_up / sig[0])
    pos = sig > 0
    if not np.any(pos):
        return np.inf
    n_up = poisson_upper_limit(obs[pos], cl)
    return abs(coupling) * float(np.min(np.sqrt(n_up / sig[pos])))


def electron_limits_vs_mass(experiment, masses, set_coefficients, rhoDM=PAPER_HALO['rhoDM'], jchi=0.5,
                            shells=None, m_e=M_ELECTRON, cl=0.9):
    '''
    Limits on a coupling as a function of the DM mass.

    Parameters
    ----------
    experiment : {"Xenon10", "Xenon1T"}
    masses : iterable of float
        DM masses [GeV].
    set_coefficients : callable
        set_coefficients(mDM) sets the coefficients for a unit coupling (after
        reset_coefficients, which is called here), e.g.
        lambda m: set_any_Ecoeff(1., 8) or lambda m: set_any_E_EFT_coeff("anapole", m, 1.).
    rhoDM, jchi, shells, m_e : see difrate_electron.
    cl : float
        Confidence level.

    Returns
    -------
    np.ndarray of the limits on the coupling (np.inf: no signal, np.nan: mass outside the table).
    '''
    from .core import reset_coefficients
    limits = []
    for m in masses:
        reset_coefficients()
        set_coefficients(m)
        sig = electron_signal(experiment, rhoDM=rhoDM, mDM=m, jchi=jchi, shells=shells, m_e=m_e)
        limits.append(electron_limit(sig, coupling=1.0, cl=cl))
    reset_coefficients()
    return np.array(limits)
