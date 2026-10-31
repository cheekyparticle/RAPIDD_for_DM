import ctypes

import numpy as np
import os
import platform


# Determine the platform and set the library name accordingly
if platform.system() == 'Linux':
    lib_name = 'libRAPIDD.so'
elif platform.system() == 'Darwin':  # Darwin is macOS
    lib_name = 'libRAPIDD.dylib'
elif platform.system() == 'Windows':
    lib_name = 'RAPIDD.dll'
else:
    raise OSError("Unsupported operating system")

#### Call the shared library #######
base_dir = os.path.dirname(os.path.realpath(__file__))

# Construct the absolute path to the shared library
lib_path = os.path.join(base_dir, '..', 'lib', 'build', lib_name)

_crapidd = ctypes.CDLL(lib_path)


# same values as lib/source/phys_consts.h
G_F = 1.1663788e-5 # Fermi constant [GeV^-2] (PDG 2024)
vev = 1/np.sqrt(np.sqrt(2)*G_F) # Higgs vev [GeV], 246.2196
mneutron = 0.93956542194 # [GeV] (CODATA 2022)
mproton = 0.93827208943 # [GeV] (CODATA 2022)

########################################################################################
# Halo related functions
########################################################################################
_define_and_write_halo = _crapidd.define_and_write_halo
_define_and_write_halo.argtypes = [ctypes.c_char_p, ctypes.c_char_p] + [ctypes.c_double] * 5 \
                                  + [ctypes.c_int] + [ctypes.c_double] * 5
_define_and_write_halo.restype = ctypes.c_void_p
_lab_frame_speed_annual_avg = _crapidd.lab_frame_speed_annual_avg  
_lab_frame_speed_annual_avg.argtypes = [ctypes.c_double, ctypes.c_double * 3]  
_lab_frame_speed_annual_avg.restype = ctypes.c_double  
_define_and_write_halo_time = _crapidd.define_and_write_halo_time  
_define_and_write_halo_time.argtypes = [
    ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double * 3, ctypes.c_double, ctypes.c_int,
    ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
]
  
_define_and_write_halo_time_path = _crapidd.define_and_write_halo_time_path
_define_and_write_halo_time_path.argtypes = [ctypes.c_char_p] + _define_and_write_halo_time.argtypes

VALID_PROFILES = ("SHM", "SHM_numeric", "SHM_beta", "SHM_wLMC", "Lisanti")  
  
def calc_halo(table_path, profile="SHM", vesc=544, v0=238., beta=0.,
               v0_lsr=238., v_pec=(11.1, 12.2, 7.3), k_lisanti=1.5, i=2,
               t0=151., T=None, w=0.006, vb=570., cosb=-0.71, sigb=100., vcut=200., table_dir=None, ve=None):
    """
    Compute the dark matter halo velocity-integral table and write it to disk
    (calls define_and_write_halo if T is None, define_and_write_halo_time
    otherwise).

    Parameters
    ----------
    table_path : str
        Destination path for the output table.
    profile : {"SHM", "SHM_beta", "SHM_wLMC", "Lisanti"}
        Halo velocity-distribution model.
    vesc : float
        Galactic escape velocity [km/s].
    v0 : float
        Velocity-dispersion parameter of the halo distribution [km/s].
    beta : float
        Smooth-cutoff parameter, used only by profile="SHM_beta".
    v0_lsr : float
        Local Standard of Rest circular speed around the Galactic center
        [km/s]. Used only to build the Sun/Earth velocity relative to the
        halo rest frame. It does not affect the shape of the velocity distribution.
    v_pec : tuple of 3 floats
        Solar peculiar velocity (v_r, v_phi, v_theta) relative to the LSR
        [km/s], used together with v0_lsr to compute the Sun's velocity ve.
    k_lisanti : float
        Shape parameter of the "Lisanti" generalized (Tsallis-like) profile.
        Ignored for "SHM"/"SHM_beta". Must be > 0. Typical range,
        following arXiv:1802.03174, is [0.5, 3.5].
    i : int
        i - 1 is the highest power of v included in the tabulated halo integrals:
        eta_i(vmin) = int_{vmin}^infty v^{i-1} f(vec{v}) d^3v.
        i=2 includes all the integrals for the usual differential-rate calculation.
    t0 : float
        Reference day offset (days since March 22, 2018) used only when
        T is not None, to phase the annual modulation.
    T : int or None
        If None, write a single annual-average (static) table. If an
        integer, write T daily-modulated tables (halo_table_0.dat ...
        halo_table_{T-1}.dat), one per day from t=0 to T-1.
    w : float
        Weight fraction (w:[0,1], eg. 0.6%->w=0.006). Only used for SHM_wLMC.
    vb : float
        Bulk velocity of the LMC component. Only used for SHM_wLMC.
    cosb : float
        Cosine of the angle between v_b and v_lab. Only used for SHM_wLMC.
    sigb : float
        Sigma_b parameter of the LMC component. Only used for SHM_wLMC.
    vcut : float
        Max velocity distance from vb. Only used for SHM_wLMC.
    table_dir : str or None
        Only used when T is not None: directory where the daily tables
        halo_table_<t>.dat are written. If None, they are written in
        halo_table/ relative to the current working directory.
    ve : float or None
        Only used when T is None: speed of the Earth in the galactic rest frame [km/s].
        If None (default), it is the annual average computed from v0_lsr and v_pec.
    """
    if profile not in VALID_PROFILES:  
        raise ValueError(f"Profile must be one of {VALID_PROFILES}")

    if profile == "Lisanti" and k_lisanti <= 0:
        raise ValueError(
            "For profile='Lisanti', k must be > 0. Typical range: k in [0.5, 3.5]."
        )
  
    v_pec_c = (ctypes.c_double * 3)(*v_pec)
    if T is None:
        if ve is None:
            ve = _lab_frame_speed_annual_avg(v0_lsr, v_pec_c)
        _define_and_write_halo(
            table_path.encode(), profile.encode(),
            vesc, v0, beta, ve, k_lisanti, i, w,
            vb, cosb, sigb, vcut
        )
    elif table_dir is None:
        _define_and_write_halo_time(
            profile.encode(), vesc, v0, beta,
            v0_lsr, v_pec_c, k_lisanti, i, t0, T,
            w, vb, cosb, sigb, vcut
        )
    else:
        os.makedirs(table_dir, exist_ok=True)
        _define_and_write_halo_time_path(
            os.path.abspath(table_dir).encode(), profile.encode(), vesc, v0, beta,
            v0_lsr, v_pec_c, k_lisanti, i, t0, T,
            w, vb, cosb, sigb, vcut
        )

_read_halo = _crapidd.read_halo
_read_halo.argtypes = [ctypes.c_char_p]
_read_halo.restype = ctypes.c_void_p 
halo_path = os.path.join(base_dir, '../lib/halo_table/halo_table_SHM.dat')

def read_halo(path = halo_path):
    _read_halo(path.encode())
    return

# halo parameters of the default table halo_path (calc_halo defaults)
default_halo_vesc = 544.
default_halo_v0 = 238.
default_halo_v_pec = (11.1, 12.2, 7.3)

def halo_path_from_card(card, table_dir):
    """
    Annual-average SHM table for the halo parameters of the MadDM card.

    Returns halo_path if vescape, vmp and (vSun_r, vSun_phi, vSun_theta) are the
    default ones, otherwise computes the table in table_dir/SHM.dat and returns
    its path. vmp is used both as v0 and as the LSR speed.
    """
    card_pars = (card['vescape'], card['vmp'], card['vSun_r'], card['vSun_phi'], card['vSun_theta'])
    if np.allclose(card_pars, (default_halo_vesc, default_halo_v0) + default_halo_v_pec, rtol=1e-9, atol=0.):
        return halo_path
    path = os.path.join(table_dir, 'SHM.dat')
    calc_halo(path, profile="SHM", vesc=card['vescape'], v0=card['vmp'], v0_lsr=card['vmp'],
              v_pec=(card['vSun_r'], card['vSun_phi'], card['vSun_theta']))
    return path

########################################################################################
# Functions used to set coefficient
########################################################################################
_set_any_coeffs = _crapidd.set_any_coeffs
_set_any_coeffs.argtypes = [ctypes.c_double, ctypes.c_int]

def set_any_coeffs(C, i):
    return _set_any_coeffs(C, i)

_set_any_Ncoeff = _crapidd.set_any_Ncoeff
_set_any_Ncoeff.argtypes= [ ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
_set_any_Ncoeff.restype = ctypes.c_void_p

_Cp = _crapidd.Cp
_Cp.argtypes = [ctypes.c_int]
_Cp.restype = ctypes.c_double

_Cn = _crapidd.Cn
_Cn.argtypes = [ctypes.c_int]
_Cn.restype = ctypes.c_double

_set_any_Ecoeff = _crapidd.set_any_Ecoeff
_set_any_Ecoeff.argtypes = [ctypes.c_double, ctypes.c_int]
_set_any_Ecoeff.restype = None

_Ce = _crapidd.Ce
_Ce.argtypes = [ctypes.c_int]
_Ce.restype = ctypes.c_double

_set_any_Ecoeff_LR = _crapidd.set_any_Ecoeff_LR
_set_any_Ecoeff_LR.argtypes = [ctypes.c_double, ctypes.c_int]
_set_any_Ecoeff_LR.restype = None

_Ce_LR = _crapidd.Ce_LR
_Ce_LR.argtypes = [ctypes.c_int]
_Ce_LR.restype = ctypes.c_double

def set_any_Ncoeff(coeff, op, nuc):
    '''coeff, op, nuc'''
    _set_any_Ncoeff(coeff, op, nuc.encode())
    return 

ELECTRON_OPERATORS = (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
ELECTRON_INTERACTIONS = ("contact", "long-range")

def set_any_Ecoeff(coeff, op, interaction="contact"):
    '''
    Set the DM-electron NREFT coefficient of an operator of arXiv:1912.08204.

    The amplitude is M = sum_i (c_i^s + (q_ref/q)^2 c_i^l) <O_i>, with q_ref = alpha m_e
    and dimensionless c_i (e.g. c_1 = 4 m_chi m_e a_SI for a contact SI coupling a_SI
    [GeV^-2]). Note that the limits on c_i^l depend on the choice of q_ref.

    Parameters
    ----------
    coeff : float
        The dimensionless coefficient c_i^s (contact) or c_i^l (long-range).
    op : int
        The operator index: 1, 3-15 (O_2, O_16 and O_17-O_20 do not enter the DM responses
        of arXiv:1912.08204 and are ignored by the rate).
    interaction : {"contact", "long-range"}
        Contact (c_i^s) or long-range (c_i^l) coefficient.
    '''
    if interaction == "contact":
        _set_any_Ecoeff(coeff, op)
    elif interaction == "long-range":
        _set_any_Ecoeff_LR(coeff, op)
    else:
        raise ValueError(f"interaction must be one of {ELECTRON_INTERACTIONS}")
    return

def reset_coefficients() :
    '''Set to zero all the DM-nucleon and DM-electron (contact and long-range) coefficients'''
    for i in range(16):
        set_any_Ncoeff(0, i, "p")
        set_any_Ncoeff(0, i, "n")    
        set_any_Ecoeff(0, i)
        set_any_Ecoeff(0, i, interaction="long-range")
    return

def Cp_val(op):
    return _Cp(op)

def Cn_val(op):
    return _Cn(op)

def Ce_val(op):
    '''Contact DM-electron coefficient c_op^s'''
    return _Ce(op)

def Ce_LR_val(op):
    '''Long-range DM-electron coefficient c_op^l'''
    return _Ce_LR(op)

def isofromneuc(cp, cn):
    ''' Simply takes coeffs from p n basis to 0 1 basis '''
    return (cp+cn)/2, (cp-cn)/2

########################################################################################
# Differential rate functions
########################################################################################
_difrate_dER_python = _crapidd.difrate_dER_python
_difrate_dER_python.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double]
_difrate_dER_python.restype = ctypes.c_double

_nuclear_model_supports_target = _crapidd.nuclear_model_supports_target
_nuclear_model_supports_target.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
_nuclear_model_supports_target.restype = ctypes.c_int

VALID_NUCLEAR_MODELS = ("iso_GCN5082", "pn_Fitz", "iso_Fitz", "pn_SN100PN")
_supported_targets = {}

def check_nuclear_model(target, basis):
    '''Raise a ValueError if basis is not a valid nuclear model or if its form factors are
    not tabulated for all the isotopes of target (difrate_dER would silently return 0)'''
    if basis not in VALID_NUCLEAR_MODELS:
        raise ValueError(f"Nuclear model must be one of {VALID_NUCLEAR_MODELS}")
    if (target, basis) not in _supported_targets:
        _supported_targets[(target, basis)] = bool(_nuclear_model_supports_target(basis.encode(), target.encode()))
    if not _supported_targets[(target, basis)]:
        raise ValueError(f"The nuclear model {basis} has no form factors for (all the isotopes of) the target {target}")

@np.vectorize
def difrate_dER(rhoDM, mDM, ER, model="None", target="Xe", basis="iso_GCN5082", delta=0.0):
    check_nuclear_model(target, basis)
    return _difrate_dER_python(rhoDM, mDM, np.log10(ER), model.encode(), target.encode(), basis.encode(), delta)

########################################################################################
# DM-electron scattering
########################################################################################
GRIDSIZE_K = 100
GRIDSIZE_Q = 200
N_SHELL_XE = 11
XE_SHELLS = ("1s", "2p", "2s", "3d", "3p", "3s", "4d", "4p", "4s", "5p", "5s")
S2_MAX_XENON10 = 203
N_BINS_XENON10 = 7
S2_BIN_MAX_XENON1T = 63
N_ATOMIC_RESP = 4
M_ELECTRON = 0.51099895000e-3 # electron mass [GeV] (CODATA 2018)

atomic_resp_dir = os.path.join(base_dir, '..', 'lib', 'atomic_responses')
experiments_dir = os.path.join(base_dir, '..', 'lib', 'experiments')

_read_atomic_responses = _crapidd.read_atomic_responses
_read_atomic_responses.argtypes = [ctypes.c_char_p, ctypes.c_int]
_read_atomic_responses.restype = ctypes.c_int

def read_atomic_responses(path=atomic_resp_dir, response_index=None):
    '''
    Read the Xenon atomic responses W_j(k', q) of arXiv:1912.08204 from
    path/atomic_response_<j>/Xe_<shell>.txt.

    Parameters
    ----------
    path : str
        Directory containing atomic_response_1 ... atomic_response_4.
    response_index : int, iterable of int or None
        Response(s) j to read (1-4). None reads all the four responses (W_1 is the only
        one needed by O_1, O_4 and the other operators that do not depend on v_el^perp).
    '''
    if response_index is None:
        indices = range(1, N_ATOMIC_RESP + 1)
    elif np.ndim(response_index) == 0:
        indices = [int(response_index)]
    else:
        indices = [int(j) for j in response_index]
    for j in indices:
        resp_dir = os.path.join(path, 'atomic_response_%i' % j)
        err = _read_atomic_responses(resp_dir.encode(), j)
        if err != 0:
            raise IOError("Could not read the atomic responses in %s (error %i)" % (resp_dir, err))
    return

_dm_response_electron = _crapidd.dm_response_electron
_dbl_p = ctypes.POINTER(ctypes.c_double)
_int_p = ctypes.POINTER(ctypes.c_int)
_dm_response_electron.argtypes = [ctypes.c_int] + [ctypes.c_double] * 5 + [_dbl_p]
_dm_response_electron.restype = ctypes.c_double

def dm_response_electron(mchi, jchi, resp=1, q=0.0, dE=0.0, m_e=M_ELECTRON, return_v2=False):
    '''
    DM response R_resp of arXiv:1912.08204 built from the electron coefficients Ce
    (contact) and Ce_LR (long-range), times the legacy MadDM normalisation
    gev_to_pb*pb_to_cm2 [GeV^-2 -> cm^2].

    Parameters
    ----------
    mchi : float
        DM mass [GeV].
    jchi : float
        DM spin (0 or 1/2).
    resp : int
        Response index (1-4).
    q : float
        Exchanged momentum [GeV]. Must be > 0 if a long-range coefficient is set.
    dE : float
        Deposited energy E_B + E_e [GeV] (enters through v_el^perp at k = 0).
    m_e : float
        Electron mass [GeV].
    return_v2 : bool
        If True, also return the coefficient of v^2 (v in units of c), not zero only for R_1.

    Returns
    -------
    float or (float, float)
        The v independent part of R_resp (and the coefficient of v^2 if return_v2).
        With q = dE = 0 and only c_1, c_4 contact, R_1 = c_1^2 + jchi(jchi+1)/4 c_4^2 as in MadDM.
    '''
    R_v2 = ctypes.c_double(0.0)
    R = _dm_response_electron(resp, q, dE, mchi, jchi, m_e, ctypes.byref(R_v2))
    if return_v2:
        return R, R_v2.value
    return R

_difrate_electron_dEr_all = _crapidd.difrate_electron_dEr_all
_difrate_electron_dEr_all.argtypes = [ctypes.c_double] * 4 + [_dbl_p] * 3
_difrate_electron_dEr_all.restype = ctypes.c_int

_set_electron_w2_sign = _crapidd.set_electron_w2_sign
_set_electron_w2_sign.argtypes = [ctypes.c_int]
_set_electron_w2_sign.restype = None
_electron_w2_sign_value = _crapidd.electron_w2_sign_value
_electron_w2_sign_value.argtypes = []
_electron_w2_sign_value.restype = ctypes.c_int

W2_CONVENTIONS = {"Liang": -1, "Catena": 1}

def set_w2_sign(convention="Liang"):
    '''
    Sign of the tabulated atomic response W_2 in the DM-electron rate.

    Parameters
    ----------
    convention : {"Liang", "Catena"} or {-1, +1}
        "Liang" (-1, default): sign of the vector atomic form factor corrected in
        arXiv:2405.04855. "Catena" (+1): W_2 as in arXiv:1912.08204 (to reproduce its
        figures). Only the operators O_7, O_8, O_12, O_13, O_14 (R_2 != 0) are affected.
    '''
    sign = W2_CONVENTIONS.get(convention, convention)
    if sign not in (1, -1):
        raise ValueError(f"convention must be one of {list(W2_CONVENTIONS)} or +1/-1")
    _set_electron_w2_sign(int(sign))

def w2_sign():
    '''Sign of W_2 currently used in the DM-electron rate (-1: Liang, +1: Catena)'''
    return _electron_w2_sign_value()

_difrate_electron_dEr_resp = _crapidd.difrate_electron_dEr_resp
_difrate_electron_dEr_resp.argtypes = [ctypes.c_double] * 4 + [ctypes.c_int] + [_dbl_p] * 3
_difrate_electron_dEr_resp.restype = ctypes.c_int

_electron_rate_errors = {1: "a needed atomic response W_j has not been read (read_atomic_responses)",
                         2: "the halo table has no halo(vmin, 2) column, needed by the operators "
                            "depending on v_el^perp (write it with i >= 2)",
                         3: "wrong resp_mask"}

def difrate_electron_dEr_all(rhoDM, mDM, jchi, m_e, resp_mask=15):
    '''Returns k [GeV] (GRIDSIZE_K+1), E [GeV] (GRIDSIZE_K+1) and the per-shell
    dR/dE [events/kg/day/GeV] (N_SHELL_XE x GRIDSIZE_K) for DM-electron scattering in Xenon.
    resp_mask selects the terms R_j W_j (bit j-1 set: included; 15: all).'''
    k = np.zeros(GRIDSIZE_K + 1)
    E = np.zeros(GRIDSIZE_K + 1)
    rate = np.zeros(N_SHELL_XE * GRIDSIZE_K)
    if resp_mask == 15:
        err = _difrate_electron_dEr_all(rhoDM, mDM, jchi, m_e, k.ctypes.data_as(_dbl_p),
                                        E.ctypes.data_as(_dbl_p), rate.ctypes.data_as(_dbl_p))
    else:
        err = _difrate_electron_dEr_resp(rhoDM, mDM, jchi, m_e, int(resp_mask), k.ctypes.data_as(_dbl_p),
                                         E.ctypes.data_as(_dbl_p), rate.ctypes.data_as(_dbl_p))
    if err != 0:
        raise RuntimeError("difrate_electron_dEr_all failed (error %i): %s"
                           % (err, _electron_rate_errors.get(err, "unknown error")))
    return k, E, rate.reshape(N_SHELL_XE, GRIDSIZE_K)

_xenon10_bins = _crapidd.xenon10_bins
_xenon10_bins.argtypes = [_int_p]
_xenon10_bins.restype = ctypes.c_int

def xenon10_bins():
    bins = np.zeros(N_BINS_XENON10 + 1, dtype=np.int32)
    _xenon10_bins(bins.ctypes.data_as(_int_p))
    return bins

_difrate_dS2 = _crapidd.difrate_dS2
_difrate_dS2.argtypes = [ctypes.c_char_p, _dbl_p, _dbl_p, ctypes.c_double, ctypes.c_char_p,
                         _dbl_p, _dbl_p, _dbl_p, _dbl_p, _dbl_p, _int_p]
_difrate_dS2.restype = ctypes.c_int

def difrate_dS2(experiment, E, rate_shell, mDM, exp_dir=experiments_dir):
    '''Detector response in S2 for experiment "Xenon10" or "Xenon1T", given the grids E and
    the per-shell dR/dE from difrate_electron_dEr_all. Returns a dict with S2, dRdS2, signal,
    n_obs (and dRdS2_shell for Xenon10, n_bkg for Xenon1T).'''
    E = np.ascontiguousarray(E, dtype=np.float64)
    rate_shell = np.ascontiguousarray(rate_shell, dtype=np.float64).ravel()
    if experiment == "Xenon10":
        nS2, nsig = S2_MAX_XENON10, N_BINS_XENON10
        shell = np.zeros(N_SHELL_XE * nS2)
        shell_p = shell.ctypes.data_as(_dbl_p)
    elif experiment == "Xenon1T":
        nS2, nsig = S2_BIN_MAX_XENON1T, 1
        shell, shell_p = None, None
    else:
        raise ValueError("experiment must be Xenon10 or Xenon1T")
    S2 = np.zeros(nS2)
    dRdS2 = np.zeros(nS2)
    signal = np.zeros(nsig)
    n_bkg = np.zeros(1)
    n_obs = np.zeros(nsig, dtype=np.int32)
    n = _difrate_dS2(experiment.encode(), E.ctypes.data_as(_dbl_p), rate_shell.ctypes.data_as(_dbl_p),
                     mDM, exp_dir.encode(), S2.ctypes.data_as(_dbl_p), dRdS2.ctypes.data_as(_dbl_p),
                     shell_p, signal.ctypes.data_as(_dbl_p), n_bkg.ctypes.data_as(_dbl_p),
                     n_obs.ctypes.data_as(_int_p))
    if n < 0:
        raise RuntimeError("difrate_dS2 failed for %s" % experiment)
    res = {'S2': S2, 'dRdS2': dRdS2, 'signal': signal, 'n_obs': n_obs}
    if experiment == "Xenon10":
        res['dRdS2_shell'] = shell.reshape(N_SHELL_XE, nS2)
    else:
        res['signal'] = signal[0]
        res['n_obs'] = int(n_obs[0])
        res['n_bkg'] = n_bkg[0]
    return res
