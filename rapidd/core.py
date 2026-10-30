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


vev = 246.2
mneutron = 0.939565
mproton = 0.938272

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
               t0=151., T=None, w=0.006, vb=570., cosb=-0.71, sigb=100., vcut=200., table_dir=None):
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
    """
    if profile not in VALID_PROFILES:  
        raise ValueError(f"Profile must be one of {VALID_PROFILES}")

    if profile == "Lisanti" and k_lisanti <= 0:
        raise ValueError(
            "For profile='Lisanti', k must be > 0. Typical range: k in [0.5, 3.5]."
        )
  
    v_pec_c = (ctypes.c_double * 3)(*v_pec)
    if T is None:
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

def set_any_Ncoeff(coeff, op, nuc):
    '''coeff, op, nuc'''
    _set_any_Ncoeff(coeff, op, nuc.encode())
    return 

def set_any_Ecoeff(coeff, op):
    '''Set the DM-electron NREFT coefficient of operator op'''
    _set_any_Ecoeff(coeff, op)
    return

def reset_coefficients() :
    for i in range(16):
        set_any_Ncoeff(0, i, "p")
        set_any_Ncoeff(0, i, "n")    
        set_any_Ecoeff(0, i)
    return

def Cp_val(op):
    return _Cp(op)

def Cn_val(op):
    return _Cn(op)

def Ce_val(op):
    return _Ce(op)

def isofromneuc(cp, cn):
    ''' Simply takes coeffs from p n basis to 0 1 basis '''
    return (cp+cn)/2, (cp-cn)/2

########################################################################################
# Differential rate functions
########################################################################################
_difrate_dER_python = _crapidd.difrate_dER_python
_difrate_dER_python.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double]
_difrate_dER_python.restype = ctypes.c_double

VALID_NUCLEAR_MODELS = ("iso_GCN5082", "pn_Fitz", "iso_Fitz", "pn_SN100PN")

@np.vectorize
def difrate_dER(rhoDM, mDM, ER, model="None", target="Xe", basis="iso_GCN5082", delta=0.0):
    if basis not in VALID_NUCLEAR_MODELS:  
        raise ValueError(f"Nuclear model must be one of {VALID_NUCLEAR_MODELS}")
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

atomic_resp_dir = os.path.join(base_dir, '..', 'lib', 'atomic_responses')
experiments_dir = os.path.join(base_dir, '..', 'lib', 'experiments')

_read_atomic_responses = _crapidd.read_atomic_responses
_read_atomic_responses.argtypes = [ctypes.c_char_p, ctypes.c_int]
_read_atomic_responses.restype = ctypes.c_int

def read_atomic_responses(path=atomic_resp_dir, response_index=1):
    '''Read the Xenon atomic response W_<response_index> from
    path/atomic_response_<response_index>/Xe_<shell>.txt'''
    resp_dir = os.path.join(path, 'atomic_response_%i' % response_index)
    err = _read_atomic_responses(resp_dir.encode(), response_index)
    if err != 0:
        raise IOError("Could not read the atomic responses in %s (error %i)" % (resp_dir, err))
    return

_dm_response_electron = _crapidd.dm_response_electron
_dm_response_electron.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double]
_dm_response_electron.restype = ctypes.c_double

def dm_response_electron(mchi, jchi, resp=1, q=0.0):
    '''DM response R_resp built from the electron coefficients Ce (legacy MadDM normalisation)'''
    return _dm_response_electron(resp, q, mchi, jchi)

_dbl_p = ctypes.POINTER(ctypes.c_double)
_int_p = ctypes.POINTER(ctypes.c_int)

_difrate_electron_dEr_all = _crapidd.difrate_electron_dEr_all
_difrate_electron_dEr_all.argtypes = [ctypes.c_double] * 4 + [_dbl_p] * 3
_difrate_electron_dEr_all.restype = ctypes.c_int

def difrate_electron_dEr_all(rhoDM, mDM, jchi, m_e):
    '''Returns k [GeV] (GRIDSIZE_K+1), E [GeV] (GRIDSIZE_K+1) and the per-shell
    dR/dE [events/kg/day/GeV] (N_SHELL_XE x GRIDSIZE_K) for DM-electron scattering in Xenon'''
    k = np.zeros(GRIDSIZE_K + 1)
    E = np.zeros(GRIDSIZE_K + 1)
    rate = np.zeros(N_SHELL_XE * GRIDSIZE_K)
    err = _difrate_electron_dEr_all(rhoDM, mDM, jchi, m_e, k.ctypes.data_as(_dbl_p),
                                    E.ctypes.data_as(_dbl_p), rate.ctypes.data_as(_dbl_p))
    if err != 0:
        raise RuntimeError("difrate_electron_dEr_all failed (error %i)" % err)
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
