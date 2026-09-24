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
  
VALID_PROFILES = ("SHM", "SHM_numeric", "SHM_beta", "SHM_wLMC", "Lisanti")  
  
def calc_halo(table_path, profile="SHM", vesc=544, v0=238., beta=0.,
               v0_lsr=238., v_pec=(11.1, 12.2, 7.3), k_lisanti=1.5, i=2,
               t0=151., T=None, w=0.006, vb=570., cosb=-0.71, sigb=100., vcut=200.):
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
    else:
        _define_and_write_halo_time(
            profile.encode(), vesc, v0, beta,
            v0_lsr, v_pec_c, k_lisanti, i, t0, T,
            w, vb, cosb, sigb, vcut
        )

_read_halo = _crapidd.read_halo
_read_halo.argtypes = [ctypes.c_char_p]
_read_halo.restype = ctypes.c_void_p 
halo_path = os.path.join(base_dir, '../lib/halo_table/halo_table.dat')

def read_halo(path = halo_path):
    _read_halo(path.encode())
    return

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

def set_any_Ncoeff(coeff, op, nuc):
    '''coeff, op, nuc'''
    _set_any_Ncoeff(coeff, op, nuc.encode())
    return 

def reset_coefficients() :
    for i in range(16):
        set_any_Ncoeff(0, i, "p")
        set_any_Ncoeff(0, i, "n")    
    return

def Cp_val(op):
    return _Cp(op)

def Cn_val(op):
    return _Cn(op)

def isofromneuc(cp, cn):
    ''' Simply takes coeffs from p n basis to 0 1 basis '''
    return (cp+cn)/2, (cp-cn)/2

########################################################################################
# Differential rate functions
########################################################################################
_difrate_dER_python = _crapidd.difrate_dER_python
_difrate_dER_python.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double]
_difrate_dER_python.restype = ctypes.c_double

VALID_NUCLEAR_MODELS = ("iso_GCN5082", "pn_BD", "pn_SN100PN")

@np.vectorize
def difrate_dER(rhoDM, mDM, ER, model="None", target="Xe", basis="iso_GCN5082", delta=0.0):
    if basis not in VALID_NUCLEAR_MODELS:  
        raise ValueError(f"Nuclear model must be one of {VALID_NUCLEAR_MODELS}")
    return _difrate_dER_python(rhoDM, mDM, np.log10(ER), model.encode(), target.encode(), basis.encode(), delta)