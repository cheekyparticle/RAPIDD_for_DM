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


vev = 246.2 # GeV
mneutron = 0.939565 # GeV
mproton = 0.938272 # GeV
amu = 0.9315  # GeV

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
halo_path = os.path.join(base_dir, '../lib/halo_table/halo_table_SHM.dat')

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
    '''
    Set the coefficient for a given operator and nucleon.
    
    Parameters
    ----------
    coeff : float
        The coefficient value to set.
    op : int
        The operator index (1-15), 2 excluded.
    nuc : str
        The nucleon type ('p' for proton, 'n' for neutron).
    '''
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

def set_any_N_EFT_coeff(lagr_num, m_chi, m_M, nucl, E, nucleon, c=1.0):
    """\
    Set the EFT coefficients from Table 1 of https://arxiv.org/pdf/1308.6288
    lagr_num: int, the number of the Lagrangian
    E: float, recoil energy in keV
    m_M: float, mediator mass in GeV
    nucl: str, target nucleus (e.g., 'xe', 'ge', 'ar')
    m_chi: float, DM mass in GeV. Required for j in {3,4,5,6,7,8,9,12}
    c: float, coupling to the nucleon (default is 1)
    nucleon: str, the type of nucleon ('p' for proton, 'n' for neutron)
    Returns:
        cNRp: dict, NREFT coefficients for protons
        cNRn: dict, NREFT coefficients for neutrons
    """
    # natural masses of nuclei in amu (atomic mass units)
    nucleus_masses = {'xe': 131.293,
                      'ge': 72.64,
                      'ar': 39.948}
    try:
        m_nucl = nucleus_masses[nucl.lower()]
    except KeyError:
        raise ValueError(f"{nucl} not implemented. Choose from {list(nucleus_masses.keys())}.")

    if nucleon.lower() == 'p':
        m_N = mproton
    elif nucleon.lower() == 'n':
        m_N = mneutron
    else:
        raise ValueError(f"Invalid nucleon type: {nucleon}. Choose 'p' for proton or 'n' for neutron.")

    q = np.sqrt(2 * m_nucl * amu * E * 1e-6)
    q_ratio2 = (q / m_M)**2
    mN_ratio1 = m_N / m_M
    mN_ratio2 = (m_N / m_M)**2

    # set the coefficients based on the Lagrangian number (c = user-defined coupling)
    if lagr_num == 1:
        set_any_Ncoeff(c, 1, nucleon)
    elif lagr_num == 2:
        set_any_Ncoeff(c, 10, nucleon)
    elif lagr_num == 3:
        set_any_Ncoeff(-c * m_N / m_chi, 11, nucleon)
    elif lagr_num == 4:
        set_any_Ncoeff(-c * m_N / m_chi, 6, nucleon)
    elif lagr_num == 5:
        set_any_Ncoeff(c * 4 * m_chi * m_N / m_M**2, 1, nucleon)
    elif lagr_num == 6:
        set_any_Ncoeff(-c * (m_chi / m_N) * q_ratio2, 1, nucleon)
        set_any_Ncoeff(c * 4 * m_chi * m_N / m_M**2, 3, nucleon)
    elif lagr_num == 7:
        set_any_Ncoeff(-c * 4 * m_chi / m_M, 7, nucleon)
    elif lagr_num == 8:
        set_any_Ncoeff(c * 4 * m_chi * m_N / m_M**2, 10, nucleon)
    elif lagr_num == 9:
        set_any_Ncoeff(c * (m_N / m_chi) * q_ratio2, 1, nucleon)
        set_any_Ncoeff(-c * 4 * mN_ratio2, 5, nucleon)
    elif lagr_num == 10:
        set_any_Ncoeff(c * 4 * q_ratio2, 4, nucleon)
        set_any_Ncoeff(-c * 4 * mN_ratio2, 6, nucleon)
    elif lagr_num == 11:
        set_any_Ncoeff(-c * 4 * mN_ratio1, 9, nucleon)
    elif lagr_num == 12:
        set_any_Ncoeff(c * (m_N / m_chi) * q_ratio2, 10, nucleon)
        set_any_Ncoeff(c * 4 * q_ratio2, 12, nucleon)
        set_any_Ncoeff(c * 4 * mN_ratio2, 15, nucleon)
    elif lagr_num == 13:
        set_any_Ncoeff(c * 4 * mN_ratio1, 8, nucleon)
    elif lagr_num == 14:
        set_any_Ncoeff(c * 4 * mN_ratio1, 9, nucleon)
    elif lagr_num == 15:
        set_any_Ncoeff(-c * 4, 4, nucleon)
    elif lagr_num == 16:
        set_any_Ncoeff(c * 4 * mN_ratio1, 13, nucleon)
    elif lagr_num == 17:
        set_any_Ncoeff(-c * 4 * mN_ratio2, 11, nucleon)
    elif lagr_num == 18:
        set_any_Ncoeff(c * q_ratio2, 11, nucleon)
        set_any_Ncoeff(c * 4 * mN_ratio2, 15, nucleon)
    elif lagr_num == 19:
        set_any_Ncoeff(c * 4 * mN_ratio1, 14, nucleon)
    elif lagr_num == 20:
        set_any_Ncoeff(-c * 4 * mN_ratio2, 6, nucleon)
    else:
        raise ValueError(f"Lagrangian number {lagr_num} not recognized (must be 1-20).")


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