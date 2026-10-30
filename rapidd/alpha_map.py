import numpy 

import ctypes
import numpy as np
import os



from .core import _crapidd, base_dir, reset_coefficients, isofromneuc, set_any_Ncoeff, set_any_Ecoeff, read_halo
from .core import mproton, mneutron
from .experiments import counts_bin_LZ

gev_to_pb = 3.893793721e+8 # (hbar c)^2: GeV^-2 -> pb
pb_to_cm2 = 1.e-36 # pb -> cm^2

### SHOULD USE VALUES FROM PARAM CARD ### 

m_u = 2.550000e-03
m_d = 5.040000e-03
m_s = 1.010000e-01
m_c = 1.270000e+00
m_b = 4.700000e+00
m_t = 1.720000e+02


def cn_diracDM_SI(alphas, mchi, q=0.0):


    result = np.zeros(15) # Initialize a list with 31 zeros

    result[0] = (0.063255 * alphas["b_even"] +
        0.0628148 * alphas["c_even"]  +
        0.036 * alphas["d_even"]  +
        2 * alphas["d_odd"]  +
        0.0413 * alphas["s_even"]  +
        0.015 * alphas["u_even"] +
        alphas["u_odd"] 
    )
    #return result* vev**2 /(4*mchi*mneutron)
    return result 


def cp_diracDM_SI(alphas, mchi,  q=0.0):
    result = np.zeros(15) # Initialize a list with 31 zeros

    result[0] = (0.0632493 * alphas["b_even"] +
        0.00507276 * alphas["b_odd"] +
        0.0628148 * alphas["c_even"] -
        0.0125906 * alphas["c_odd"] +
        0.032 * alphas["d_even"] +
        1.0063 * alphas["d_odd"] +
        0.00629531 * alphas["e_odd"] +
        0.00629531 * alphas["mu_odd"] +
        0.0413 * alphas["s_even"] +
        0.00629531 * alphas["s_odd"] +
        0.00629531 * alphas["tau_odd"] +
        0.017 * alphas["u_even"] +
        1.98741 * alphas["u_odd"]
    )
    #return result* vev**2 /(4*mchi*mneutron)
    return result 

def cn_diracDM_SD(alphas, mchi, q=0.0):
    
    result = np.zeros(15) # Initialize a list with 31 zeros
    qsq = q**2
    result[3] = (
        -3.588 * alphas["d_even"]
        + 0.0170931 * alphas["d_odd"]
        + 0.124 * alphas["s_even"]
        + 0.00014071 * alphas["s_odd"]
        + 1.504 * alphas["u_even"]
        - 0.00205568 * alphas["u_odd"]
    )

    result[5] = (
        (0.679556 + 2.58587 * qsq) * alphas["d_even"]
        + (-0.0124848 - 0.685272 * qsq) * alphas["s_even"]
        + (-0.667072 - 1.9006 * qsq) * alphas["u_even"]
        / ((0.0182187 + qsq) * (0.300153 + qsq))
    )

    #return result * vev**2 /(4*mchi*mneutron)
    return result 


def cp_diracDM_SD(alphas, mchi, q=0.0):
    result = np.zeros(15) # Initialize a list with 31 zeros
    qsq = q**2
    result[3] = (
        1.504 * alphas["d_even"]
        - 0.00439169 * alphas["d_odd"]
        + 0.124 * alphas["s_even"]
        + 0.00014071 * alphas["s_odd"]
        - 3.588 * alphas["u_even"]
        + 0.00800105 * alphas["u_odd"]
    )

    result[5] = (
        (-0.667072 - 1.9006 * qsq) * alphas["d_even"]
        + (-0.0124848 - 0.685272 * qsq) * alphas["s_even"]
        + (0.679556 + 2.58587 * qsq) * alphas["u_even"]
    ) / ((0.0182187 + qsq) * (0.300153 + qsq))

    #return result* vev**2 /(4*mchi*mneutron)
    return result

alphasSI = {
    "d_even": 0.0,
    "d_odd": 0.0,
    "u_even": 0.0,
    "u_odd": 0.0,
    "s_even": 0.0,
    "s_odd": 0.0,
    "c_even": 0.0,
    "c_odd": 0.0,
    "b_odd": 0.0, 
    "b_even": 0.0, 
    "e_even": 0.0,
    "e_odd" : 0.0, 
    "mu_odd": 0.0,
    "tau_odd":0.0
}

alphasSD = {
    "d_even": 0.0,
    "d_odd": 0.0,
    "u_even": 0.0,
    "u_odd": 0.0,
    "s_even": 0.0,
    "s_odd": 0.0,
    "c_even": 0.0,
    "c_odd": 0.0,
    "b_odd": 0.0, 
    "b_even": 0.0, 
    "e_even": 0.0,
    "e_odd" : 0.0, 
    "mu_odd": 0.0,
    "tau_odd":0.0
}

def set_alphas_rpd(alphasSI, alphasSD, mchi, q=0.0):
    #print(alphasSI)
    cp = cp_diracDM_SI(alphasSI, mchi, q) + cp_diracDM_SD(alphasSD, mchi, q)
    cn = cn_diracDM_SI(alphasSI, mchi, q) + cn_diracDM_SD(alphasSD, mchi, q)

    c0,c1=isofromneuc(cp,cn)
    
    #print(c0)
    #print(c1)
    for i in range(15):
        set_any_Ncoeff(c0[i], i+1, "p")
        set_any_Ncoeff(c1[i], i+1, "n")
    
    return


###################### 
# Codes for the on-the fly sigma calculation
######################

def gen_Pff_mdm(card):


    Pff = {# the quark -> proton form factor in SIe SIo SDe SDo order
            "alpha_d": np.array([mproton*card['SPd']/m_d, card['VPd'], card['AVPd'], card['SigPd'] ]),
            "alpha_u": np.array([mproton*card['SPu']/m_u, card['VPu'], card['AVPu'], card['SigPu'] ]) ,
            "alpha_s": np.array([mproton*card['SPs']/m_s, 0.0, card['AVPs'], card['SigPs'] ]),
            "alpha_c": np.array([(2/27)*(mproton/m_c)*card['SPg'], 0.0, 0.0, 0.0 ]),##
            "alpha_b": np.array([(2/27)*(mproton/m_b)*card['SPg'], 0.0, 0.0, 0.0 ]),
            "alpha_t": np.array([(2/27)*(mproton/m_t)*card['SPg'], 0.0, 0.0, 0.0 ]) #
    }

    return Pff


def gen_Nff_mdm(card):
    Nff = {# the quark -> neutron form factor in SIe SIo SDe SDo order
            "alpha_d": np.array([mneutron*card['SNd']/m_d, card['VNd'], card['AVNd'], card['SigNd'] ]),
            "alpha_u": np.array([mneutron*card['SNu']/m_u, card['VNu'], card['AVNu'], card['SigNu'] ]) ,
            "alpha_s": np.array([mneutron*card['SNs']/m_s, 0.0, card['AVNs'], card['SigNs'] ]),
            "alpha_c": np.array([(2/27)*(mneutron/m_c)*card['SNg'], 0.0, 0.0, 0.0 ]),##
            "alpha_b": np.array([(2/27)*(mneutron/m_b)*card['SNg'], 0.0, 0.0, 0.0 ]),
            "alpha_t": np.array([(2/27)*(mneutron/m_t)*card['SNg'], 0.0, 0.0, 0.0 ]) # 
    }

    return Nff


QUARK_ALPHAS = ("alpha_d", "alpha_u", "alpha_s", "alpha_c", "alpha_b", "alpha_t")


def _alphas_missing(values):
    '''The Fortran sets all the alphas to -1 when the DM-target process does not exist'''
    return np.all(np.asarray(values, dtype=float) == -1.0)


def nucleon_alphas_valid(alphas):
    '''True if at least one quark has a DM-quark process (alphas != -1)'''
    return any(q in alphas and len(alphas[q]) == 4 and not _alphas_missing(alphas[q]) for q in QUARK_ALPHAS)


def electron_alphas_valid(alphas):
    '''True if the DM-electron process exists (alpha_e != -1)'''
    return "alpha_e" in alphas and len(alphas["alpha_e"]) == 4 and not _alphas_missing(alphas["alpha_e"])


def c1_dirac_mdm(card, alphaq):
    '''Non-relativistic NREFT coefficients c1 (proton, neutron) [GeV^-2] of O1 = 1,
    in RAPIDD's convention (sigma_SI = mu^2 c1^2 / pi).
    The alphas are the coefficients of the EffOperators reference operators, which have
    the Dirac normalisation for both COMPLEX and REAL DM: the scalar (SI even) and vector
    (SI odd) DM-quark couplings give f_N = sum_q alpha FF and chibar chi Nbar N -> c1 = f_N.'''

    PFF = gen_Pff_mdm(card)
    NFF = gen_Nff_mdm(card)

    c1p, c1n = 0.0, 0.0 
    for quark in QUARK_ALPHAS:
        # quarks without a DM-quark process (alphas = -1) contribute 0
        if quark not in alphaq or len(alphaq[quark]) != 4 or _alphas_missing(alphaq[quark]):
            continue
        values = np.asarray(alphaq[quark], dtype=float)
        c1p += (values[0:2] * PFF[quark][0:2]).sum()
        c1n += (values[0:2] * NFF[quark][0:2]).sum()

    return c1p, c1n

def c4_dirac_mdm(card, alphaq):
    '''Non-relativistic NREFT coefficients c4 (proton, neutron) [GeV^-2] of O4 = S_chi.S_N,
    in RAPIDD's convention (sigma_SD = 3 mu^2 c4^2 / (16 pi)).
    The axial (SD even) and tensor (SD odd) alphas, projected on the Dirac-normalised
    axial-vector reference operator, give a_N = sum_q alpha FF, the coefficient of
    (chibar gamma^mu gamma^5 chi)(Nbar gamma_mu gamma^5 N). In the NR limit
    ubar gamma^i gamma^5 u -> 2m xi^dag sigma^i xi, so M -> -16 a_N m_chi m_N S_chi.S_N
    and c4 = -4 a_N.'''
    PFF = gen_Pff_mdm(card)
    NFF = gen_Nff_mdm(card)

    c4p, c4n = 0.0, 0.0 
    for quark in QUARK_ALPHAS:
        # quarks without a DM-quark process (alphas = -1) contribute 0
        if quark not in alphaq or len(alphaq[quark]) != 4 or _alphas_missing(alphaq[quark]):
            continue
        values = np.asarray(alphaq[quark], dtype=float)
        c4p += (values[2:] * PFF[quark][2:]).sum()
        c4n += (values[2:] * NFF[quark][2:]).sum()

    return -4*c4p, -4*c4n




def sigmaSI_nucleon_mdm(card, alphaq, dmmass):
    '''SI DM-proton and DM-neutron cross sections [cm^2]'''
    c1p, c1n = c1_dirac_mdm(card, alphaq)
    #print('c1p = %.2e c1n: %.2e' % (c1p, c1n))
    
    mup = dmmass * mproton / (dmmass + mproton)
    mun = dmmass * mneutron / (dmmass + mneutron)


    sigp = mup**2 * c1p**2 / np.pi * gev_to_pb * pb_to_cm2
    sign = mun**2 * c1n**2 / np.pi * gev_to_pb * pb_to_cm2

    return sigp, sign

def sigmaSD_nucleon_mdm(card, alphaq, dmmass):
    '''SD DM-proton and DM-neutron cross sections [cm^2]'''
    c4p, c4n = c4_dirac_mdm(card, alphaq)
    #print('c4p = %.2e c4n: %.2e' % (c4p, c4n))
    mup = dmmass * mproton / (dmmass + mproton)
    mun = dmmass * mneutron / (dmmass + mneutron)


    sigp = 3 * mup**2 * c4p**2 / (16 * np.pi) * gev_to_pb * pb_to_cm2
    sign = 3 * mun**2 * c4n**2 / (16 * np.pi) * gev_to_pb * pb_to_cm2

    return sigp, sign






###################### 
# DM-electron scattering
######################

def jchi_from_dm_spin(dm_spin):
    '''DM spin from the MadDM dm_spin = 2s+1 (1: scalar, 2: fermion). Vector DM is
    not implemented: returns -1.'''
    dm_spin = int(round(dm_spin))
    if dm_spin == 1:
        return 0.0
    if dm_spin == 2:
        return 0.5
    return -1.0


def ce_dm_electron(alphas, mchi, m_e, dm_spin):
    '''NREFT DM-electron coefficients c1, c4 [GeV units of the legacy MadDM
    get_dm_response] from alphas["alpha_e"] = (aSIe, aSIo, aSDe, aSDo).
    The unit factor gev_to_pb*pb_to_cm2 is applied in the C function dm_response_electron.
    Vector DM is not implemented: returns -1, -1.'''
    aSIe, aSIo, aSDe, aSDo = [float(a) for a in alphas["alpha_e"]]
    dm_spin = int(round(dm_spin))
    if dm_spin == 1:
        # The scalar effective operator already contains the DM mass (2*mchi*phi^dag phi ebar e),
        # so in terms of the alphas the matching is the same as for the fermion SI case
        c1 = 4*mchi*m_e*(aSIo + aSIe)
        c4 = 0.0
    elif dm_spin == 2:
        c1 = 4*mchi*m_e*(aSIo + aSIe)
        # The alphas are projections on the axial x axial operator: the DM-e- amplitude is
        # (aSDe + aSDo) times it, as for c1 (a tensor term gives aSDo = -2 C_T)
        c4 = 16*mchi*m_e*(aSDo + aSDe)
    else:
        return -1.0, -1.0
    return c1, c4


def set_alphas_electron_rpd(alphas, mchi, m_e, dm_spin):
    '''Set the DM-electron coefficients Ce(1), Ce(4) in RAPIDD. Returns False for the
    DM spins not implemented.'''
    if jchi_from_dm_spin(dm_spin) < 0:
        return False
    c1, c4 = ce_dm_electron(alphas, mchi, m_e, dm_spin)
    set_any_Ecoeff(c1, 1)
    set_any_Ecoeff(c4, 4)
    return True




if __name__== '__main__':

    rhochi_p = 0.3
    mchi_p = 50
    alphasSI["u_odd"]=(1e-20)**(1/4)
    alphasSI["d_odd"]=(1e-20)**(1/4)

    set_alphas_rpd(alphasSI, alphasSD, mchi_p) 
    
    read_halo()
  
    result = counts_bin_LZ(rhochi_p, mchi_p, 3.0, 60.0)* (60/1000)*2 *0.9

    print(result)