import numpy 

import ctypes
import numpy as np
import os



from rapidd.core import _crapidd, base_dir, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo
from rapidd.experiments import counts_bin_LZ

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

alphas = {
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

def set_alphas_rpd(alphas, mchi, q=0.0):
    cp = cp_diracDM_SI(alphas, mchi, q) + cp_diracDM_SD(alphas, mchi, q)
    cn = cn_diracDM_SI(alphas, mchi, q) + cn_diracDM_SD(alphas, mchi, q)

    c0,c1=isofromneuc(cp,cn)
    
    #print(c0)
    
    for i in range(15):
        set_any_Ncoeff(c0[i], i+1, "p")
        set_any_Ncoeff(c1[i], i+1, "n")
    
    return


if __name__== '__main__':

    rhochi_p = 0.3
    mchi_p = 50
    alphas["u_odd"]=(1e-20)**(1/4)
    alphas["d_odd"]=(1e-20)**(1/4)

    set_alphas_rpd(alphas, mchi_p) 
    
    read_halo()
  
    result = counts_bin_LZ(rhochi_p, mchi_p, 3.0, 60.0)* (60/1000)*2 *0.9

    print(result)