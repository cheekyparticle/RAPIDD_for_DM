import ctypes
import numpy as np
import os

from rapidd.core import _crapidd
from rapidd.core import base_dir
from scipy.interpolate import interp1d




######## efficiency tables ########

_read_efficiency = _crapidd.read_efficiency
_read_efficiency.argtypes = [ctypes.c_char_p] 
_read_efficiency.restype = ctypes.c_void_p

def read_efficiency(path):
    _read_efficiency(path.encode())
    return 


######## response ######

def lindhard (Er) :
    k=0.166
    Z=54 # so only for xenon
    e=11.5*Er*(Z**(-7/3))
    g=3*e**0.15 + 0.7*e**0.6 + e
    return k*g/(1+k*g)


########## LZ experiment ##########

LZ22_eff_path = os.path.join(base_dir, '..', 'lib', 'build', 'efficiency_tables','LZ_NR_2022.csv')

def read_LZ_eff(path = LZ22_eff_path):
    read_efficiency(path)
    return 



_counts_effres_bin_LZ = _crapidd.counts_effres_bin_LZ

_counts_effres_bin_LZ.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p]
_counts_effres_bin_LZ.restype = ctypes.c_double


def counts_bin_LZ(rhoDM, mDM, E1, E2, effpath = LZ22_eff_path, model="None", basis="ISO"):
    read_efficiency(effpath)
    return _counts_effres_bin_LZ(rhoDM, mDM, E1, E2, model.encode(), basis.encode())

 




######### Xenon1T ############### 

Xe1T_eff_path = os.path.join(base_dir, '..', 'lib', 'build', 'efficiency_tables','Xenon1t.dat')


def read_Xe1T_eff(path = Xe1T_eff_path):
    read_efficiency(path)
    return 






