import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d
from scipy.stats import poisson
import csv


def get_simple_limit (coeff,mchi,counts) :
    if counts==0 : 
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        return np.sqrt(coeff* 2.3/counts)


def test_statistic(c,dm,bkgrd,obs) :
    #print(dm, bkgrd, obs)
    return np.sum( -2 * obs * np.log(((dm*c)+bkgrd)/bkgrd) + 2 * (dm * c) )

def binned_poisson_likelihood_limit (coeff,mchi,counts,data, background, cl=2.706) :
    '''cl is the desired confidence level, 0.2707 is for Cl table 40.2 pdg statistics'''
    observed = data
    bkgrd = background
    
    if np.sum(counts)==0:
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        
        c_90 = float(optimize.root(lambda c: test_statistic(c,counts,bkgrd,observed)-cl, 1).x)
        return np.sqrt(c_90 * coeff**2)