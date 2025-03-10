import numpy as np
import pycoeffs_eft as rapidd
from scipy import optimize
from scipy.interpolate import interp1d
from scipy.stats import poisson
import csv
rapidd.py_read_halo(b"lib/halo_table.dat") ### read the halo function 
pi = 3.1415


def calc_xsec(mchi,coeff) :
    mN = 0.938
    mv = 246
    pi = 3.14159
    u  = (mchi*mN)/(mchi+mN)
    return (float(u**2 * coeff**2 * 0.0389e-26/ (pi * mv**4)))

def reset_coefficients() :
    for i in range(16):
        rapidd.py_set_any_Ncoeff(0, i, b"p")
        rapidd.py_set_any_Ncoeff(0, i, b"n")
        
def isofromneuc(cp, cn):
    ''' Simply takes coeffs from p n basis to 0 1 basis '''
    return (cp+cn)/2, (cp-cn)/2

def lindhard (Er) :
    k=0.166
    Z=54
    e=11.5*Er*(Z**(-7/3))
    g=3*e**0.15 + 0.7*e**0.6 + e
    return k*g/(1+k*g)

def test_statistic(c,dm,bkgrd,obs) :
    #print(dm, bkgrd, obs)
    return np.sum( -2 * obs * np.log(((dm*c)+bkgrd)/bkgrd) + 2 * (dm * c) )

def binned_poisson_likelihood_limit (coeff,mchi,counts,data, background) :
    
    observed = data
    bkgrd = background
    
    if np.sum(counts)==0:
        return np.inf
    else:
        crosssec = calc_xsec(mchi,coeff)
        
        c_90 = float(optimize.root(lambda c: test_statistic(c,counts,bkgrd,observed)-2.706, 1).x)
        return c_90 * crosssec
    
def binned_poisson_likelihood_limit_glam (glam4,mchi,counts,data, background) :
    
    observed = data
    bkgrd = background
    
    if np.sum(counts)==0:
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        
        c_90 = float(optimize.root(lambda c: test_statistic(c,counts,bkgrd,observed)-2.706, 1).x)
        return c_90 * glam4
    
def get_simple_limit (coeff,mchi,counts) :
    if counts==0 : 
        return np.inf
    else:
        crosssec = calc_xsec(mchi,coeff)
        return crosssec* 2.3/counts
    
def get_simple_limit_glam (glam4,mchi,counts) :
    if counts==0 : 
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        #return glam* (2.3/counts)**(1/4)
        return glam4* (2.3/counts) 
    
def get_1count_limit (coeff,mchi,counts) :
    if counts==0 : 
        return np.inf
    else:
        crosssec = calc_xsec(mchi,coeff)
        return crosssec* 1/counts
    
def poisson_likelihood_limit (coeff,mchi,counts,bkgrd,observed) :
    if counts==0:
        return np.inf
    else:
        crosssec = calc_xsec(mchi,coeff)
        
        N_90 = float(optimize.root(lambda mu: poisson.cdf(observed, mu) - 0.1, 5).x)-bkgrd
        return crosssec* (N_90/counts)
    
def poisson_likelihood_limit_glam (glam4,mchi,counts,bkgrd,observed) :
    if counts==0:
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        
        N_90 = float(optimize.root(lambda mu: poisson.cdf(observed, mu) - 0.1, 5).x)-bkgrd

        return glam4* (N_90/counts) #**(1/4)

@np.vectorize
def counts_LZ_efficiency(rho, mass, E1, E2):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/LZ_NR_2022.csv') ### read efficiency table 
    return rapidd.py_bin_LZ( rho, mass, E1, E2, b"None" , b"ISO" )


@np.vectorize
def counts_LZ_efficiency_res(rho, mass, E1, E2, eff_file = b'lib/efficiency_tables/LZ_NR_2022.csv'):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/LZ_NR_2022.csv') ### read efficiency table 
    return rapidd.py_bin_LZ_res( rho, mass, E1, E2, b"None" , b"ISO" )

@np.vectorize
def counts_DS50_efficiency(rho, mass, E1, E2):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/DS50.dat') ### read efficiency table 
    return rapidd.py_bin_DS50( rho, mass, E1, E2, b"None" , b"ISO" )

@np.vectorize
def counts_deap_efficiency(rho, mass, E1, E2):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/DEAP3600.dat') ### read efficiency table 
    rapidd.py_read_DS50_LEFF(b'lib/efficiency_tables/LeffDS50.dat')

    return rapidd.py_bin_DS50_res( rho, mass, E1, E2, b"None" , b"ISO" )

@np.vectorize
def counts_DS50_efficiency_res(rho, mass, E1, E2):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/DS50.dat') ### read efficiency table 
    rapidd.py_read_DS50_LEFF(b'lib/efficiency_tables/LeffDS50.dat')
    return rapidd.py_bin_DS50_res( rho, mass, E1, E2, b"None" , b"ISO" )


@np.vectorize
def counts_Xe1T_efficiency(rho, mass, E1, E2):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/Xenon1t.dat') ### read efficiency table 
    return rapidd.py_bin_Xe1t( rho, mass, E1, E2, b"None" , b"ISO" )

@np.vectorize
def counts_DS20k_efficiency(rho, mass, E1, E2):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/DS20k_NR.dat') ### read efficiency table 
    return rapidd.py_bin_DS20k( rho, mass, E1, E2, b"None" , b"ISO" )

@np.vectorize
def counts_DS20k_efficiency_res(rho, mass, E1, E2,eff_file='lib/efficiency_tables/DS20k_NR.dat'):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/DS20k_NR.dat') ### read efficiency table 
    rapidd.py_read_DS50_LEFF(b'lib/efficiency_tables/LeffDS50.dat')
    return rapidd.py_bin_DS20k_res( rho, mass, E1, E2, b"None" , b"ISO" )

@np.vectorize
def counts_DS20k_noefficiency(rho, mass, E1, E2):
    rapidd.py_read_efficiency(b'lib/efficiency_tables/DS20k_NRdoesntexist.dat') ### read efficiency table 
    #rapidd.py_read_DS50_LEFF(b'lib/efficiency_tables/LeffDS50.dat')
    return rapidd.py_bin_DS20k( rho, mass, E1, E2, b"None" , b"ISO" )


def lzlimit(mchi, op=1, fnfp=1., coeff=1e-3) :
    reset_coefficients()
    
    ## import bkgrds and sort out binning in kevee units
    e_kevee, data = np.loadtxt('lznew/lz2022-data.csv',delimiter=',',unpack=True)
    e_kevee, bkgrd = np.loadtxt('lznew/lz2022-bkgrd.csv',delimiter=',',unpack=True)
    E1 = np.linspace(0,17,52)[0:51]
    E2 = np.linspace(0,17,52)[1:]
    bins = (E1+E2)/2
    spacing = bins[1]-bins[0]
    
    ## want to convert binning to kevnr to calculate the dm
    x=np.linspace(0,200,1000)
    y=lindhard(x)

    f_lind = interp1d(x*y,y)
    E1_lind = E1/f_lind(E1)
    E2_lind = E2/f_lind(E2)
    E1_lind[0]=0
    
    ## set the coefficients for the correct operator
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    ## scale the bkgrds and data down so that there are 11 bkgrd events
    totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))
    
    ## calculate the dm
    limarray = []
    #for mchi in masses: 
    dm=(counts_LZ_efficiency_res(0.3,mchi,E1_lind,E2_lind) * (60/1000)*2 *0.9) 
    lim=(binned_poisson_likelihood_limit(coeff, mchi, dm , totaldata*scaling, totalbkgrd*scaling) )
        
    return lim

def lzlimitproj(mchi, op=1, fnfp=1., coeff=1e-3, eff_file=b'lib/efficiency_tables/LZ_NR_2022.csv') :
    reset_coefficients()
    
    ## import bkgrds and sort out binning in kevee units
    #e_kevee, data = np.loadtxt('lznew/lz2022-data.csv',delimiter=',',unpack=True)
    e_kevee, bkgrd = np.loadtxt('lznew/lz2022-bkgrd.csv',delimiter=',',unpack=True)
    E1 = np.linspace(0,17,52)[0:51]
    E2 = np.linspace(0,17,52)[1:]
    bins = (E1+E2)/2
    spacing = bins[1]-bins[0]
    
    ## want to convert binning to kevnr to calculate the dm
    x=np.linspace(0,200,1000)
    y=lindhard(x)

    f_lind = interp1d(x*y,y)
    E1_lind = E1/f_lind(E1)
    E2_lind = E2/f_lind(E2)
    E1_lind[0]=0
    
    ## set the coefficients for the correct operator
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    ## scale the bkgrds and data down so that there are 11 bkgrd events
    #totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))
    
    ## calculate the dm
    limarray = []
    #for mchi in masses: 
    dm=(counts_LZ_efficiency_res(0.3,mchi,E1_lind,E2_lind,eff_file) * 2 ) #(60/1000)*2 ) 
    lim=(binned_poisson_likelihood_limit(coeff, mchi, dm , totalbkgrd*scaling*1000/60, totalbkgrd*scaling*1000/60) )
        
    return lim
        

def DS50Limits_res (masses, op=1, fnfp=1., coeff=1e-3) :
    reset_coefficients()
    
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron
    
    totalcounts = np.array([counts_DS50_efficiency_res(0.3, m, 40, 200) for m in masses])
    ds50limit = np.vectorize(get_simple_limit)(cp,masses,totalcounts)
    return ds50limit

def DS50Limits_res_roi(masses, op=1, fnfp=1., coeff=1e-3, start=40,end=200) :
    reset_coefficients()
    
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron
    
    totalcounts = np.array([counts_DS50_efficiency_res(0.3, m, start, end) for m in masses])
    ds50limit = np.vectorize(get_simple_limit)(cp,masses,totalcounts)
    return ds50limit

def deapLimits_res (masses, op=1, fnfp=1., coeff=1e-3) :
    reset_coefficients()
    
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)
    #print(cp)
    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron
    
    totalcounts = np.array([counts_deap_efficiency(0.3, m, 50, 100) for m in masses])
    #print(totalcounts)
    totalcounts2 = totalcounts*(758*1e3/16660) *0.25 # deap exposure compared to ds50
    #print(totalcounts)
    ds50limit = np.vectorize(get_simple_limit)(cp,masses,totalcounts2)
    return ds50limit

## Xenon-1T ##
def Xe1TLimits (masses, op=1, fnfp=1., coeff=1e-3) :
    reset_coefficients()
    
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)
    
    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron
    totalcounts = np.array([counts_Xe1T_efficiency(0.3, m, 1, 100) for m in masses])
    bkgrd = 1.62
    observed = 2.
    xe1tlimit = (np.vectorize(poisson_likelihood_limit)(cp,masses,totalcounts,bkgrd,observed) )

    return xe1tlimit

def get_neutrino_background (endroi = 200) :
    exposure_tonneyear = 200
    neutrino_ar = np.loadtxt('./Ar-spectra-DSWP.dat')
    energies = np.logspace(-1,3, 200)
    ds20keff = './DS20k-NREfficiency-TDR.dat'

    ## read in the efficiency curve ##
    en, eff = [], []
    with open(ds20keff) as csvfile:
        f = csv.reader(csvfile,delimiter=' ')
        for row in f:
            en.append(float(row[0]))
            eff.append(float(row[1]))
    e = interp1d(en, eff)

    ## split into 1keV bins, apply efficiency and scale ##
    bins = np.linspace(20,endroi-1,endroi-20)+0.5
    g = interp1d(energies,neutrino_ar)
    counts = []

    for i in bins :
        counts.append( g(i) * e(i) * exposure_tonneyear)
        
    return np.array(counts)


## Darkside-20k ##
def DS20kLimits_res_masses (masses, op=1, fnfp=1., coeff=1e-3,eff_file='lib/efficiency_tables/DS20k_NR.dat') :
    reset_coefficients()
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    E1, E2 = np.linspace(20.0, 199.0, 180), np.linspace(21.0, 200.0, 180)
    Width = E2 - E1
    ds20klimit = []
    for m in masses :
        counts = ( np.vectorize(counts_DS20k_efficiency_res)(0.3,m, E1, E2,eff_file) ) 
        ds20klimit.append(binned_poisson_likelihood_limit(cp,m,counts,get_neutrino_background(),get_neutrino_background()) )
    return ds20klimit

## Darkside-20k ##
def DS20kLimit_res (mchi, op=1, fnfp=1., coeff=1e-3,eff_file='lib/efficiency_tables/DS20k_NR.dat') :
    reset_coefficients()
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    E1, E2 = np.linspace(20.0, 199.0, 180), np.linspace(21.0, 200.0, 180)
    Width = E2 - E1
    #ds20klimit = []
    counts = ( np.vectorize(counts_DS20k_efficiency_res)(0.3,mchi, E1, E2,eff_file) ) 
    ds20klimit = (binned_poisson_likelihood_limit(cp,mchi,counts,get_neutrino_background(),get_neutrino_background()) )
    return ds20klimit

## Darkside-20k ##
def DS20kLimit_noeff (mchi, op=1, fnfp=1., coeff=1e-3) :
    reset_coefficients()
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    E1, E2 = np.linspace(20.0, 199.0, 180), np.linspace(21.0, 200.0, 180)
    Width = E2 - E1
    #ds20klimit = []
    counts = ( np.vectorize(counts_DS20k_noefficiency)(0.3,mchi, E1, E2) ) 
    ds20klimit = (binned_poisson_likelihood_limit(cp,mchi,counts,get_neutrino_background(),get_neutrino_background()) )
    return ds20klimit

## Darkside-20k ##
def DS20kLimit_noeff_roi (mchi, op=1, fnfp=1., coeff=1e-3, roiend = 200) :
    reset_coefficients()
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    E1, E2 = np.linspace(20.0, roiend-1, roiend-20), np.linspace(21.0, roiend, roiend-20)
    Width = E2 - E1
    #ds20klimit = []
    counts = ( np.vectorize(counts_DS20k_noefficiency)(0.3,mchi, E1, E2) ) 
    ds20klimit = (binned_poisson_likelihood_limit(cp,mchi,counts,get_neutrino_background(roiend),get_neutrino_background(roiend)) )
    return ds20klimit



## Darkside-20k ##
def DS20kLimits_res_1bin (masses, op=1, fnfp=1., coeff=1e-3) :
    reset_coefficients()
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.py_set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    bkgrd = 3.2 #np.sum ( get_neutrino_background() )
    observed = bkgrd
        
    counts = np.array([counts_DS20k_efficiency(0.3,m, 1, 200) for m in masses])
    ds20klimit = (np.vectorize(poisson_likelihood_limit)(cp,masses,counts,bkgrd,observed) )

    return ds20klimit



def C_6_1(glam, theta):
    
    return glam**2 * np.cos(theta)

def C_7_1(glam, theta, phi):
    return glam**2*np.sin(theta)*np.cos(phi)

def C_8_1(glam, theta, phi):
    return glam**2 * np.sin(theta)*np.sin(phi)

def CAV_u ( glam, theta, phi):
    return C_7_1(glam, theta, phi) + C_6_1(glam, theta)

def CAV_d (glam, theta, phi): 
    return C_8_1(glam, theta, phi) + C_6_1(glam, theta)

def CAA_u ( glam, theta, phi):
    return C_7_1(glam, theta, phi) - C_6_1(glam, theta)

def CAA_d (glam, theta, phi): 
    return C_8_1(glam, theta, phi) - C_6_1(glam, theta)

def set_gaugeinv_AX_coeffs_ISO(glam, theta, phi):
    reset_coefficients()
    vev = 246.2

    ### Axial vector hadron factors for protons, use p -> n and u -> d to obtain the neutron values
    F_up_A = 0.89665 
    F_dp_A = -0.37565
    F_sp_A = -0.031
    
    cp4 = 4* vev**2 * (CAA_d(glam, theta, phi) * F_dp_A + CAA_u(glam, theta, phi) *F_up_A )
    cn4 = 4* vev**2 * (CAA_d(glam, theta, phi) * F_up_A + CAA_u(glam, theta, phi) *F_dp_A )
    
    
    
    #### Op 6 stuff 
    mn = 0.93956542052
    mp = 0.93827208816
    mpi = 0.1396
    meta = 0.547862
    F_up_pp = mp**2/(mpi**2) * 2*(F_up_A-F_dp_A)  + mp**2/(meta**2) * 2*(F_up_A + F_dp_A - 2*F_sp_A)/3
    F_dp_pp = mp**2/(mpi**2) * 2*(-1)*(F_up_A-F_dp_A)  + mp**2/(meta**2) * 2*(F_up_A + F_dp_A - 2*F_sp_A)/3
    
    cp6 = vev**2 * (CAA_d(glam, theta, phi) * F_dp_pp + CAA_u(glam, theta, phi) * F_up_pp )
    cn6 = vev**2 * (CAA_d(glam, theta, phi) * F_up_pp + CAA_u(glam, theta, phi) * F_dp_pp )
    
    
    ### Vector hadron factors for protons, use p -> n and u -> d to obtain the neutron values
    
    F_up_1 = 2
    F_dp_1 = 1 
    
    cp8 = 2* vev**2 * ( CAV_d(glam, theta, phi) * F_dp_1 + CAV_u(glam, theta, phi) * F_up_1)
    cn8 = 2* vev**2 * ( CAV_d(glam, theta, phi) * F_up_1 + CAV_u(glam, theta, phi) * F_dp_1)
    
    
    ### Pauli Form factors 
    
    F_up_2 = 1.609
    F_dp_2 = -2.097

    cp9 = 2* vev**2 * (CAV_d(glam, theta, phi) * ( F_dp_1 + F_dp_2) +
                       CAV_u(glam, theta, phi) * ( F_up_1 + F_up_2))
    cn9 = 2*vev**2 * (CAV_d(glam, theta, phi) * ( F_up_1 + F_up_2) +
                       CAV_u(glam, theta, phi) * ( F_dp_1 + F_dp_2))
    
    c04 = 0.5*(cp4 + cn4)
    c14 = 0.5*(cp4-cn4)
    
    c06 = 0.5*(cp6 + cn6)
    c16 = 0.5*(cp6-cn6)
    
    c08 = 0.5*(cp8 + cn8)
    c18 = 0.5*(cp8-cn8)
    
    c09 = 0.5*(cp9+cn9)
    c19 = 0.5*(cp9-cn9)
    
    #print(c04,c14,c08,c18,c09,c19)

    
    rapidd.py_set_any_Ncoeff(c04, 4, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c14, 4, b"n") # ci, i (operator number), p: proton and n:neutron
    
    rapidd.py_set_any_Ncoeff(c06, 6, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c16, 6, b"n") # ci, i (operator number), p: proton and n:neutron
    
    rapidd.py_set_any_Ncoeff(c08, 8, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c18, 8, b"n") # ci, i (operator number), p: proton and n:neutron
    
    rapidd.py_set_any_Ncoeff(c09, 9, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c19, 9, b"n") # ci, i (operator number), p: proton and n:neutron
    
    return

def C_2_1(glam, theta):
    
    return glam**2 * np.cos(theta)

def C_3_1(glam, theta, phi):
    return glam**2*np.sin(theta)*np.cos(phi)

def C_4_1(glam, theta, phi):
    return glam**2 * np.sin(theta)*np.sin(phi)
    
def C_1_u ( glam, theta, phi):
    return (C_2_1(glam, theta) + C_3_1(glam, theta, phi) ) * 0.5

def C_1_d (glam, theta, phi): 
    return 0.5*(C_2_1(glam, theta) + C_4_1(glam, theta, phi))

def C_3_u ( glam, theta, phi):
    return 0.5*(C_3_1(glam, theta, phi) - C_2_1(glam, theta))

def C_3_d (glam, theta, phi): 
    return (C_4_1(glam, theta, phi) - C_2_1(glam, theta))*0.5
    
def reset_coeffs():
    for i in range(16): 
        rapidd.py_set_any_coeffs(0.0,i)
        

    
def set_gaugeinv_coeffs_ISO(glam, theta, phi, mchi):
    reset_coeffs()
    vev = 246.2
    
    mn = 0.93956542052
    mp = 0.93827208816
    
    #### vector hadron form factors 
    
    F_up_1 = 2
    F_dp_1 = 1 

    cp1 = vev**2 * ( C_1_d(glam, theta, phi) * F_dp_1 + C_1_u(glam, theta, phi) * F_up_1)
    cn1 = vev**2 * ( C_1_d(glam, theta, phi) * F_up_1 + C_1_u(glam, theta, phi) * F_dp_1)
    
    
    ### Axial vector hadron factors for protons, use p -> n and u -> d to obtain the neutron values
    F_up_A = 0.89665 
    F_dp_A = -0.37565
    F_sp_A = -0.031
    
    cp7 = -2* vev**2 * (C_3_d(glam, theta, phi) * F_dp_A + C_3_u(glam, theta, phi) *F_up_A )
    cn7 = -2* vev**2 * (C_3_d(glam, theta, phi) * F_up_A + C_3_u(glam, theta, phi) *F_dp_A )
    
    cp9 = -2*(mp/mchi)*vev**2*(C_3_d(glam, theta, phi) * F_dp_A + C_3_u(glam, theta, phi) *F_up_A)
    cn9 = -2*(mn/mchi)*vev**2*(C_3_d(glam, theta, phi) * F_up_A + C_3_u(glam, theta, phi) *F_dp_A)
    
    #### change to isospin basis
    
    c01 = 0.5*(cp1 + cn1)
    c11 = 0.5*(cp1-cn1)
    
    c07 = 0.5*(cp7 + cn7)
    c17 = 0.5*(cp7-cn7)
    
    c09 = 0.5*(cp9 + cn9)
    c19 = 0.5*(cp9-cn9)
    
    
    rapidd.py_set_any_Ncoeff(c01, 1, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c11, 1, b"n") # ci, i (operator number), p: proton and n:neutron
    
    rapidd.py_set_any_Ncoeff(c07, 7, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c17, 7, b"n") # ci, i (operator number), p: proton and n:neutron
    
    rapidd.py_set_any_Ncoeff(c09, 9, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.py_set_any_Ncoeff(c19, 9, b"n") # ci, i (operator number), p: proton and n:neutron
    
    
    return cn1/cp1


def lzlimit_isovec(mchi, op=1, glam4=1e-20,theta=0., phi=0, dd=False) :
    pi = 3.1415
    reset_coefficients()
    glam=glam4**(1/4)
    
    ## import bkgrds and sort out binning in kevee units
    e_kevee, data = np.loadtxt('lznew/lz2022-data.csv',delimiter=',',unpack=True)
    e_kevee, bkgrd = np.loadtxt('lznew/lz2022-bkgrd.csv',delimiter=',',unpack=True)
    E1 = np.linspace(0,17,52)[0:51]
    E2 = np.linspace(0,17,52)[1:]
    bins = (E1+E2)/2
    spacing = bins[1]-bins[0]
    
    ## want to convert binning to kevnr to calculate the dm
    x=np.linspace(0,200,1000)
    y=lindhard(x)

    f_lind = interp1d(x*y,y)
    E1_lind = E1/f_lind(E1)
    E2_lind = E2/f_lind(E2)
    E1_lind[0]=0
    
    ## set the coefficients for the correct operator
    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)


    ## scale the bkgrds and data down so that there are 11 bkgrd events
    totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))
    
    ## calculate the dm
    limarray = []
    #for mchi in masses: 
    dm=(counts_LZ_efficiency_res(0.3,mchi,E1_lind,E2_lind) * (60/1000)*2 ) 
    lim=(binned_poisson_likelihood_limit_glam(glam4, mchi, dm , totaldata*scaling, totalbkgrd*scaling) )
        
    return lim

def lzlimitproj_isovec(mchi, op=1, glam4=1e-3,theta=0., phi=0, dd=False) :
    
    pi = 3.1415
    reset_coefficients()
    glam=glam4**(1/4)

    ## import bkgrds and sort out binning in kevee units
    #e_kevee, data = np.loadtxt('lznew/lz2022-data.csv',delimiter=',',unpack=True)
    e_kevee, bkgrd = np.loadtxt('lznew/lz2022-bkgrd.csv',delimiter=',',unpack=True)
    E1 = np.linspace(0,17,52)[0:51]
    E2 = np.linspace(0,17,52)[1:]
    bins = (E1+E2)/2
    spacing = bins[1]-bins[0]
    
    ## want to convert binning to kevnr to calculate the dm
    x=np.linspace(0,200,1000)
    y=lindhard(x)

    f_lind = interp1d(x*y,y)
    E1_lind = E1/f_lind(E1)
    E2_lind = E2/f_lind(E2)
    E1_lind[0]=0
    
    ## set the coefficients
    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)


    ## scale the bkgrds and data down so that there are 11 bkgrd events
    #totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))
    
    ## calculate the dm
    limarray = []
    #for mchi in masses: 
    dm=(counts_LZ_efficiency_res(0.3,mchi,E1_lind,E2_lind) * 2 ) #(60/1000)*2 ) 
    lim=(binned_poisson_likelihood_limit_glam(glam4, mchi, dm , totalbkgrd*scaling*1000/60, totalbkgrd*scaling*1000/60) )
        
    return lim

## Darkside-20k ##
def DS20kLimit_res_isovec (mchi, op=1, glam4=1e-3,theta=0., phi=0,dd=False) :
    reset_coefficients()
    glam=glam4**(1/4)

    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)
    

    E1, E2 = np.linspace(20.0, 199.0, 180), np.linspace(21.0, 200.0, 180)
    Width = E2 - E1
    #ds20klimit = []
    counts = ( np.vectorize(counts_DS20k_efficiency_res)(0.3,mchi, E1, E2) ) 
    ds20klimit = (binned_poisson_likelihood_limit_glam(glam4,mchi,counts,get_neutrino_background(),get_neutrino_background()) )
    return ds20klimit

def DS50Limits_res_isovec (masses, op=1, glam4=1e-3,theta=0., phi=0, dd=False) :
    reset_coefficients()
    glam=glam4**(1/4)
    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)

    
    totalcounts = np.array([counts_DS50_efficiency_res(0.3, m, 40, 200) for m in masses])
    ds50limit = np.vectorize(get_simple_limit_glam)(glam4,masses,totalcounts)
    return ds50limit

def deapLimits_res_isovec (masses, op=1, glam4=1e-3,theta=0., phi=0, dd=False) :
    reset_coefficients()
    glam=glam4**(1/4)
    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)

    totalcounts = np.array([counts_deap_efficiency(0.3, m, 50, 100) for m in masses])
    totalcounts2 = totalcounts*(758*1e3/16660) *0.25
    
    ds50limit = np.vectorize(get_simple_limit_glam)(glam4,masses,totalcounts2)
    return ds50limit

## Xenon-1T ##
def Xe1TLimits_isovec (masses, op=1, glam4=1e-3,theta=0., phi=0, dd=False) :
    reset_coefficients()
    glam=glam4**(1/4)
    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)

    totalcounts = np.array([counts_Xe1T_efficiency(0.3, m, 1, 100) for m in masses])
    bkgrd = 1.62
    observed = 2.
    xe1tlimit = (np.vectorize(poisson_likelihood_limit_glam)(glam4,masses,totalcounts,bkgrd,observed) )

    return xe1tlimit

def DS50Limits_res_isovec_1mass (mchi, op=1, glam4=1e-3,theta=0., phi=0, dd=False) :
    reset_coefficients()
    glam=glam4**(1/4)
    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)

    
    totalcounts = counts_DS50_efficiency_res(0.3, mchi, 40, 200) 
    ds50limit = get_simple_limit_glam(glam4,mchi,totalcounts)
    return ds50limit

def deapLimits_res_isovec_1mass (mchi, op=1, glam4=1e-3,theta=0., phi=0, dd=False) :
    reset_coefficients()
    glam=glam4**(1/4)
    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)

    totalcounts = counts_deap_efficiency(0.3, mchi, 50, 100) 
    totalcounts2 = totalcounts*(758*1e3/16660) *0.25
    
    ds50limit = get_simple_limit_glam(glam4,mchi,totalcounts2)
    return ds50limit

## Xenon-1T ##
def Xe1TLimits_isovec_1mass (mchi, op=1, glam4=1e-3,theta=0., phi=0, dd=False) :
    reset_coefficients()
    glam=glam4**(1/4)
    if dd: 
        set_gaugeinv_coeffs_ISO(glam, theta*pi, phi*2*pi, mchi)
    else:
        set_gaugeinv_AX_coeffs_ISO(glam, theta*pi, phi*2*pi)

    totalcounts = (counts_Xe1T_efficiency(0.3, mchi, 1, 100) )
    bkgrd = 1.62
    observed = 2.
    xe1tlimit = poisson_likelihood_limit_glam(glam4,mchi,totalcounts,bkgrd,observed) 

    return xe1tlimit