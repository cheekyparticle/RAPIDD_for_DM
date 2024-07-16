import numpy 

import ctypes
import numpy as np
import os
from scipy.interpolate import interp1d


from rapidd.core import _crapidd, base_dir, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo
from rapidd.experiments import counts_bin_LZ, lindhard, LZ22_eff_path, counts_bin_Xe1T, counts_bin_DS50, DS50_LEFF
from rapidd.stats import poisson_likelihood_limit, binned_poisson_likelihood_limit, get_simple_limit



########### LZ DATA + Background ##########

LZ22datapath = os.path.join(base_dir, '..', 'lib', 'experiments', 'lz2022','lz2022-data.csv')

LZ22backgroundpath = os.path.join(base_dir, '..', 'lib', 'experiments', 'lz2022','lz2022-bkgrd.csv')


e_kevee22, data22 = np.loadtxt(LZ22datapath,delimiter=',',unpack=True)  ### Need to add where these come from
e_kevee22, bkgrd22 = np.loadtxt(LZ22backgroundpath,delimiter=',',unpack=True)

def lzlimit22(rhoDM, mchi, op=1, fnfp=1., coeff=1e-3, e_kevee=e_kevee22, data=data22, bkgrd=bkgrd22) :
    reset_coefficients()
    
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

    set_any_Ncoeff(c0, op, "p") # ci, i (operator number), p: proton and n:neutron  
    set_any_Ncoeff(c1, op, "n") # ci, i (operator number), p: proton and n:neutron

    ## scale the bkgrds and data down so that there are 11 bkgrd events
    totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))
    
    ## calculate the dm
    limarray = []
    #for mchi in masses: 
    dm=(np.vectorize(counts_bin_LZ)(rhoDM,mchi,E1_lind,E2_lind) * (60/1000) *0.9) 
    lim=(binned_poisson_likelihood_limit(coeff, mchi, dm , totaldata*scaling, totalbkgrd*scaling) )
        
    return lim


def lzlimitproj(rhoDM, mchi, op=1, fnfp=1., coeff=1e-3, eff_file=LZ22_eff_path,
               e_kevee=e_kevee22, bkgrd=bkgrd22) :
    reset_coefficients()
    
    ## import bkgrds and sort out binning in kevee units
    #e_kevee, data = np.loadtxt('lznew/lz2022-data.csv',delimiter=',',unpack=True)
    #e_kevee, bkgrd = np.loadtxt('lznew/lz2022-bkgrd.csv',delimiter=',',unpack=True)
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

    set_any_Ncoeff(c0, op, "p") # ci, i (operator number), p: proton and n:neutron  
    set_any_Ncoeff(c1, op, "n") # ci, i (operator number), p: proton and n:neutron

    ## scale the bkgrds and data down so that there are 11 bkgrd events
    #totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))
    
    ## calculate the dm
    limarray = []
    #for mchi in masses: 
    dm=(np.vectorize(counts_bin_LZ)(rhoDM,mchi,E1_lind,E2_lind,eff_file)  )    #### 
    lim=(binned_poisson_likelihood_limit(coeff, mchi, dm , totalbkgrd*scaling*1000/60, totalbkgrd*scaling*1000/60) ) ### 
        
    return lim





############ Xenon-1T ############
def Xe1TLimits (rhoDM, masses, op=1, fnfp=1., coeff=1e-3) :
    reset_coefficients()
    
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)
    
    set_any_Ncoeff(c0, op, "p") # ci, i (operator number), p: proton and n:neutron  
    set_any_Ncoeff(c1, op, "n") # ci, i (operator number), p: proton and n:neutron
    totalcounts = counts_bin_Xe1T(rhoDM, masses, 1, 100)
    bkgrd = 1.62
    observed = 2.
    xe1tlimit = poisson_likelihood_limit(cp,masses,totalcounts,bkgrd,observed) 

    return xe1tlimit




######### DS50 ################

def DS50Limits_res(rhoDM, mass, op=1, fnfp=1., coeff=1e-3) :
    reset_coefficients()
    
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    set_any_Ncoeff(c0, op, "p") # ci, i (operator number), p: proton and n:neutron  
    set_any_Ncoeff(c1, op, "n") # ci, i (operator number), p: proton and n:neutron
    
    totalcounts = counts_bin_DS50(rhoDM, mass, 40, 200)
    print(totalcounts, mass)
    ds50limit = get_simple_limit(cp,mass,totalcounts)
    return ds50limit




###############################################################################
#                                 MAIN                                       ##
###############################################################################


if __name__== '__main__':
    import matplotlib.pyplot as plt 
    from rapidd.core import calc_xsec_SI, calc_xsec_SD 
    mspace= np.geomspace(1e0,1e3,50)

    rhoDM = 0.3 
    read_halo()

    LZresultSI = np.zeros(np.shape(mspace))
    LZresultSD = np.zeros(np.shape(mspace))
    LZprojSI = np.zeros(np.shape(mspace))
    LZprojSD = np.zeros(np.shape(mspace))

    Xe1TresultSI = np.zeros(np.shape(mspace))
    Xe1TresultSD = np.zeros(np.shape(mspace))


    DS50resSI = np.zeros(np.shape(mspace))
    
    


    for i in range(len(mspace)):
        # LZresultSI[i] = calc_xsec_SI(mspace[i], lzlimit22(rhoDM, mspace[i], op=1))
        # LZresultSD[i] = calc_xsec_SD(mspace[i], lzlimit22(rhoDM, mspace[i], op=4, coeff=1e1))
        
        # LZprojSI[i] = calc_xsec_SI(mspace[i], lzlimitproj(rhoDM, mspace[i], op=1))
        # LZprojSD[i] = calc_xsec_SD(mspace[i], lzlimitproj(rhoDM, mspace[i], op=4, coeff=1e1))
            
        # Xe1TresultSI[i] = calc_xsec_SI(mspace[i], Xe1TLimits(rhoDM, mspace[i], op=1))
        # Xe1TresultSD[i] = calc_xsec_SD(mspace[i], Xe1TLimits(rhoDM, mspace[i], op=4, coeff=1e1))

        DS50resSI[i] = calc_xsec_SI(mspace[i], DS50Limits_res(rhoDM, mspace[i], op=1)) 
        
        
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns

    # ax1.loglog(mspace, LZresultSI, label='LZ-2022')
    # ax1.loglog(mspace, Xe1TresultSI, label='Xenon1T')

    # ax1.loglog(mspace, LZprojSI, ls='--', label='LZ future')

    ax1.loglog(mspace, DS50resSI, label='DS50')

    print(DS50resSI)
    
    
    ax1.set_xlabel(r'$m_{\rm DM}\,\,\left[{\rm GeV}\right]$')
    ax1.set_ylabel(r'$\sigma_{N}^{\rm SI}\,\,\left[{\rm cm}^2\right]$')

    
    
    # ax2.loglog(mspace, LZresultSD)
    
    # ax2.loglog(mspace, Xe1TresultSD)

    # ax2.loglog(mspace, LZprojSD, ls='--')

    ax2.set_xlabel(r'$m_{\rm DM}\,\,\left[{\rm GeV}\right]$')
    ax2.set_ylabel(r'$\sigma_{N}^{\rm SD}\,\,\left[{\rm cm}^2\right]$')

    ax1.legend()


    plt.show()
    