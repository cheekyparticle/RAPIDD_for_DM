
import ctypes
import numpy as np
import os
from scipy.interpolate import interp1d




from rapidd.core import _crapidd, base_dir, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo

# mproton = 0.938 # GeV 
# mneutron = 0.939565

mneutron = 0.94  # madDM values
mproton = 0.938


#### BELOW IS MY QUARK VALUES #####
# # Light quarks (GeV)
# m_u = 2.16e-3   # Up
# m_d = 4.67e-3   # Down
# m_s = 93.0e-3   # Strange

# # Heavy quarks (GeV)
# m_c = 1.27      # Charm
# m_b = 4.18      # Bottom
# m_t = 172.76    # Top (pole mass)

### SHOULD USE VALUES FROM PARAM CARD ### 

m_u = 2.550000e-03
m_d = 5.040000e-03
m_s = 1.010000e-01
m_c = 1.270000e+00
m_b = 4.700000e+00
m_t = 1.720000e+02




pff = {# the quark -> proton form factor in SIe SIo SDe SDo order
        "alpha_d": np.array([mproton*0.0191/m_d, 1.0, -0.427, -0.23 ]),
        "alpha_u": np.array([mproton*0.0153/m_u, 2.0, 0.842, 0.84 ]) ,
        "alpha_s": np.array([mproton*0.0447/m_s, 0.0, -0.085, -0.046 ]),
        "alpha_c": np.array([(2/27)*(mproton/m_c)*0.9209, 0.0, 0.0, 0.0 ]),##
        "alpha_b": np.array([(2/27)*(mproton/m_b)*0.9209, 0.0, 0.0, 0.0 ]),
        "alpha_t": np.array([(2/27)*(mproton/m_t)*0.9209, 0.0, 0.0, 0.0 ]) # 0.9021 comes from eq 8 of maddm 2.0
}

nff = {# the quark -> neutron form factor in SIe SIo SDe SDo order
        "alpha_d": np.array([mneutron*0.0273/m_d, 2.0, 0.842, 0.84 ]),
        "alpha_u": np.array([mneutron*0.0110/m_u, 1.0, -0.427, -0.23 ]) ,
        "alpha_s": np.array([mneutron*0.0447/m_s, 0.0, -0.085, -0.046 ]),
        "alpha_c": np.array([(2/27)*(mneutron/m_c)*0.917, 0.0, 0.0, 0.0 ]),##
        "alpha_b": np.array([(2/27)*(mneutron/m_b)*0.917, 0.0, 0.0, 0.0 ]),
        "alpha_t": np.array([(2/27)*(mneutron/m_t)*0.917, 0.0, 0.0, 0.0 ]) # 0.9021 comes from eq 8 of maddm 2.0
}

def read_alpha_values(filename):
    import numpy as np

    # Initialize dictionaries to store the values for each quark
    alpha_values = {
        "alpha_d": [],
        "alpha_u": [],
        "alpha_s": [],
        "alpha_c": [],
        "alpha_b": [],
        "alpha_t": []
    }

    # Read the file and extract alpha values
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split(':')
            
            if parts[0].strip() in alpha_values:
                # Convert the four numerical values into a list of floats
                values = list(map(float, parts[1:5]))
                alpha_values[parts[0].strip()] = np.array(values)  # Store as NumPy array

    return alpha_values

def read_mass(filename):

    with open(filename, 'r') as file:
        for line in file:
            parts = line.split(':')
            #print(parts)
            if parts[0].strip() == 'Wimp_Mass':
                mass = float(parts[1].split()[0])
                # Convert the four numerical values into a list of floats
    return mass

def c1_dirac(filename):

    alpha_data = read_alpha_values(filename)


    c1p, c1n = 0.0, 0.0 
    for quark, values in alpha_data.items():
        #print(values[0:2] * pff[quark][0:2])
        c1p += (values[0:2] * pff[quark][0:2]).sum()
        c1n += (values[0:2] * nff[quark][0:2]).sum()

    return c1p, c1n

def c4_dirac(filename):

    alpha_data = read_alpha_values(filename)


    c4p, c4n = 0.0, 0.0 
    for quark, values in alpha_data.items():
        #print(values[0:2] * pff[quark][0:2])
        c4p += (values[2:] * pff[quark][2:]).sum()
        c4n += (values[2:] * nff[quark][2:]).sum()

    return c4p, c4n




def sigmaSI_nucleon(filename):
    c1p, c1n = c1_dirac(filename)
    print('c1p = %.2e c1n: %.2e' % (c1p, c1n))
    mass = read_mass(filename)
    mup = mass * mproton / (mass + mproton)
    mun = mass * mproton / (mass + mneutron)


    sigp = 4 *mup**2 * c1p**2 / np.pi
    sign = 4 *mun**2 * c1n**2 / np.pi

    print('SI sigma_dmp: %.2e GeV^{-2} : %.2e cm^2' % (sigp, sigp*0.0389e-26))
    print('SI sigma_dmn: %.2e GeV^{-2} : %.2e cm^2' % (sign, sign*0.0389e-26))

    return sigp, sign



def sigmaSD_nucleon(filename):
    c4p, c4n = c4_dirac(filename)
    mass = read_mass(filename)
    mup = mass * mproton / (mass + mproton)
    mun = mass * mproton / (mass + mneutron)


    sigp = 12 *mup**2 * c4p**2 / np.pi
    sign = 12 *mun**2 * c4n**2 / np.pi

    print('SD sigma_dmp: %.2e GeV^{-2} : %.2e cm^2' % (sigp, sigp*0.0389e-26))
    print('SD sigma_dmn: %.2e GeV^{-2} : %.2e cm^2' % (sign, sign*0.0389e-26))

    return sigp, sign



########### LZ DATA + Background ##########

LZ22datapath = os.path.join(base_dir, '..', 'lib', 'experiments', 'lz2022','lz2022-data.csv')

LZ22backgroundpath = os.path.join(base_dir, '..', 'lib', 'experiments', 'lz2022','lz2022-bkgrd.csv')


e_kevee22, data22 = np.loadtxt(LZ22datapath,delimiter=',',unpack=True)  ### Need to add where these come from
e_kevee22, bkgrd22 = np.loadtxt(LZ22backgroundpath,delimiter=',',unpack=True)

def lz22_likelihood(rhoDM, filename, e_kevee=e_kevee22, data=data22, bkgrd=bkgrd22):
    from rapidd.experiments import counts_bin_LZ, lindhard, LZ22_eff_path, read_LZ_eff

    from rapidd.stats import test_statistic

    reset_coefficients()
    read_halo()
    
    vev = 246.2
    
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
    c1p, c1n = c1_dirac(filename)
    c4p, c4n = c4_dirac(filename)


    c10,c11=isofromneuc(c1p,c1n)
    c40,c41=isofromneuc(c4p,c4n)

    mchi = read_mass(filename)

    set_any_Ncoeff(c10*vev**2, 1, "p") # ci, i (operator number), p: proton and n:neutron  
    set_any_Ncoeff(c11*vev**2, 1, "n") # ci, i (operator number), p: proton and n:neutron
    set_any_Ncoeff(c40*vev**2*(-4), 4, "p") # ci, i (operator number), p: proton and n:neutron  
    set_any_Ncoeff(c41*vev**2*(-4), 4, "n") # ci, i (operator number), p: proton and n:neutron
    ## scale the bkgrds and data down so that there are 11 bkgrd events
    totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))

    dm=(np.vectorize(counts_bin_LZ)(rhoDM,mchi,E1_lind,E2_lind) * (60/1000) *0.9)
    print('here %.2e' % dm.sum())
    lik = test_statistic(1.0, dm , totaldata*scaling, totalbkgrd*scaling)
    return lik

if __name__ == "__main__":
    # Example Usage
    filename = "testoutput/maddm_output.out"  # Replace with actual file path
    alpha_data = read_alpha_values(filename)
    mass = read_mass(filename)
    print('dark matter mass %.2e GeV\n' % (mass))
    # Print the results
    for quark, values in alpha_data.items():
        print(f"{quark}: {values}")

    sigmaSI_nucleon(filename)
    sigmaSD_nucleon(filename)
    ll = lz22_likelihood(0.3, filename)
    print(ll)