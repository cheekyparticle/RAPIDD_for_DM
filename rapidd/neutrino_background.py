import numpy as np 

from rapidd.core import base_dir
import os
import csv
from scipy.interpolate import interp1d


neutrino_ar_path = os.path.join(base_dir, '..', 'lib', 'experiments', 'DS20', 'Ar-spectra-DSWP.dat')

DS20_TDReff = os.path.join(base_dir, '..', 'lib', 'experiments', 'DS20', 'DS20k-NREfficiency-TDR.dat')



def get_neutrino_background_DS20K (endroi = 200, nupath = neutrino_ar_path, effpath = DS20_TDReff) :
    exposure_tonneyear = 200
    neutrino_ar = np.loadtxt(nupath)
    energies = np.logspace(-1,3, 200)
    ds20keff = effpath

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