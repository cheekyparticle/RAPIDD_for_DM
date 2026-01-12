from core import base_dir, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo, difrate_dER, vev
import numpy as np

def DDrate_save(card, alphaq, dmmass, output_path, halo_path=base_dir+"/SHM.dat"):
    

    from alpha_map import c1_dirac_mdm, c4_dirac_mdm 
    reset_coefficients()

    rhoDM = card['rhoDM']  
    read_halo(halo_path)
    c1p, c1n = c1_dirac_mdm(card, alphaq)
    c4p, c4n = c4_dirac_mdm(card, alphaq) 

    set_any_Ncoeff(c1p*vev**2, 1, "p")
    set_any_Ncoeff(c1n*vev**2, 1, "n")
    set_any_Ncoeff(c4p*vev**2, 4, "p")
    set_any_Ncoeff(c4n*vev**2, 4, "n")
    
    energies = np.linspace(0.01, 40, 1000)

    drde_Xe = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Xe", basis="All")
    drde_Ar = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ar", basis="All")
    drde_Ge = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ge", basis="All")

    np.savez(output_path+"/DDrates.npz", energies=energies, drde_Xe=drde_Xe, drde_Ar=drde_Ar, drde_Ge=drde_Ge)

    np.savetxt(output_path+"/DDrates.txt", np.column_stack((energies,drde_Xe, drde_Ar, drde_Ge)),
           header="Er [keV] drde_Xe [events/kg/day/keV] drde_Ar [events/kg/day/keV] drde_Ge [events/kg/day/keV]")

    return





    