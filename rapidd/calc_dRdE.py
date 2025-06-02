from core import base_dir, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo, difrate_dER
import numpy as np

def DDrate_save(card, alphaq, dmmass, output_path, halo_path=base_dir+"/SHM.dat"):
    

    from alpha_map import c1_dirac_mdm, c4_dirac_mdm 
    reset_coefficients()

    rhoDM = card['rhoDM']  
    read_halo(halo_path)
    c1p, c1n = c1_dirac_mdm(card, alphaq)
    c4p, c4n = c4_dirac_mdm(card, alphaq) 

    set_any_Ncoeff(c1p, 1, "p")
    set_any_Ncoeff(c1n, 1, "n")
    set_any_Ncoeff(c4p, 4, "p")
    set_any_Ncoeff(c4n, 4, "n")
    
    energies = np.linspace(0.01, 40, 100)

    drde_Xe = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Xe", basis="All")
    drde_Ar = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ar", basis="All")
    drde_Ge = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ge", basis="All")

    np.savez(output_path+"/DDrates.npz", energies=energies, drde_Xe=drde_Xe, drde_Ar=drde_Ar, drde_Ge=drde_Ge)

    np.savetxt(output_path+"/DDrates.txt", np.column_stack((energies,drde_Xe, drde_Ar, drde_Ge)),
           header="Er drde_Xe drde_Ar drde_Ge")

    return





    