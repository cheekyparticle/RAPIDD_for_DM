from .core import halo_path, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo, difrate_dER, vev
import numpy as np

def DDrate_save(card, alphaq, dmmass, output_path, halo_path=halo_path):
    

    from .alpha_map import c1_dirac_mdm, c4_dirac_mdm 
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

    drde_Xe = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Xe", basis="iso_GCN5082")
    drde_Ar = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ar", basis="pn_Fitz")
    drde_Ge = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ge", basis="pn_Fitz")

    np.savez(output_path+"/DDrates.npz", energies=energies, drde_Xe=drde_Xe, drde_Ar=drde_Ar, drde_Ge=drde_Ge)

    np.savetxt(output_path+"/DDrates.txt", np.column_stack((energies,drde_Xe, drde_Ar, drde_Ge)),
           header="Er [keV] drde_Xe [events/kg/day/keV] drde_Ar [events/kg/day/keV] drde_Ge [events/kg/day/keV]")

    return





    ########################################################################################
# DM-electron scattering
########################################################################################

def _f77(*items):
    '''Format a record like the gfortran list-directed write (write(unit,*) ...), so that
    the DM-electron output files keep the format of the legacy MadDM Fortran.'''
    out = ' '
    first = True
    for it in items:
        if isinstance(it, str):
            out += it
        else:
            if not first:
                out += ' '
            if isinstance(it, (int, np.integer)):
                out += '%11d' % it
            else:
                x = float(it)
                ax = abs(x)
                if x == 0.0:
                    out += '%20s     ' % '0.0000000000000000'
                elif 0.1 <= ax < 1e16:
                    ndig = len('%d' % int(ax)) if ax >= 1 else 0
                    out += '%20s     ' % ('%.*f' % (17 - ndig, x))
                else:
                    mant, exp = ('%.16E' % x).split('E')
                    out += '%25s' % ('%sE%s%03d' % (mant, exp[0], abs(int(exp))))
        first = False
    return out + '\n'


def _electron_rates(rhoDM, dmmass, jchi, m_e):
    from .core import difrate_electron_dEr_all
    return difrate_electron_dEr_all(rhoDM, dmmass, jchi, m_e)


def DDrate_electron_save(card, alphas, dm_spin, m_e, dmmass, output_path, halo_path=halo_path,
                         modulation=False, modulation_halo_dir=None, days=None):
    '''Compute and save the DM-electron scattering rates in Xenon and the XENON10 / XENON1T
    S2 signals (port of the legacy MadDM electron_recoil_signal.f).

    Parameters
    ----------
    card : maddm_card (rhoDM, and vmp, vescape, vSun_* for the modulation tables)
    alphas : dict with "alpha_e" = (aSIe, aSIo, aSDe, aSDo)
    dm_spin : MadDM dm_spin = 2s+1 (1: scalar, 2: fermion)
    m_e : electron mass [GeV]
    dmmass : DM mass [GeV]
    output_path : directory where the files are written
    halo_path : annual-average halo table (eta_0 = halo(vmin, 0))
    modulation : if True, also compute the rate for each day of the year
    modulation_halo_dir : directory with the daily tables halo_table_<t>.dat. If None they are
        generated in output_path/halo_table_modulation from the card (t = days since March 22).
    days : days for the modulation (default range(365))

    Returns a dict with DM_response, sigma_e, Xenon10_bins, Xenon10_signal, Xenon10_obs,
    Xenon1T_signal, Xenon1T_bkg, Xenon1T_obs (DM_response = -1 if the DM spin is not implemented).
    '''
    import os
    from .core import (read_atomic_responses, dm_response_electron, difrate_dS2, xenon10_bins,
                       XE_SHELLS, calc_halo)
    from .alpha_map import set_alphas_electron_rpd, jchi_from_dm_spin

    result = {'DM_response': -1}
    jchi = jchi_from_dm_spin(dm_spin)
    reset_coefficients()
    if not set_alphas_electron_rpd(alphas, dmmass, m_e, dm_spin):
        return result

    rhoDM = card['rhoDM']
    read_halo(halo_path)
    read_atomic_responses()

    dm_resp = dm_response_electron(dmmass, jchi)
    result['DM_response'] = dm_resp
    result['sigma_e'] = dm_resp/(16.*np.pi*(dmmass + m_e)**2)

    k, E, rate_shell = _electron_rates(rhoDM, dmmass, jchi, m_e)
    rate_tot = rate_shell.sum(axis=0)

    x10 = difrate_dS2("Xenon10", E, rate_shell, dmmass)
    x1t = difrate_dS2("Xenon1T", E, rate_shell, dmmass)
    bins = xenon10_bins()

    result['Xenon10_bins'] = [int(b) for b in bins]
    result['Xenon10_signal'] = list(x10['signal'])
    result['Xenon10_obs'] = [int(n) for n in x10['n_obs']]
    result['Xenon1T_signal'] = x1t['signal']
    result['Xenon1T_bkg'] = x1t['n_bkg']
    result['Xenon1T_obs'] = x1t['n_obs']

    with open(os.path.join(output_path, 'signal_e_recoil.dat'), 'w') as f:
        f.write(_f77('Xenon10_bins', *[int(b) for b in bins]))
        f.write(_f77('Xenon10_signal', *x10['signal']))

    with open(os.path.join(output_path, 'dRdE_e_recoil.dat'), 'w') as f:
        f.write(_f77('##    Momentum (keV)        Energy(keV)               dR/dE[events/kg/day/keV]'))
        for ik in range(len(rate_tot)):
            f.write(_f77(k[ik]*1e6, E[ik]*1e6, rate_tot[ik]/1e6))

    with open(os.path.join(output_path, 'dRdS2_Xenon10_e_recoil.dat'), 'w') as f:
        f.write(_f77('##       S2   dN/dS2[events]'))
        for S2, val in zip(x10['S2'], x10['dRdS2']):
            f.write(_f77(int(S2), val))

    with open(os.path.join(output_path, 'dRdS2_Xenon1T_e_recoil.dat'), 'w') as f:
        f.write(_f77('##                S2        dN/dS2[events]'))
        for S2, val in zip(x1t['S2'], x1t['dRdS2']):
            f.write(_f77(S2, val))

    npz = dict(k_keV=k[:-1]*1e6, E_keV=E[:-1]*1e6, drde_tot=rate_tot/1e6,
               drde_shell=rate_shell/1e6, shells=np.array(XE_SHELLS),
               S2_Xenon10=x10['S2'], dRdS2_Xenon10=x10['dRdS2'], dRdS2_shell_Xenon10=x10['dRdS2_shell'],
               S2_Xenon1T=x1t['S2'], dRdS2_Xenon1T=x1t['dRdS2'])

    if modulation:
        if days is None:
            days = range(365)
        days = list(days)
        if modulation_halo_dir is None:
            modulation_halo_dir = os.path.join(output_path, 'halo_table_modulation')
            calc_halo(None, profile="SHM", vesc=card['vescape'], v0=card['vmp'], v0_lsr=card['vmp'],
                      v_pec=(card['vSun_r'], card['vSun_phi'], card['vSun_theta']),
                      t0=0., T=max(days) + 1, table_dir=modulation_halo_dir)
        drde_day = np.zeros((len(days), len(rate_tot)))
        rate_day = np.zeros(len(days))
        dE = np.diff(E)
        for i, t in enumerate(days):
            read_halo(os.path.join(modulation_halo_dir, 'halo_table_%i.dat' % t))
            _, _, rs = _electron_rates(rhoDM, dmmass, jchi, m_e)
            drde_day[i] = rs.sum(axis=0)
            rate_day[i] = (drde_day[i]*dE).sum()   # events/kg/day
        read_halo(halo_path)
        with open(os.path.join(output_path, 'rate_vs_time_e_recoil.dat'), 'w') as f:
            f.write('## Days since 22-03-2018    Rate[events/kg/day]\n')
            for t, r in zip(days, rate_day):
                f.write('%6i  %.10E\n' % (t, r))
        npz.update(days=np.array(days), drde_day=drde_day/1e6, rate_day=rate_day)

    np.savez(os.path.join(output_path, 'DDrates_e.npz'), **npz)
    np.savetxt(os.path.join(output_path, 'DDrates_e.txt'),
               np.column_stack((k[:-1]*1e6, E[:-1]*1e6, rate_tot/1e6, (rate_shell/1e6).T)),
               header="k [keV] E [keV] drde_tot [events/kg/day/keV] " +
                      " ".join("drde_%s" % s for s in XE_SHELLS))

    return result
