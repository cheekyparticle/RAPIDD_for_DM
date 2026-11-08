from .core import halo_path, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo, difrate_dER, vev
import numpy as np

NUCLEAR_TARGETS = (("Xe", "iso_GCN5082"), ("Ar", "pn_Fitz"), ("Ge", "pn_Fitz"))


def _nuclear_rates(rhoDM, dmmass, energies):
    '''dR/dE [events/kg/day/keV] of Xe, Ar, Ge with the halo table currently loaded'''
    return {target: difrate_dER(rhoDM, dmmass, energies, target=target, basis=basis)
            for target, basis in NUCLEAR_TARGETS}


def DDrate_save(card, alphaq, dmmass, output_path, halo_path=halo_path,
                modulation=False, modulation_halo_dir=None, day_bins=12, verbose=False,
                save_rates_per_bin=True):
    '''Compute and save the DM-nucleus scattering spectra in Xenon, Argon and Germanium
    (DDrates.npz, DDrates.txt) for the MadDM couplings c1, c4.

    Parameters
    ----------
    card : maddm_card (rhoDM; vmp, vescape, vSun_*, vEarth_mod for the modulation tables)
    alphaq : MadDM quark alphas
    dmmass : DM mass [GeV]
    output_path : directory where the files are written
    halo_path : annual-average halo table (Earth velocity vector of March 9, Eq. 12 of arXiv:2105.00599)
    modulation : if True, also compute the annual modulation (see rapidd.modulation): the year is
        divided in day_bins bins of width dt = 365.24/day_bins days and the halo and the spectrum are
        evaluated at the centers t_i = (i + 0.5) dt (days since March 22, 2018)
    modulation_halo_dir : directory of the modulation tables halo_table_<i>.dat. If None they are
        generated (or reused, if already there for the same parameters) in output_path/halo_table_modulation
    day_bins : number of time bins
    verbose : print the progress of the computation of the modulation tables
    save_rates_per_bin : save the spectrum dR/dE of each time bin (npz key drde_<target>_bin)

    With modulation it also writes rate_vs_time.dat, dRdE_modulation.dat and, in DDrates.npz,
    t_bin, dt_bin, rate_<target>_bin [events/kg/dt], drde_<target>_year [events/kg/yr/keV],
    drde_<target>_avg [events/kg/day/keV] (and drde_<target>_bin), target = Xe, Ar, Ge.
    '''
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
    names = [t for t, _ in NUCLEAR_TARGETS]

    def spectrum_fn():
        return _nuclear_rates(rhoDM, dmmass, energies)

    drde = spectrum_fn()
    npz = dict(energies=energies, **{"drde_" + t: drde[t] for t in names})

    if modulation:
        from .modulation import run_modulation, trapezoid_weights
        npz.update(run_modulation(spectrum_fn, card, output_path, energies, trapezoid_weights(energies),
                                  names, "", modulation_halo_dir=modulation_halo_dir, day_bins=day_bins,
                                  verbose=verbose, save_rates_per_bin=save_rates_per_bin,
                                  restore_halo=halo_path))

    np.savez(output_path+"/DDrates.npz", **npz)

    np.savetxt(output_path+"/DDrates.txt", np.column_stack([energies] + [drde[t] for t in names]),
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
                         modulation=False, modulation_halo_dir=None, day_bins=12, verbose=False,
                         save_rates_per_bin=True):
    '''Compute and save the DM-electron scattering rates in Xenon and the XENON10 / XENON1T
    S2 signals (port of the legacy MadDM electron_recoil_signal.f).

    Parameters
    ----------
    card : maddm_card (rhoDM; vmp, vescape, vSun_*, vEarth_mod for the modulation tables)
    alphas : dict with "alpha_e" = (aSIe, aSIo, aSDe, aSDo)
    dm_spin : MadDM dm_spin = 2s+1 (1: scalar, 2: fermion)
    m_e : electron mass [GeV]
    dmmass : DM mass [GeV]
    output_path : directory where the files are written
    halo_path : annual-average halo table (Earth velocity vector of March 9, Eq. 12 of arXiv:2105.00599)
    modulation, modulation_halo_dir, day_bins, verbose, save_rates_per_bin : annual modulation,
        as in DDrate_save. With modulation it also writes rate_vs_time_e_recoil.dat,
        dRdE_modulation_e_recoil.dat and, in DDrates_e.npz, t_bin, dt_bin, rate_tot_bin [events/kg/dt],
        drde_tot_year [events/kg/yr/keV], drde_tot_avg [events/kg/day/keV] (and drde_tot_bin),
        for the spectrum summed over the shells.

    Returns a dict with DM_response, sigma_e, Xenon10_bins, Xenon10_signal, Xenon10_obs,
    Xenon1T_signal, Xenon1T_bkg, Xenon1T_obs (DM_response = -1 if the DM spin is not implemented).
    '''
    import os
    from .core import (read_atomic_responses, dm_response_electron, difrate_dS2, xenon10_bins,
                       XE_SHELLS)
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
        from .modulation import run_modulation
        E_keV = E[:-1]*1e6

        def spectrum_fn():
            return {'tot': _electron_rates(rhoDM, dmmass, jchi, m_e)[2].sum(axis=0)/1e6}

        # same energy integral as the legacy rate: sum_k dR/dE(E_k) (E_{k+1} - E_k)
        npz.update(run_modulation(spectrum_fn, card, output_path, E_keV, np.diff(E)*1e6, ['tot'],
                                  "_e_recoil", modulation_halo_dir=modulation_halo_dir, day_bins=day_bins,
                                  verbose=verbose, save_rates_per_bin=save_rates_per_bin,
                                  restore_halo=halo_path))

    np.savez(os.path.join(output_path, 'DDrates_e.npz'), **npz)
    np.savetxt(os.path.join(output_path, 'DDrates_e.txt'),
               np.column_stack((k[:-1]*1e6, E[:-1]*1e6, rate_tot/1e6, (rate_shell/1e6).T)),
               header="k [keV] E [keV] drde_tot [events/kg/day/keV] " +
                      " ".join("drde_%s" % s for s in XE_SHELLS))

    return result
