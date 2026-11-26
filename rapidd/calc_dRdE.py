import os

import warnings

from .core import (halo_path, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo, difrate_dER, vev,
                   ELECTRON_TARGETS, CRYSTAL_TARGETS)
import numpy as np

# DM-nucleus targets and nuclear models (NREFT form factors) of the MadDM rates
NUCLEAR_TARGETS = (("Xe", "McElvain2022"), ("Ar", "Fitz2012"), ("Ge", "Fitz2012"), ("F", "McElvain2022"))
NUCLEAR_TARGET_NAMES = tuple(t for t, _ in NUCLEAR_TARGETS)
# nuclear structure of the DM-nucleus spectra: chiral EFT (rapidd.chiral_eft) or the one-body NREFT responses above
NUCLEAR_MODELS = ("chiral_eft", "nreft")
# format of the text tables (7 significant digits)
TXT_FMT = '%.6e'


def select_targets(targets, allowed, kind="target"):
    '''Targets of the rate files: a name or a list of names among allowed (case insensitive,
    e.g. "xe" -> "Xe"), without repetitions and in the order of allowed.'''
    if isinstance(targets, str):
        targets = [targets]
    canonical = {t.lower(): t for t in allowed}
    chosen = set()
    for t in targets:
        if str(t).lower() not in canonical:
            raise ValueError("Unknown %s %s: choose from %s" % (kind, t, ", ".join(allowed)))
        chosen.add(canonical[str(t).lower()])
    if not chosen:
        raise ValueError("No %s selected: choose from %s" % (kind, ", ".join(allowed)))
    return tuple(t for t in allowed if t in chosen)


def _nuclear_rates(rhoDM, dmmass, energies, targets=NUCLEAR_TARGET_NAMES):
    '''dR/dE [events/kg/day/keV] of the targets (Xe, Ar, Ge, F) with the halo table currently loaded'''
    basis = dict(NUCLEAR_TARGETS)
    return {target: difrate_dER(rhoDM, dmmass, energies, target=target, basis=basis[target])
            for target in targets}


def _chiral_rates(rhoDM, dmmass, energies, couplings, targets=NUCLEAR_TARGET_NAMES):
    '''SI and SD dR/dE [events/kg/day/keV] of the targets with the chiral EFT structure factors (HKMS for SI, HMS
    central for SD) and the halo table currently loaded'''
    from .chiral_eft import difrate_chiral
    return {target: difrate_chiral(rhoDM, dmmass, energies, couplings, target=target, parts=True)
            for target in targets}


def equivalent_sigma_mdm(card, alphaq, dmmass, halo_path=halo_path):
    '''DM-nucleon cross sections [cm^2] equivalent to the MadDM model in LZ (rapidd.chiral_eft.equivalent_sigma),
    to be compared with the published LZ (and XLZD) SI, SD-n and SD-p limits: sigmaN_SI_eq, sigmaN_SD_n_eq,
    sigmaN_SD_p_eq. The couplings are those of alpha_map.chiral_couplings_mdm (Dirac DM, j_chi = 1/2) and the
    expected LZ signal is computed with the chiral EFT structure factors and the halo table halo_path, whatever
    the nuclear model of the saved spectra. Loads the halo table halo_path.'''
    from .alpha_map import chiral_couplings_mdm
    from .chiral_eft import equivalent_sigma
    read_halo(halo_path)
    sig = equivalent_sigma(chiral_couplings_mdm(card, alphaq, dmmass), dmmass)
    return {'sigmaN_%s_eq' % k: float(v) for k, v in sig.items()}


def DDrate_save(card, alphaq, dmmass, output_path, targets=NUCLEAR_TARGET_NAMES, halo_path=halo_path,
                modulation=False, modulation_halo_dir=None, day_bins=12, verbose=False,
                save_rates_per_bin=True, save_npz=False, nuclear_model="chiral_eft"):
    '''Compute and save the DM-nucleus scattering spectra in Xenon, Argon, Germanium and/or Fluorine
    (dRdE_nucleon.dat: E_R [keV] and one column dR/dE [events/kg/day/keV] per target) for the
    MadDM couplings.

    nuclear_model = "chiral_eft" (default): chiral EFT structure factors (rapidd.chiral_eft) for the four channels
    of the MadDM alphas, SS and VV (SI, Hoferichter, Klos, Menendez, Schwenk 2019, with the two-body pion couplings
    of SS) and AA and TT (SD, Hoferichter, Menendez, Schwenk 2020 central, with the two-body currents of AA only),
    see alpha_map.chiral_couplings_mdm. Argon has no SD rate (J = 0).
    nuclear_model = "nreft": the one-body NREFT responses of NUCLEAR_TARGETS (McElvain2022 for Xe and F, Fitz2012
    for Ar and Ge) with the coefficients c1, c4 of O_1 and O_4 (c1_dirac_mdm, c4_dirac_mdm).

    Parameters
    ----------
    card : maddm_card (rhoDM; vmp, vescape, vSun_*, vEarth_mod for the modulation tables)
    alphaq : MadDM quark alphas
    dmmass : DM mass [GeV]
    output_path : directory where the files are written
    targets : targets among "Xe", "Ar", "Ge", "F" (case insensitive; columns in this order)
    halo_path : annual-average halo table (Earth velocity vector of March 9, Eq. 12 of arXiv:2105.00599)
    modulation : if True, also compute the annual modulation (see rapidd.modulation): the year is
        divided in day_bins bins of width dt = 365.24/day_bins days and the halo and the spectrum are
        evaluated at the centers t_i = (i + 0.5) dt (days since March 22, 2018)
    modulation_halo_dir : directory of the modulation tables halo_table_<i>.dat. If None they are
        generated (or reused, if already there for the same parameters) in output_path/halo_table_modulation
    day_bins : number of time bins
    verbose : print the progress of the computation of the modulation tables
    save_rates_per_bin : save the spectrum dR/dE of each time bin (npz key drde_<target>_bin, only
        with save_npz)
    save_npz : also save all the arrays in dRdE_nucleon.npz (energies, drde_<target>; with chiral_eft also
        the SI and SD parts drde_<target>_SI, drde_<target>_SD)
    nuclear_model : "chiral_eft" or "nreft" (see above)

    With modulation it also writes rate_vs_time_nucleon.dat, dRdE_modulation_nucleon.dat and, in
    dRdE_nucleon.npz, t_bin, dt_bin, rate_<target>_bin [events/kg/dt], drde_<target>_year
    [events/kg/yr/keV], drde_<target>_avg [events/kg/day/keV] (and drde_<target>_bin).
    '''
    from .alpha_map import c1_dirac_mdm, c4_dirac_mdm, chiral_couplings_mdm
    if nuclear_model not in NUCLEAR_MODELS:
        raise ValueError("nuclear_model must be one of %s" % (NUCLEAR_MODELS,))
    names = select_targets(targets, NUCLEAR_TARGET_NAMES, "DM-nucleus target")
    reset_coefficients()

    rhoDM = card['rhoDM']
    read_halo(halo_path)
    energies = np.linspace(0.01, 40, 1000)
    parts = {}

    if nuclear_model == "nreft":
        c1p, c1n = c1_dirac_mdm(card, alphaq)
        c4p, c4n = c4_dirac_mdm(card, alphaq)
        set_any_Ncoeff(c1p*vev**2, 1, "p")
        set_any_Ncoeff(c1n*vev**2, 1, "n")
        set_any_Ncoeff(c4p*vev**2, 4, "p")
        set_any_Ncoeff(c4n*vev**2, 4, "n")

        def spectrum_fn():
            return _nuclear_rates(rhoDM, dmmass, energies, names)
    else:
        couplings = chiral_couplings_mdm(card, alphaq, dmmass)

        def spectrum_fn():
            rates = _chiral_rates(rhoDM, dmmass, energies, couplings, names)
            parts.update({"drde_%s_%s" % (t, p): rates[t][p] for t in names for p in ("SI", "SD")})
            return {t: rates[t]["total"] for t in names}

    drde = spectrum_fn()
    npz = dict(energies=energies, **{"drde_" + t: drde[t] for t in names})
    npz.update(parts)

    if modulation:
        from .modulation import run_modulation, trapezoid_weights
        npz.update(run_modulation(spectrum_fn, card, output_path, energies, trapezoid_weights(energies),
                                  names, "_nucleon", modulation_halo_dir=modulation_halo_dir, day_bins=day_bins,
                                  verbose=verbose, save_rates_per_bin=save_rates_per_bin,
                                  restore_halo=halo_path))

    if save_npz:
        np.savez(os.path.join(output_path, "dRdE_nucleon.npz"), **npz)

    np.savetxt(os.path.join(output_path, "dRdE_nucleon.dat"),
               np.column_stack([energies] + [drde[t] for t in names]), fmt=TXT_FMT,
               header="E_R [keV] " + " ".join("drde_%s [events/kg/day/keV]" % t for t in names))

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


def _electron_rates(rhoDM, dmmass, jchi, m_e, target="Xe"):
    from .core import difrate_electron_dEr_all
    return difrate_electron_dEr_all(rhoDM, dmmass, jchi, m_e, target=target)


def electron_cross_section(alphas, dm_spin, m_e, dmmass):
    '''Set the MadDM DM-electron coefficients c1, c4 in RAPIDD and return a dict with the DM response
    DM_response and the reference DM-electron cross section sigma_e [cm^2], without computing any rate
    (DM_response = -1 and no sigma_e if the DM spin is not implemented).'''
    from .core import dm_response_electron
    from .alpha_map import set_alphas_electron_rpd, jchi_from_dm_spin

    reset_coefficients()
    if not set_alphas_electron_rpd(alphas, dmmass, m_e, dm_spin):
        return {'DM_response': -1}
    dm_resp = dm_response_electron(dmmass, jchi_from_dm_spin(dm_spin))
    return {'DM_response': dm_resp, 'sigma_e': dm_resp/(16.*np.pi*(dmmass + m_e)**2)}


# DM-electron targets of the rate files: atoms (dRdE_electron.dat) and crystals (dRdE_electron_<target>.dat)
ELECTRON_RATE_TARGETS = ELECTRON_TARGETS + CRYSTAL_TARGETS
# numbers of electron-hole pairs of dRdQ_electron_<crystal>.dat
CRYSTAL_Q_MAX = 10


class CrystalTargetSkippedWarning(UserWarning):
    '''A crystal target is skipped because the couplings set need crystal responses that are not available'''


def DDrate_electron_save(card, alphas, dm_spin, m_e, dmmass, output_path, targets=("Xe",), halo_path=halo_path,
                         modulation=False, modulation_halo_dir=None, day_bins=12, verbose=False,
                         save_rates_per_bin=True, save_npz=False):
    '''Compute and save the DM-electron scattering rates in Xenon, Argon and/or silicon and, for Xenon, the
    XENON10 / XENON1T S2 signals (port of the legacy MadDM electron_recoil_signal.f).

    Parameters
    ----------
    card : maddm_card (rhoDM; vmp, vescape, vSun_*, vEarth_mod for the modulation tables)
    alphas : dict with "alpha_e" = (aSIe, aSIo, aSDe, aSDo)
    dm_spin : MadDM dm_spin = 2s+1 (1: scalar, 2: fermion)
    m_e : electron mass [GeV]
    dmmass : DM mass [GeV]
    output_path : directory where the files are written
    targets : targets among "Xe", "Ar", "Si" (case insensitive; columns in this order). "Si" (crystal, QCDark2
        response, RPA screening) is skipped with a CrystalTargetSkippedWarning if the coefficients set need the
        crystal responses W_2 ... W_5 (rapidd.electron.check_crystal_couplings)
    halo_path : annual-average halo table (Earth velocity vector of March 9, Eq. 12 of arXiv:2105.00599)
    modulation, modulation_halo_dir, day_bins, verbose, save_rates_per_bin : annual modulation,
        as in DDrate_save. With modulation it also writes rate_vs_time_electron.dat and
        dRdE_modulation_electron.dat (and, in dRdE_electron.npz, t_bin, dt_bin, rate_<target>_bin
        [events/kg/dt], drde_<target>_year [events/kg/yr/keV], drde_<target>_avg [events/kg/day/keV]
        and drde_<target>_bin), for the spectra summed over the shells.
    save_npz : also save all the arrays in dRdE_electron.npz

    Files: dRdE_electron.dat if an atom is a target (k [keV], E [keV], dR/dE [events/kg/day/keV] summed over
    the shells for each atom, then per shell for each atom); for silicon, dRdE_electron_Si.dat (deposited energy
    Delta E [keV], dR/dDelta E [events/kg/day/keV], on its own grid 0.1 ... 50 eV) and dRdQ_electron_Si.dat
    (number Q = 1 ... CRYSTAL_Q_MAX of electron-hole pairs, dR/dQ [events/kg/day], Ramanathan-Kurinsky
    100 K); if Xenon is a target, dRdS2_electron_XENON10.dat
    (S2, dN/dS2 [events] and its shell contributions), dRdS2_electron_XENON1T.dat (S2, dN/dS2 [events])
    and signal_e_recoil.dat (signal per XENON10 S2 bin).

    Returns a dict with DM_response, sigma_e, skipped_targets (target -> reason), if Xenon is a target,
    Xenon10_bins, Xenon10_signal, Xenon10_obs, Xenon1T_signal, Xenon1T_bkg, Xenon1T_obs and, for silicon,
    dRdQ_Si (rates per e-h pair, Q = 1 ... CRYSTAL_Q_MAX) (DM_response = -1 if the DM spin is not implemented).
    '''
    from .core import read_atomic_responses, difrate_dS2, xenon10_bins, ELECTRON_TARGET_SHELLS, XE_SHELLS
    from .alpha_map import jchi_from_dm_spin
    from .electron import difrate_electron, difrate_dQ, check_crystal_couplings, UnsupportedCrystalCouplingError

    selected = select_targets(targets, ELECTRON_RATE_TARGETS, "DM-electron target")
    # sets the DM-electron coefficients
    result = electron_cross_section(alphas, dm_spin, m_e, dmmass)
    if result['DM_response'] == -1:
        return result
    jchi = jchi_from_dm_spin(dm_spin)
    result['skipped_targets'] = {}
    crystals = []
    for t in selected:
        if t not in CRYSTAL_TARGETS:
            continue
        try:
            check_crystal_couplings(t)
            crystals.append(t)
        except UnsupportedCrystalCouplingError as err:
            result['skipped_targets'][t] = str(err)
            warnings.warn("%s skipped: %s" % (t, err), CrystalTargetSkippedWarning, stacklevel=2)
    names = tuple(t for t in selected if t in ELECTRON_TARGETS)

    rhoDM = card['rhoDM']
    read_halo(halo_path)
    for t in names:
        read_atomic_responses(target=t)

    npz = {}
    for t in crystals:
        r = difrate_electron(rhoDM, dmmass, jchi, m_e=m_e, target=t)
        dq = difrate_dQ(r, Q_max=CRYSTAL_Q_MAX)
        result['dRdQ_' + t] = [float(x) for x in dq['dRdQ']]
        np.savetxt(os.path.join(output_path, 'dRdE_electron_%s.dat' % t), np.column_stack([r['E'], r['dRdE']]),
                   fmt=TXT_FMT, header="DeltaE [keV] drde_%s [events/kg/day/keV] (deposited energy, %s crystal, "
                                       "RPA screening)" % (t, t))
        np.savetxt(os.path.join(output_path, 'dRdQ_electron_%s.dat' % t), np.column_stack([dq['Q'], dq['dRdQ']]),
                   fmt=['%d', TXT_FMT], header="Q [e-h pairs] drdq_%s [events/kg/day] (Ramanathan-Kurinsky 100 K)" % t)
        npz.update({'DeltaE_keV_' + t: r['E'], 'drde_' + t: r['dRdE'], 'Q_' + t: dq['Q'], 'drdq_' + t: dq['dRdQ']})
        if modulation:
            from .modulation import run_modulation

            def crystal_fn(t=t):
                return {t: difrate_electron(rhoDM, dmmass, jchi, m_e=m_e, target=t)['dRdE']}

            npz.update(run_modulation(crystal_fn, card, output_path, r['E'], np.full(r['E'].size, r['E'][1] - r['E'][0]),
                                      (t,), "_electron_" + t, modulation_halo_dir=modulation_halo_dir,
                                      day_bins=day_bins, verbose=verbose, save_rates_per_bin=save_rates_per_bin,
                                      restore_halo=halo_path))
    if not names:
        if save_npz:
            np.savez(os.path.join(output_path, 'dRdE_electron.npz'), **npz)
        return result

    rate_shell = {}
    for t in names:
        k, E, rate_shell[t] = _electron_rates(rhoDM, dmmass, jchi, m_e, t)
    k_keV, E_keV = k[:-1]*1e6, E[:-1]*1e6
    drde = {t: rate_shell[t].sum(axis=0)/1e6 for t in names}

    npz.update(k_keV=k_keV, E_keV=E_keV)
    for t in names:
        npz['drde_' + t] = drde[t]
        npz['drde_%s_shell' % t] = rate_shell[t]/1e6
        npz['shells_' + t] = np.array(ELECTRON_TARGET_SHELLS[t])

    np.savetxt(os.path.join(output_path, 'dRdE_electron.dat'),
               np.column_stack([k_keV, E_keV] + [drde[t] for t in names] +
                               [(rate_shell[t]/1e6).T for t in names]), fmt=TXT_FMT,
               header="k [keV] E [keV] " + " ".join("drde_%s [events/kg/day/keV]" % t for t in names) + " " +
                      " ".join("drde_%s_%s" % (t, s) for t in names for s in ELECTRON_TARGET_SHELLS[t]))

    if "Xe" in names:
        x10 = difrate_dS2("Xenon10", E, rate_shell["Xe"], dmmass)
        x1t = difrate_dS2("Xenon1T", E, rate_shell["Xe"], dmmass)
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

        np.savetxt(os.path.join(output_path, 'dRdS2_electron_XENON10.dat'),
                   np.column_stack([x10['S2'], x10['dRdS2'], x10['dRdS2_shell'].T]),
                   fmt=['%d'] + [TXT_FMT]*(1 + len(XE_SHELLS)),
                   header="S2 dN/dS2 [events] " + " ".join("dN/dS2_%s" % s for s in XE_SHELLS))
        np.savetxt(os.path.join(output_path, 'dRdS2_electron_XENON1T.dat'),
                   np.column_stack([x1t['S2'], x1t['dRdS2']]), fmt=TXT_FMT, header="S2 dN/dS2 [events]")

        npz.update(S2_XENON10=x10['S2'], dRdS2_XENON10=x10['dRdS2'], dRdS2_shell_XENON10=x10['dRdS2_shell'],
                   S2_XENON1T=x1t['S2'], dRdS2_XENON1T=x1t['dRdS2'])

    if modulation:
        from .modulation import run_modulation

        def spectrum_fn():
            return {t: _electron_rates(rhoDM, dmmass, jchi, m_e, t)[2].sum(axis=0)/1e6 for t in names}

        # same energy integral as the legacy rate: sum_k dR/dE(E_k) (E_{k+1} - E_k)
        npz.update(run_modulation(spectrum_fn, card, output_path, E_keV, np.diff(E)*1e6, names,
                                  "_electron", modulation_halo_dir=modulation_halo_dir, day_bins=day_bins,
                                  verbose=verbose, save_rates_per_bin=save_rates_per_bin,
                                  restore_halo=halo_path))

    if save_npz:
        np.savez(os.path.join(output_path, 'dRdE_electron.npz'), **npz)

    return result
