'''Annual modulation of the direct-detection rates, shared by DM-nucleus (calc_dRdE.DDrate_save)
and DM-electron (calc_dRdE.DDrate_electron_save) scattering.

The year (YEAR_DAYS = 365.24 days) is divided in day_bins bins of width dt = YEAR_DAYS/day_bins.
For each bin the halo table and the spectrum dR/dE are evaluated at the center of the bin,
t_i = (i + 0.5) dt days since March 22, 2018 (the reference date of Eq. 11 of arXiv:2105.00599).
The Earth velocity has modulus vEarth_mod and the time-dependent direction of Eq. 11; the
annual-average table (core.halo_path, core.halo_path_from_card) instead uses the fixed vector of
March 9 (Eq. 12) and does not depend on vEarth_mod.

The only part specific to the target is spectrum_fn (see modulated_rates): a function returning
the spectra dR/dE [events/kg/day/keV] computed with the halo table currently loaded.
'''
import json
import os

import numpy as np

from .core import YEAR_DAYS, calc_halo, read_halo, halo_path, default_halo_v_earth

MODULATION_HALO_DIR = 'halo_table_modulation'
_PARAMS_FILE = 'modulation_halo_params.json'


def modulation_times(day_bins):
    '''Centers t_i = (i + 0.5) dt [days since March 22, 2018] of the day_bins time bins and
    their width dt = YEAR_DAYS/day_bins [days]'''
    day_bins = int(day_bins)
    if day_bins < 1:
        raise ValueError("day_bins must be >= 1")
    dt = YEAR_DAYS / day_bins
    return (np.arange(day_bins) + 0.5) * dt, dt


def _card_value(card, key, default):
    try:
        return card[key]
    except KeyError:
        return default


def modulation_halo_tables(table_dir, day_bins=12, profile="SHM", vesc=544., v0=238., v0_lsr=238.,
                           v_pec=(11.1, 12.2, 7.3), v_earth=default_halo_v_earth, verbose=False, **halo_kwargs):
    '''
    Halo tables for the annual modulation: one table table_dir/halo_table_<i>.dat per time bin,
    at t_i of modulation_times(day_bins), computed by core.calc_halo(times=t_i). The halo
    parameters are those of calc_halo (halo_kwargs: beta, k_lisanti, i, w, vb, cosb, sigb, vcut);
    v_earth is the modulus of the Earth velocity, whose direction changes during the year
    (Eq. 11 of arXiv:2105.00599).
    The tables are not recomputed if table_dir already contains the tables for the same
    parameters (e.g. when the DM-nucleus and DM-electron rates of the same run share them).

    Returns the list of the paths of the tables (index i = time bin).
    '''
    times, _ = modulation_times(day_bins)
    pars = dict(profile=profile, vesc=float(vesc), v0=float(v0), v0_lsr=float(v0_lsr),
                v_pec=[float(v) for v in v_pec], v_earth=float(v_earth),
                times=[float(t) for t in times], **{k: halo_kwargs[k] for k in sorted(halo_kwargs)})
    paths = [os.path.join(table_dir, 'halo_table_%i.dat' % i) for i in range(len(times))]
    params_file = os.path.join(table_dir, _PARAMS_FILE)
    if os.path.isfile(params_file) and all(os.path.isfile(p) for p in paths):
        with open(params_file) as f:
            try:
                if json.load(f) == pars:
                    return paths
            except ValueError:
                pass
    calc_halo(None, profile=profile, vesc=vesc, v0=v0, v0_lsr=v0_lsr, v_pec=tuple(v_pec), times=times,
              v_earth=v_earth, table_dir=table_dir, verbose=verbose, **halo_kwargs)
    with open(params_file, 'w') as f:
        json.dump(pars, f)
    return paths


def halo_pars_from_card(card):
    '''SHM parameters of the MadDM card for modulation_halo_tables: vescape, vmp (as v0 and as the
    LSR speed), (vSun_r, vSun_phi, vSun_theta) and vEarth_mod (core.default_halo_v_earth if missing)'''
    return dict(profile="SHM", vesc=card['vescape'], v0=card['vmp'], v0_lsr=card['vmp'],
                v_pec=(card['vSun_r'], card['vSun_phi'], card['vSun_theta']),
                v_earth=_card_value(card, 'vEarth_mod', default_halo_v_earth))


def modulation_halo_from_card(card, table_dir, day_bins=12, verbose=False):
    '''modulation_halo_tables for the halo parameters of the MadDM card (halo_pars_from_card)'''
    return modulation_halo_tables(table_dir, day_bins=day_bins, verbose=verbose, **halo_pars_from_card(card))


def modulated_rates(spectrum_fn, halo_paths, dt, energy_weights, restore_halo=halo_path):
    '''
    Loop over the time bins: load each halo table and compute the spectra.

    Parameters
    ----------
    spectrum_fn : callable
        spectrum_fn() returns a dict {name: dR/dE [events/kg/day/keV]} (arrays on the energy grid),
        computed with the halo table currently loaded (DM-nucleus: Xe, Ar, Ge; DM-electron: the
        sum over the Xenon shells).
    halo_paths : list of str
        Halo table of each time bin (modulation_halo_from_card).
    dt : float
        Width of the time bins [days].
    energy_weights : array
        Weights w_k [keV] of the energy integral, int dR/dE dE = sum_k w_k dR/dE(E_k).
    restore_halo : str or None
        Halo table read at the end (None: leave the last one loaded).

    Returns a dict {name: dict(drde_bin, rate_bin, drde_year, drde_avg)} with
      drde_bin  : dR/dE(t_i)                           [events/kg/day/keV]  (day_bins x n_E)
      rate_bin  : R_i = dt int dR/dE(t_i) dE           [events/kg/dt]       (day_bins)
      drde_year : sum_i dR/dE(t_i) dt                  [events/kg/yr/keV]
      drde_avg  : drde_year / YEAR_DAYS (time average) [events/kg/day/keV], to be compared with
                  the spectrum computed with the annual-average halo
    '''
    spectra = {}
    try:
        for path in halo_paths:
            read_halo(path)
            for name, s in spectrum_fn().items():
                spectra.setdefault(name, []).append(np.asarray(s, dtype=float))
    finally:
        if restore_halo is not None:
            read_halo(restore_halo)
    w = np.asarray(energy_weights, dtype=float)
    out = {}
    for name, s in spectra.items():
        drde_bin = np.array(s)
        drde_year = drde_bin.sum(axis=0) * dt
        out[name] = dict(drde_bin=drde_bin, rate_bin=dt * (drde_bin @ w),
                         drde_year=drde_year, drde_avg=drde_year / YEAR_DAYS)
    return out


def trapezoid_weights(E):
    '''Weights of the trapezoidal rule on the grid E'''
    E = np.asarray(E, dtype=float)
    w = np.zeros_like(E)
    d = np.diff(E)
    w[:-1] += d / 2
    w[1:] += d / 2
    return w


def modulation_npz(res, times, dt, save_rates_per_bin=True):
    '''npz entries of the modulation: t_bin, dt_bin and, for each spectrum name,
    rate_<name>_bin, drde_<name>_year, drde_<name>_avg (and drde_<name>_bin if save_rates_per_bin)'''
    npz = dict(t_bin=np.asarray(times), dt_bin=np.float64(dt))
    for name, r in res.items():
        npz['rate_%s_bin' % name] = r['rate_bin']
        npz['drde_%s_year' % name] = r['drde_year']
        npz['drde_%s_avg' % name] = r['drde_avg']
        if save_rates_per_bin:
            npz['drde_%s_bin' % name] = r['drde_bin']
    return npz


def write_modulation_files(res, names, times, dt, E, rate_path, spectrum_path):
    '''
    Text files of the modulation (same format for DM-nucleus and DM-electron):
      rate_path     : t_i [days since 22-03-2018], dt [days], R_<name> [events/kg/dt] for each name
      spectrum_path : E [keV], drde_<name>_avg [events/kg/day/keV], drde_<name>_year [events/kg/yr/keV]
    '''
    np.savetxt(rate_path, np.column_stack([times, np.full(len(times), dt)] + [res[n]['rate_bin'] for n in names]),
               fmt='%.10E', header="t [days since 22-03-2018] dt [days] " +
                      " ".join("R_%s [events/kg/dt]" % n for n in names))
    np.savetxt(spectrum_path, np.column_stack([E] + [res[n]['drde_avg'] for n in names] +
                                              [res[n]['drde_year'] for n in names]),
               fmt='%.10E', header="E [keV] " + " ".join("drde_%s_avg [events/kg/day/keV]" % n for n in names) + " " +
                      " ".join("drde_%s_year [events/kg/yr/keV]" % n for n in names))


def modulation_spectra(spectrum_fn, energy_weights, table_dir, day_bins=12, verbose=False,
                       restore_halo=halo_path, **halo_pars):
    '''
    Generic annual modulation of any spectrum: halo tables (modulation_halo_tables with halo_pars,
    in table_dir) and loop over the time bins (modulated_rates).
    Returns times t_i [days since March 22, 2018], dt [days] and the dict of modulated_rates.
    '''
    times, dt = modulation_times(day_bins)
    paths = modulation_halo_tables(table_dir, day_bins=day_bins, verbose=verbose, **halo_pars)
    return times, dt, modulated_rates(spectrum_fn, paths, dt, energy_weights, restore_halo=restore_halo)


def run_modulation(spectrum_fn, card, output_path, E, energy_weights, names, file_suffix,
                   modulation_halo_dir=None, day_bins=12, verbose=False, save_rates_per_bin=True,
                   restore_halo=halo_path):
    '''
    Full modulation workflow shared by DDrate_save and DDrate_electron_save: modulation_spectra
    with the halo of the MadDM card (tables in output_path/halo_table_modulation if
    modulation_halo_dir is None), text files rate_vs_time<file_suffix>.dat and
    dRdE_modulation<file_suffix>.dat (write_modulation_files).
    Returns the npz entries (modulation_npz).
    '''
    if modulation_halo_dir is None:
        modulation_halo_dir = os.path.join(output_path, MODULATION_HALO_DIR)
    times, dt, res = modulation_spectra(spectrum_fn, energy_weights, modulation_halo_dir, day_bins=day_bins,
                                        verbose=verbose, restore_halo=restore_halo, **halo_pars_from_card(card))
    write_modulation_files(res, names, times, dt, E,
                           os.path.join(output_path, 'rate_vs_time%s.dat' % file_suffix),
                           os.path.join(output_path, 'dRdE_modulation%s.dat' % file_suffix))
    return modulation_npz(res, times, dt, save_rates_per_bin=save_rates_per_bin)
