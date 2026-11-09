'''Regression of the DM-electron SI/SD (Ce(1), Ce(4)) results against reference outputs
(tests/data/regression_direct_electron_port, see generate_reference.py): MadDM only uses this path.

The references were first produced with the branch direct_electron_port and regenerated with the
current code after the following documented fixes:
- The XENON1T S2 response is evaluated at E_e + E_sec instead of the kinetic energy E_e of the
  outgoing electron, with E_sec the energy of the de-excitation photon of the shell ("electron"
  convention, default: 63.1, 87.9, 50.3 eV for 4d, 4p, 4s, from which the secondary quanta of
  Table II of arXiv:1703.00910 come; 0 for the other shells), as XENON1T and the XENONnT
  DM-electron limits; the "deposited" convention (E_sec = E_b) is checked too. The XENON1T
  outputs are also checked against an independent fold of the spectra with the S2-only response
  (xenon1t_fold).
- The q integral is Simpson's rule in ln q with W interpolated between the tabulated points
  (the left sum on the tabulated q of direct_electron_port overestimated all the rates by ~1.8%),
  and the Xe_4d atomic responses were regenerated with the current DarkART (the old tables were
  produced before DarkART fixed the binding energy of Xe 4d, Dec. 2022, which entered the
  effective charge of the continuum states; see lib/atomic_responses/README.md). The spectra are
  also checked against an independent evaluation of the q integral (si_q_integral).
- The XENON10 response integrates dR/dE (interpolated as a power law between the grid points) over
  E in each bin of the number of primary quanta n_1 = floor(E/W), instead of the left sum
  dR/dE(E_i) (E_{i+1} - E_i) at n_1(E_i) of direct_electron_port (~7.5% more rate). The outputs are
  also checked against an independent fold (xenon10_fold).'''
import os
import re
import shutil

import numpy as np
import pytest

import rapidd
from rapidd.calc_dRdE import DDrate_electron_save

from scipy.interpolate import CubicSpline
from scipy.stats import binom

from conftest import DATA_DIR
from test_electron_rate import simpson_lnq

REF_DIR = os.path.join(DATA_DIR, "regression_direct_electron_port")
RTOL = 1e-12
M_E = 0.000510998950
RHO = 0.3
COEFFS = (("c1", 1.0e-3, 0.0), ("c4", 0.0, 1.0e-3), ("c1c4", 1.0e-3, 2.0e-3))
CASES = [(label, c1, c4, m, j) for label, c1, c4 in COEFFS for m in (0.05, 0.1, 1.0) for j in (0.5, 0.0)]
CARD = {'rhoDM': RHO, 'vmp': 238., 'vescape': 544., 'vSun_r': 11.1, 'vSun_phi': 12.2, 'vSun_theta': 7.3}
ALPHAS = {"alpha_e": (1.0e-4, 2.0e-4, 1.5e-4, 0.5e-4)}
DD_CASES = ((2, 0.1, False), (2, 1.0, False), (1, 0.05, False), (2, 0.1, True))
DAY_BINS = 4  # time bins of the modulation case (as in generate_reference.py)
X1T_DIR = os.path.join(rapidd.core.experiments_dir, "Xenon1T_S2only")
XE_N_SECONDARY = {"4d": 4, "4p": 6, "4s": 3}  # secondary quanta of XENON10 (Table II of arXiv:1703.00910)
# binding energies [eV] (as in lib/source/atomic_resp.c) and de-excitation photons: the hole of 4d, 4p,
# 4s is filled by an electron of 5p, 4d, 4p (Table II of arXiv:1703.00910)
XE_EB_EV = {"4d": 75.58972072252894, "4p": 163.49493390058458, "4s": 213.7805688618793, "5p": 12.44330433672413}
XE_PHOTON_GEV = {"4d": (XE_EB_EV["4d"] - XE_EB_EV["5p"]) * 1e-9, "4p": (XE_EB_EV["4p"] - XE_EB_EV["4d"]) * 1e-9,
                 "4s": (XE_EB_EV["4s"] - XE_EB_EV["4p"]) * 1e-9}


def test_deexcitation_energies():
    '''de-excitation photons of the C library, and the secondary quanta of XENON10 that they give'''
    for s, shell in enumerate(rapidd.XE_SHELLS):
        np.testing.assert_allclose(rapidd.XE_DEEXCITATION_ENERGIES[s], XE_PHOTON_GEV.get(shell, 0.), rtol=1e-14, atol=0)
        assert int(rapidd.XE_DEEXCITATION_ENERGIES[s] / 13.8e-9) == XE_N_SECONDARY.get(shell, 0)
    assert rapidd.er_energy_convention() == "electron"


def xenon1t_fold(E, rate_shell, convention="electron"):
    '''dN/dS2 of XENON1T (63 bins) from the per-shell dR/dE [events/kg/day/GeV] on the grid E [GeV],
    with the response of the S2-only data release evaluated at E_e + E_sec (E_sec: de-excitation
    photon, convention "electron", or binding energy, "deposited") and its raw exposure, written
    independently of difrate_dS2_Xenon1T.'''
    d = np.loadtxt(os.path.join(X1T_DIR, "s2_response_er.csv"), delimiter=',', skiprows=1)
    E_resp, dE = d[:, 0] * 1e-6, (d[:, 2] - d[:, 1]) * 1e-6
    resp = d[:, 3:3 + rapidd.S2_BIN_MAX_XENON1T]  # the C library uses the first 63 S2 bins
    nk = rate_shell.shape[1]
    tot = np.zeros_like(E_resp)
    for s in range(len(rapidd.XE_SHELLS)):
        shell = rapidd.XE_SHELLS[s]
        Esec = rapidd.XE_BINDING_ENERGIES[s] if convention == "deposited" else XE_PHOTON_GEV.get(shell, 0.)
        x, y = E[:nk] + Esec, rate_shell[s]
        # last energy moved to the last one of the response if closer than 5% (as in the C code)
        if abs((E_resp[-1] - x[-1]) / x[-1]) <= 0.05:
            x[-1] = E_resp[-1]
        tot += np.interp(E_resp, x, y, left=0., right=0.)
    return (tot * dE * 356770.) @ resp


def xenon10_fold(E, rate_shell):
    '''dN/dS2 of XENON10 (S2 = 1 ... 203) from the per-shell dR/dE [events/kg/day/GeV] on the grid E
    [GeV], with the detector model of arXiv:1703.00910, written independently of difrate_dS2_Xenon10:
    dR/dE interpolated as a power law between the grid points (linearly where it vanishes),
    integrated exactly in each bin of n_1 = floor(E/W) from E[0] to E[nk-1], then binomial number
    of electrons (n_1 + n_2 quanta, f_e = 1/1.2), gaussian S2 (27 PE/electron, sigma 6.7 sqrt(n_e)),
    trigger and flat efficiency (not for S2 = 1, as in the legacy code) and 15 kg day.'''
    W, nk = 13.8e-9, rate_shell.shape[1]
    trig = np.loadtxt(os.path.join(rapidd.core.experiments_dir, "Xenon10", "Xenon10_TriggerEfficiency.txt"))
    S2, ne = np.arange(1, 204), np.arange(1, 91)
    gauss = np.exp(-0.5 * ((S2[:, None] - 27 * ne) / (6.7 * np.sqrt(ne)))**2) / (6.7 * np.sqrt(2 * np.pi * ne))
    eff = np.concatenate(([1.], 15. * 0.92 * trig[:202]))
    out = np.zeros(len(S2))
    for s, shell in enumerate(rapidd.XE_SHELLS):
        r = rate_shell[s]
        if not np.any(r):
            continue
        R = np.zeros(int(E[nk - 1] / W) + 1)
        for i in range(nk - 1):
            a, b, ra, rb = E[i], E[i + 1], r[i], r[i + 1]
            edges = np.concatenate(([a], W * np.arange(np.floor(a / W) + 1, np.ceil(b / W)), [b]))
            x0, x1 = edges[:-1], edges[1:]
            if ra > 0 and rb > 0:
                p1 = np.log(rb / ra) / np.log(b / a) + 1
                v = ra * a / p1 * ((x1 / a)**p1 - (x0 / a)**p1)
            else:
                v = (ra + (rb - ra) / (b - a) * (0.5 * (x0 + x1) - a)) * (x1 - x0)
            np.add.at(R, np.floor(x0 / W).astype(int), v)
        n_tot = np.arange(len(R)) + XE_N_SECONDARY.get(shell, 0)
        dRdne = (R[:, None] * binom.pmf(ne[None, :] - 1, n_tot[:, None], 1 / 1.2)).sum(axis=0)
        out += gauss @ dRdne
    return out * eff


def xenon1t_roi_signal(dRdS2, m):
    t = np.loadtxt(os.path.join(X1T_DIR, "5d_results_dmelectron.dat"))
    lo, hi = (int(round(np.interp(m, t[:, 0], t[:, j]))) for j in (1, 2))
    return dRdS2[lo - 1:hi].sum()


def si_q_integral(E, m):
    '''int dq q eta_0(vmin) W_1(k', q) [GeV^2 (km/s)^-1] for each shell (rows) and k' of the grid,
    with E [GeV] the outgoing energies and m the DM mass [GeV], written independently of the C code
    (for contact c_1, c_4 R_1 does not depend on q), with the q quadrature of simpson_lnq. eta_0 is
    the natural cubic spline of the default halo table (as GSL), zero outside the table.'''
    h = np.loadtxt(rapidd.halo_path, skiprows=1)
    eta = CubicSpline(h[:, 0], h[:, 1], bc_type='natural')
    q = np.exp(np.log(1e-6) + np.log(1e3) * np.arange(200) / 199)
    nk = rapidd.GRIDSIZE_K
    out = np.zeros((len(rapidd.XE_SHELLS), nk))
    for s, shell in enumerate(rapidd.XE_SHELLS):
        path = os.path.join(rapidd.core.atomic_resp_dir, "atomic_response_1", "Xe_%s.txt" % shell)
        W = np.loadtxt(path, comments="//")
        for ik in range(nk):
            qs, wq, (Wq,) = simpson_lnq(q, [W[ik]])
            vmin = ((E[ik] + rapidd.XE_BINDING_ENERGIES[s]) / qs + qs / (2 * m)) * 299792.458
            e = np.where((vmin >= h[0, 0]) & (vmin <= h[-1, 0]), eta(vmin), 0.)
            out[s, ik] = np.sum(wq * qs**2 * e * Wq)
    return out


@pytest.fixture(scope="module")
def reference():
    return np.load(os.path.join(REF_DIR, "reference.npz"))


@pytest.fixture
def w2_convention(request):
    '''Sign of W_2 for the test, restored to the default afterwards'''
    rapidd.set_w2_sign(request.param)
    yield request.param
    rapidd.set_w2_sign("Liang")


# the SI/SD path of MadDM uses only W_1: it must not depend on the sign convention of W_2
@pytest.mark.parametrize("w2_convention", ("Liang", "Catena"), indirect=True)
@pytest.mark.parametrize("label,c1,c4,m,jchi", CASES)
def test_rates_and_signals(reference, default_halo, w2_convention, label, c1, c4, m, jchi):
    rapidd.set_any_Ecoeff(c1, 1)
    rapidd.set_any_Ecoeff(c4, 4)
    key = "%s_m%g_j%g/" % (label, m, jchi)
    k, E, rate_shell = rapidd.difrate_electron_dEr_all(RHO, m, jchi, M_E)
    x10 = rapidd.difrate_dS2("Xenon10", E, rate_shell, m)
    x1t = rapidd.difrate_dS2("Xenon1T", E, rate_shell, m)
    ref_shell = reference[key + "rate_shell"]
    np.testing.assert_allclose(k, reference[key + "k"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(E, reference[key + "E"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(rapidd.dm_response_electron(m, jchi), reference[key + "dm_response"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(rate_shell, ref_shell, rtol=RTOL, atol=0)
    # spectra: dR/dE E = const x si_q_integral (within 1e-5: a few points close to the kinematic
    # end points differ by up to 6e-6, the others by < 1e-9)
    nk = rate_shell.shape[1]
    if not np.any(ref_shell):  # R_1 = 0 (c_4 only, jchi = 0)
        assert not np.any(rate_shell)
    else:
        new = si_q_integral(E, m)
        assert np.array_equal(rate_shell == 0, new == 0)
        ratio = (rate_shell * E[:nk])[new != 0] / new[new != 0]
        np.testing.assert_allclose(ratio, np.median(ratio), rtol=1e-5)
    # XENON10: reference, and independent fold of the spectra
    np.testing.assert_allclose(x10['dRdS2'], xenon10_fold(E, rate_shell), rtol=1e-9, atol=0)
    np.testing.assert_allclose(x10['dRdS2'], reference[key + "X10_dRdS2"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x10['dRdS2_shell'], reference[key + "X10_dRdS2_shell"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x10['signal'], reference[key + "X10_signal"], rtol=RTOL, atol=0)
    assert np.all(x10['signal'] > 0) or not np.any(ref_shell)
    # XENON1T: reference, and independent fold of the spectra with E_e + de-excitation photon
    np.testing.assert_allclose(x1t['dRdS2'], reference[key + "X1T_dRdS2"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x1t['signal'], reference[key + "X1T_signal"], rtol=RTOL, atol=0)
    fold = xenon1t_fold(E, rate_shell)
    np.testing.assert_allclose(x1t['dRdS2'], fold, rtol=1e-10, atol=0)
    np.testing.assert_allclose(x1t['signal'], xenon1t_roi_signal(fold, m), rtol=1e-10, atol=0)
    np.testing.assert_allclose(x1t['n_bkg'], reference[key + "X1T_bkg"], rtol=RTOL, atol=0)
    assert x1t['n_obs'] == int(reference[key + "X1T_obs"])
    # "deposited" convention (E_e + E_b): independent fold, and more signal than "electron"
    rapidd.set_er_energy_convention("deposited")
    try:
        x1t_dep = rapidd.difrate_dS2("Xenon1T", E, rate_shell, m)
    finally:
        rapidd.set_er_energy_convention("electron")
    fold_dep = xenon1t_fold(E, rate_shell, "deposited")
    np.testing.assert_allclose(x1t_dep['dRdS2'], fold_dep, rtol=1e-10, atol=0)
    assert x1t_dep['signal'] >= x1t['signal']


def _numbers(path, skip=()):
    '''All the numbers of a text file, in order. The numpy scalar wrappers written by repr with
    numpy >= 2 (e.g. np.float64(0.1) in result.txt) are removed first, so that the files written
    with numpy 1 and 2 give the same numbers.'''
    out = []
    with open(path) as f:
        for line in f:
            if skip and line.startswith(tuple(skip)):
                continue
            line = re.sub(r"np\.(float|int)\d*\(([^()]*)\)", r"\2", line)
            for tok in line.replace(',', ' ').replace('[', ' ').replace(']', ' ').split():
                try:
                    out.append(float(tok))
                except ValueError:
                    pass
    return np.array(out)


@pytest.mark.parametrize("dm_spin,m,modulation", DD_CASES)
def test_DDrate_electron_save(tmp_path, dm_spin, m, modulation):
    '''Same files as the ddrate_* references,
    with the alphas of MadDM (smoke test of the MadDM entry point)'''
    name = "ddrate_spin%i_m%g%s" % (dm_spin, m, "_mod" if modulation else "")
    ref = os.path.join(REF_DIR, name)
    out = tmp_path / name
    out.mkdir()
    res = DDrate_electron_save(CARD, ALPHAS, dm_spin, M_E, m, str(out), halo_path=rapidd.halo_path,
                               modulation=modulation, day_bins=DAY_BINS)
    shutil.rmtree(out / "halo_table_modulation", ignore_errors=True)
    with open(out / "result.txt", "w") as f:
        for key in sorted(res):
            f.write("%s %r\n" % (key, res[key]))
    files = sorted(os.listdir(ref))
    assert sorted(os.listdir(out)) == files
    for fname in files:
        if fname == "result.txt":
            np.testing.assert_allclose(_numbers(out / fname), _numbers(os.path.join(ref, fname)),
                                       rtol=RTOL, atol=0, err_msg=fname)
            continue
        if fname.endswith(".npz"):
            a, b = np.load(out / fname), np.load(os.path.join(ref, fname))
            assert sorted(a.files) == sorted(b.files)
            for key in a.files:
                if a[key].dtype.kind in "US":
                    assert (a[key] == b[key]).all()
                else:
                    np.testing.assert_allclose(a[key], b[key], rtol=RTOL, atol=0, err_msg=key)
        else:
            np.testing.assert_allclose(_numbers(out / fname), _numbers(os.path.join(ref, fname)),
                                       rtol=RTOL, atol=0, err_msg=fname)
    rapidd.read_halo(rapidd.halo_path)
