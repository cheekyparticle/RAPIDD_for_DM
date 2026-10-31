'''Regression of the DM-electron SI/SD (Ce(1), Ce(4)) results against the outputs of the branch
direct_electron_port (tests/data/regression_direct_electron_port, see generate_reference.py):
MadDM only uses this path, which must not change.'''
import os
import re
import shutil

import numpy as np
import pytest

import rapidd
from rapidd.calc_dRdE import DDrate_electron_save

from conftest import DATA_DIR

REF_DIR = os.path.join(DATA_DIR, "regression_direct_electron_port")
RTOL = 1e-12
M_E = 0.000510998950
RHO = 0.3
COEFFS = (("c1", 1.0e-3, 0.0), ("c4", 0.0, 1.0e-3), ("c1c4", 1.0e-3, 2.0e-3))
CASES = [(label, c1, c4, m, j) for label, c1, c4 in COEFFS for m in (0.05, 0.1, 1.0) for j in (0.5, 0.0)]
CARD = {'rhoDM': RHO, 'vmp': 238., 'vescape': 544., 'vSun_r': 11.1, 'vSun_phi': 12.2, 'vSun_theta': 7.3}
ALPHAS = {"alpha_e": (1.0e-4, 2.0e-4, 1.5e-4, 0.5e-4)}
DD_CASES = ((2, 0.1, False), (2, 1.0, False), (1, 0.05, False), (2, 0.1, True))


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
    np.testing.assert_allclose(k, reference[key + "k"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(E, reference[key + "E"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(rate_shell, reference[key + "rate_shell"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(rapidd.dm_response_electron(m, jchi), reference[key + "dm_response"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x10['dRdS2'], reference[key + "X10_dRdS2"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x10['dRdS2_shell'], reference[key + "X10_dRdS2_shell"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x10['signal'], reference[key + "X10_signal"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x1t['dRdS2'], reference[key + "X1T_dRdS2"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x1t['signal'], reference[key + "X1T_signal"], rtol=RTOL, atol=0)
    np.testing.assert_allclose(x1t['n_bkg'], reference[key + "X1T_bkg"], rtol=RTOL, atol=0)
    assert x1t['n_obs'] == int(reference[key + "X1T_obs"])


def _numbers(path):
    '''All the numbers of a text file, in order. The numpy scalar wrappers written by repr with
    numpy >= 2 (e.g. np.float64(0.1) in result.txt) are removed first, so that the files written
    with numpy 1 and 2 give the same numbers.'''
    out = []
    with open(path) as f:
        for line in f:
            line = re.sub(r"np\.(float|int)\d*\(([^()]*)\)", r"\2", line)
            for tok in line.replace(',', ' ').replace('[', ' ').replace(']', ' ').split():
                try:
                    out.append(float(tok))
                except ValueError:
                    pass
    return np.array(out)


@pytest.mark.parametrize("dm_spin,m,modulation", DD_CASES)
def test_DDrate_electron_save(tmp_path, dm_spin, m, modulation):
    '''Same files as direct_electron_port, with the alphas of MadDM (smoke test of the MadDM entry point)'''
    name = "ddrate_spin%i_m%g%s" % (dm_spin, m, "_mod" if modulation else "")
    ref = os.path.join(REF_DIR, name)
    out = tmp_path / name
    out.mkdir()
    res = DDrate_electron_save(CARD, ALPHAS, dm_spin, M_E, m, str(out), halo_path=rapidd.halo_path,
                               modulation=modulation, days=[0, 100] if modulation else None)
    shutil.rmtree(out / "halo_table_modulation", ignore_errors=True)
    with open(out / "result.txt", "w") as f:
        for key in sorted(res):
            f.write("%s %r\n" % (key, res[key]))
    files = sorted(os.listdir(ref))
    assert sorted(os.listdir(out)) == files
    for fname in files:
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
