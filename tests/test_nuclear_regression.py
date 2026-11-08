'''Regression of DDrate_save (MadDM entry point for the DM-nucleus spectra), without and with the
annual modulation, against tests/data/regression_nuclear_modulation (see generate_reference.py).
Same structure as test_electron_regression.test_DDrate_electron_save.'''
import os
import shutil
import sys

import numpy as np
import pytest

import rapidd
from rapidd.calc_dRdE import DDrate_save

from conftest import DATA_DIR

REF_DIR = os.path.join(DATA_DIR, "regression_nuclear_modulation")
sys.path.insert(0, REF_DIR)
from generate_reference import CARD, ALPHAQ, DD_CASES, DAY_BINS, case_name  # noqa: E402
sys.path.pop(0)

RTOL = 1e-12


@pytest.mark.parametrize("m,modulation", DD_CASES)
def test_DDrate_save(tmp_path, m, modulation):
    name = case_name(m, modulation)
    ref = os.path.join(REF_DIR, name)
    out = tmp_path / name
    out.mkdir()
    DDrate_save(CARD, ALPHAQ, m, str(out), halo_path=rapidd.halo_path, modulation=modulation, day_bins=DAY_BINS)
    shutil.rmtree(out / "halo_table_modulation", ignore_errors=True)
    files = sorted(os.listdir(ref))
    assert sorted(os.listdir(out)) == files
    for fname in files:
        if fname.endswith(".npz"):
            a, b = np.load(out / fname), np.load(os.path.join(ref, fname))
            assert sorted(a.files) == sorted(b.files)
            for key in a.files:
                np.testing.assert_allclose(a[key], b[key], rtol=RTOL, atol=0, err_msg=key)
        else:
            np.testing.assert_allclose(np.loadtxt(out / fname), np.loadtxt(os.path.join(ref, fname)),
                                       rtol=RTOL, atol=0, err_msg=fname)
    rapidd.read_halo(rapidd.halo_path)
