import os

import numpy as np
import pytest

import rapidd

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(TESTS_DIR, "data")
PAPER_DIR = os.path.join(DATA_DIR, "paper_1912.08204")
PAPER_HALO = os.path.join(rapidd.core.base_dir, "..", "lib", "halo_table", "halo_table_SHM_1912.08204.dat")
GEV2_TO_CM2 = 3.893793721e+8 * 1.e-36  # unit factor of dm_response_electron


@pytest.fixture
def default_halo():
    '''Default RAPIDD halo (used by MadDM) and all the atomic responses'''
    rapidd.read_halo(rapidd.halo_path)
    rapidd.read_atomic_responses()
    rapidd.reset_coefficients()
    yield
    rapidd.reset_coefficients()


@pytest.fixture
def paper_halo():
    '''Halo of arXiv:1912.08204 (v0 = 220, vesc = 544, vE = 244 km/s, up to halo(vmin, 2))'''
    rapidd.read_halo(PAPER_HALO)
    rapidd.read_atomic_responses()
    rapidd.reset_coefficients()
    yield
    rapidd.reset_coefficients()
    rapidd.read_halo(rapidd.halo_path)


def load_csv(name):
    return np.loadtxt(os.path.join(PAPER_DIR, name + ".csv"), delimiter=",", comments="#")
