'''Reading of the nuclear responses (lib/nuclear_responses): the errors of read_nuclear_responses
keep the tables read before (and the rates unchanged).'''
import os
import shutil

import numpy as np
import pytest

import rapidd

E = np.array([1., 10., 50.])


def xenon_rate(basis):
    rapidd.read_halo(rapidd.halo_path)
    rapidd.reset_coefficients()
    rapidd.set_any_Ncoeff(1., 1, "p")
    rapidd.set_any_Ncoeff(0.5, 4, "n")
    rate = rapidd.difrate_dER(0.3, 100., E, target="Xe", basis=basis)
    rapidd.reset_coefficients()
    return rate


def corrupt(tmp_path, name, old, new):
    '''Copy of lib/nuclear_responses with the first occurrence of old replaced by new in name'''
    path = tmp_path / "nuclear_responses"
    shutil.copytree(rapidd.nuclear_resp_dir, path)
    text = (path / name).read_text()
    assert old in text
    (path / name).write_text(text.replace(old, new, 1))
    return str(path)


@pytest.mark.parametrize("basis", rapidd.VALID_NUCLEAR_MODELS)
def test_read_default(basis):
    ref = xenon_rate(basis)
    assert np.all(ref > 0)
    rapidd.read_nuclear_responses()
    assert np.array_equal(xenon_rate(basis), ref)


def test_missing_directory():
    ref = xenon_rate("iso_GCN5082")
    with pytest.raises(IOError, match="error 2"):
        rapidd.read_nuclear_responses("/nonexistent")
    assert np.array_equal(xenon_rate("iso_GCN5082"), ref)


@pytest.mark.parametrize("name,old,new", [
    ("iso_GCN5082.txt", "  128  54  M    M    +     +  ", "  128  54  M    M    p     +  "),  # pn label in iso
    ("pn_Fitz.txt", "   16   8  M    M    p   p  ", "   16   8  M    M    +   p  "),       # iso label in pn
    ("pn_SN100PN.txt", "\n  128  54  M", "\n  128  54  M  1"),                              # shifted columns
    ("iso_GCN5082_MEN.txt", "  128  54  M    +  ", "  128  54  M    +  1"),                  # extra column
    ("iso_Fitz.txt", "  128  54", "  200  54"),                                              # A > 199
    ("pn_Fitz.txt", "   16   8  M    M    n   n  ", "   16   8  M    M    p   p  "),       # repeated row
    ("iso_GCN5082_MEN.txt", "  128  54  M    -  ", "  128  54  M    +  "),                  # repeated row
])
def test_malformed(tmp_path, name, old, new):
    ref = xenon_rate("iso_Fitz")
    with pytest.raises(IOError, match="error 3"):
        rapidd.read_nuclear_responses(corrupt(tmp_path, name, old, new))
    assert np.array_equal(xenon_rate("iso_Fitz"), ref)
