'''Nuclear responses (lib/nuclear_responses): reading (the errors of read_nuclear_responses keep the tables read
before and the rates unchanged), physical checks of the tables of all the nuclear models at q = 0, and coverage
of the DM-nucleus targets by the nuclear models (check_nuclear_model).'''
import math
import os
import shutil

import numpy as np
import pytest

import rapidd

E = np.array([1., 10., 50.])

# models in the isospin basis (nucleon labels +, -; column c = j_N), the others are in the p, n basis
ISOSPIN_MODELS = ("Anand2013", "McElvain2022")
# tolerance of F_M(0) = Z^2, Z N, N^2: rounded default OBDMs of DMFormFactor (Anand2013), exact OBDMs
# (McElvain2022), digits of the papers (Fitz2012, Khaleq2023)
M0_TOL = {"Anand2013": 1e-3, "McElvain2022": 1e-5, "Fitz2012": 0.03, "Khaleq2023": 0.03}
# F_M(0) of arXiv:1203.3542 rounded to 2 digits beyond the 3% tolerance: 70Ge n n 1400 (N^2 = 1444, -3.05%),
# 74Ge p n and n p 1300 (Z N = 1344, -3.27%)
M0_TOL_EXCEPTIONS = {("Fitz2012", 70, 32, "n", "n"): 0.035, ("Fitz2012", 74, 32, "p", "n"): 0.035,
                     ("Fitz2012", 74, 32, "n", "p"): 0.035}
SPIN_RESPONSES = ("Sp", "Spp", "D", "Pp", "Pt")

# targets covered by each nuclear model (check_nuclear_model accepts them, rejects the others)
COVERAGE = {
    "McElvain2022": ("Xe", "Ge", "Ar", "F", "C", "O", "Na", "Al", "Si", "Ca", "I", "Cs", "CaWO4"),
    "Fitz2012": ("Xe", "Ge", "Ar", "F", "Na", "Ca", "I", "CaWO4"),
    "Anand2013": ("Xe", "Ge", "F", "Na", "Si", "I"),
    "Khaleq2023": ("Xe",),
}


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
    ref = xenon_rate("McElvain2022")
    with pytest.raises(IOError, match="error 2"):
        rapidd.read_nuclear_responses("/nonexistent")
    assert np.array_equal(xenon_rate("McElvain2022"), ref)


@pytest.mark.parametrize("name,old,new", [
    ("McElvain2022.txt", "  128  54  M    M    +     +  ", "  128  54  M    M    p     +  "),  # pn label in iso
    ("Fitz2012.txt", "   16   8  M    M    p   p  ", "   16   8  M    M    +   p  "),       # iso label in pn
    ("Khaleq2023.txt", "\n  128  54  M", "\n  128  54  M  1"),                              # shifted columns
    ("Anand2013.txt", "   19   9  M    M    +     +  ", "   19   9  M    M    +     +  1  "),  # extra column
    ("Anand2013.txt", "  128  54", "  200  54"),                                             # A > 199
    ("Fitz2012.txt", "   16   8  M    M    n   n  ", "   16   8  M    M    p   p  "),       # repeated row
    ("McElvain2022.txt", "   19   9  M    M    -     -  ", "   19   9  M    M    +     +  "),  # repeated row
])
def test_malformed(tmp_path, name, old, new):
    ref = xenon_rate("Anand2013")
    with pytest.raises(IOError, match="error 3"):
        rapidd.read_nuclear_responses(corrupt(tmp_path, name, old, new))
    assert np.array_equal(xenon_rate("Anand2013"), ref)


def responses_at_q0(model):
    '''{(A, Z): {(O1, O2, N1, N2): F(0)}} of lib/nuclear_responses/<model>.txt in the p, n basis
    (F(0) = norm a0; the isospin responses are converted: F^pp = (F^++ + F^-- + F^+- + F^-+)/4,
    F^nn = (F^++ + F^-- - F^+- - F^-+)/4, F^pn + F^np = (F^++ - F^--)/2), and the raw rows {(A, Z): {(O1, O2,
    N1, N2): F(0)}} in the basis of the file.'''
    raw = {}
    with open(os.path.join(rapidd.nuclear_resp_dir, model + ".txt")) as f:
        for line in f:
            s = line.split()
            if not s or s[0].startswith('#'):
                continue
            c = float(s[6])
            norm = 16 * math.pi / (2 * c + 1) if model in ISOSPIN_MODELS else c
            raw.setdefault((int(s[0]), int(s[1])), {})[(s[2], s[3], s[4], s[5])] = norm * float(s[7])
    if model not in ISOSPIN_MODELS:
        return raw, raw
    pn = {}
    for nucleus, rows in raw.items():
        out = pn.setdefault(nucleus, {})
        for o1, o2 in {(k[0], k[1]) for k in rows}:
            g = lambda t1, t2: rows.get((o1, o2, t1, t2), 0.)
            out[(o1, o2, "p", "p")] = (g("+", "+") + g("-", "-") + g("+", "-") + g("-", "+")) / 4
            out[(o1, o2, "n", "n")] = (g("+", "+") + g("-", "-") - g("+", "-") - g("-", "+")) / 4
            out[(o1, o2, "p", "n+n p")] = (g("+", "+") - g("-", "-")) / 2
    return pn, raw


@pytest.mark.parametrize("model", rapidd.VALID_NUCLEAR_MODELS)
def test_M_at_q0(model):
    '''F_M(0) = Z^2, N^2 and 2 Z N (F^pn + F^np) for every nucleus; in the isospin basis F^++ = A^2,
    F^-- = (Z - N)^2, F^+- = F^-+ = A (Z - N). Tolerance relative to A^2 (Z^2 etc. for the p, n tables).'''
    pn, raw = responses_at_q0(model)
    tol = M0_TOL[model]
    for (A, Z), rows in pn.items():
        N = A - Z
        if model in ISOSPIN_MODELS:
            r = raw[(A, Z)]
            for (t1, t2), ref in {("+", "+"): A**2, ("-", "-"): (Z - N)**2, ("+", "-"): A * (Z - N),
                                  ("-", "+"): A * (Z - N)}.items():
                assert r.get(("M", "M", t1, t2), 0.) == pytest.approx(ref, abs=tol * A**2), (A, Z, t1, t2)
            assert rows[("M", "M", "p", "p")] == pytest.approx(Z**2, abs=tol * A**2), (A, Z)
            assert rows[("M", "M", "n", "n")] == pytest.approx(N**2, abs=tol * A**2), (A, Z)
            assert rows[("M", "M", "p", "n+n p")] == pytest.approx(2 * Z * N, abs=tol * A**2), (A, Z)
        else:
            for (n1, n2), ref in {("p", "p"): Z**2, ("n", "n"): N**2, ("p", "n"): Z * N, ("n", "p"): Z * N}.items():
                t = M0_TOL_EXCEPTIONS.get((model, A, Z, n1, n2), tol)
                assert rows[("M", "M", n1, n2)] == pytest.approx(ref, rel=t), (A, Z, n1, n2)


@pytest.mark.parametrize("model", rapidd.VALID_NUCLEAR_MODELS)
def test_sigma_prime_is_twice_sigma_double_prime_at_q0(model):
    '''At q = 0 only the L = 1 multipole contributes, with two transverse (Sigma') and one longitudinal
    (Sigma'') components: F_Sigma'(0) = 2 F_Sigma''(0) for every nucleus and nucleon pair. DMFormFactor models
    (Anand2013, McElvain2022): to 1e-8; fits of the papers (Fitz2012, Khaleq2023): 5% for the values above 1e-3
    of the largest one of the nucleus (the coefficients of arXiv:1203.3542 have 2-3 digits).'''
    _, raw = responses_at_q0(model)
    checked = 0
    for nucleus, rows in raw.items():
        sp = {k[2:]: v for k, v in rows.items() if k[:2] == ("Sp", "Sp")}
        spp = {k[2:]: v for k, v in rows.items() if k[:2] == ("Spp", "Spp")}
        assert set(sp) == set(spp), nucleus
        scale = max([abs(v) for v in sp.values()] + [0.])
        for key, v in sp.items():
            if model in ISOSPIN_MODELS:
                assert v == pytest.approx(2. * spp[key], rel=1e-8, abs=1e-12 * scale), (nucleus, key)
                checked += 1
            elif abs(v) > 1e-3 * scale and abs(v) > 1e-3:
                assert v == pytest.approx(2. * spp[key], rel=0.05), (nucleus, key)
                checked += 1
    assert checked >= (8 if model != "Khaleq2023" else 4)   # Khaleq2023: Sigma' rows of 129Xe and 131Xe only


@pytest.mark.parametrize("model", rapidd.VALID_NUCLEAR_MODELS)
def test_no_spin_responses_for_even_even_nuclei(model):
    '''The even-even nuclei (j_N = 0) have no spin responses (Sigma', Sigma'', Delta, Phi', Phi-tilde'): their rows
    are omitted. The nuclei with odd A have the Sigma' rows.'''
    _, raw = responses_at_q0(model)
    for (A, Z), rows in raw.items():
        spin = [k for k in rows if k[0] in SPIN_RESPONSES or k[1] in SPIN_RESPONSES]
        if Z % 2 == 0 and (A - Z) % 2 == 0:
            assert not spin, ((A, Z), spin)
        else:
            assert any(k[:2] == ("Sp", "Sp") for k in spin), (A, Z)


def test_targets():
    '''Targets of the C table nuclear_targets: isotopes and CIAAW abundances'''
    t = rapidd.NUCLEAR_TARGET_ISOTOPES
    assert set(t) == {"Xe", "Ge", "Ar", "F", "C", "O", "Na", "Al", "Si", "Ca", "I", "Cs", "CaWO4"}
    assert [A for _, A, _ in t["Xe"]] == [128, 129, 130, 131, 132, 134, 136]
    assert t["Si"] == ((14, 28, 0.92223), (14, 29, 0.04685), (14, 30, 0.03092))
    assert t["C"] == ((6, 12, 0.9893), (6, 13, 0.0107))
    assert t["O"] == ((8, 16, 0.99757), (8, 18, 0.00205))
    assert t["Ca"] == ((20, 40, 0.96941),)
    for name in ("F", "Na", "Al", "I", "Cs", "Ar"):
        assert len(t[name]) == 1 and t[name][0][2] == 1.0
    for name, isotopes in t.items():
        assert sum(f for _, _, f in isotopes) <= 1. + 1e-12, name
        for Z, A, _ in isotopes:
            # measured masses (AME2020) of all the isotopes: within 0.2% of A u, not exactly A u (except 12C)
            m_u = (rapidd.nucleus_mass(A, Z) + Z * 0.51099895069e-3) / 0.93149410372
            assert m_u == pytest.approx(A, rel=2e-3), (A, Z)
            assert (A, Z) == (12, 6) or abs(m_u - A) > 1e-4, (A, Z)


@pytest.mark.parametrize("model", rapidd.VALID_NUCLEAR_MODELS)
@pytest.mark.parametrize("target", sorted(rapidd.NUCLEAR_TARGET_ISOTOPES))
def test_check_nuclear_model(target, model):
    '''check_nuclear_model accepts the targets whose isotopes (with nonzero abundance) are all in the model and
    rejects the others with a clear error; the accepted ones have a positive rate.'''
    if target in COVERAGE[model]:
        rapidd.check_nuclear_model(target, model)
        rapidd.read_halo(rapidd.halo_path)
        rapidd.reset_coefficients()
        rapidd.set_any_Ncoeff(1., 1, "p")
        rate = rapidd.difrate_dER(0.3, 100., 10., target=target, basis=model)
        rapidd.reset_coefficients()
        assert rate > 0
    else:
        with pytest.raises(ValueError, match="has no form factors"):
            rapidd.check_nuclear_model(target, model)


def test_default_model_covers_all_targets():
    for target in rapidd.NUCLEAR_TARGET_ISOTOPES:
        rapidd.check_nuclear_model(target, rapidd.DEFAULT_NUCLEAR_MODEL)


def test_unknown_target_or_model():
    with pytest.raises(ValueError, match="Unknown DM-nucleus target W"):
        rapidd.check_nuclear_model("W", "McElvain2022")
    with pytest.raises(ValueError, match="Fitz2012"):
        rapidd.check_nuclear_model("Xe", "Fitz")
