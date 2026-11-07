'''DM-nucleus rate with (v_T^perp)^2 = v^2 - vmin^2 (arXiv:1308.6288, elastic, and
arXiv:1409.0536, inelastic) and complex NREFT coefficients (inelastic matching, Table I of
arXiv:1409.0536).

40Ar has spin 0, so only the M (and Phi'') responses are not zero: the pure O8 rate is
(j(j+1)/3) M <v^2 - vmin^2> and the pure O1 rate is M <1>, and their ratio depends only on the
halo integrals and on vmin.'''
import ctypes

import numpy as np
import pytest

import rapidd

C_KMS = 2.99792458e+5            # km/s, as in phys_consts.h
CONV_V0 = 1.698603e+14           # conv_factor_v0 of difrate_isotope_dEr
CONV_V2 = 1.889947e+3            # conv_factor_v2 of difrate_isotope_dEr
JCHI = 0.5                       # DM spin used by difrate_dER_python
A_AR, Z_AR = 40, 18
RHO = 0.3

_halo = rapidd.core._crapidd.halo
_halo.argtypes = [ctypes.c_double, ctypes.c_int]
_halo.restype = ctypes.c_double
_mass_nucleus = rapidd.core._crapidd.approx_mass_nucleus
_mass_nucleus.argtypes = [ctypes.c_double, ctypes.c_double]
_mass_nucleus.restype = ctypes.c_double


def vmin_kms(mchi, E, delta, A=A_AR, Z=Z_AR):
    '''vmin = |q/(2 mu_T) + delta/q| in km/s, E and delta in keV'''
    mT = _mass_nucleus(A, Z)
    mu = mchi * mT / (mchi + mT)
    q = np.sqrt(2. * mT * E * 1e-6)
    return C_KMS * abs(q / (2. * mu) + delta * 1e-6 / q)


def rate(mchi, E, target="Ar", basis="pn_Fitz", delta=0.):
    return float(rapidd.difrate_dER(RHO, mchi, E, model="SHM", target=target, basis=basis, delta=delta))


def set_op(op, cp, cn=None):
    rapidd.reset_coefficients()
    rapidd.set_any_Ncoeff(cp, op, "p")
    rapidd.set_any_Ncoeff(cp if cn is None else cn, op, "n")


CASES = [(mchi, delta, E) for mchi, delta in [(20., 0.), (100., 0.), (1000., 0.), (1000., 30.), (300., -20.)]
         for E in (1., 5., 15., 40.)]


@pytest.mark.parametrize("mchi, delta, E", CASES)
def test_O8_over_O1_is_vperp2(default_halo, mchi, delta, E):
    '''rate(O8)/rate(O1) = (j(j+1)/3) <v^2 - vmin^2>/c^2, with <.> = int f/v (.) / int f/v'''
    vm = vmin_kms(mchi, E, delta)
    eta0, eta2 = _halo(vm, 0), _halo(vm, 2)
    if eta0 <= 0.:
        pytest.skip("vmin above the maximal velocity")
    set_op(1, 1.)
    r1 = rate(mchi, E, delta=delta)
    set_op(8, 1.)
    r8 = rate(mchi, E, delta=delta)
    expected = JCHI * (JCHI + 1.) / 3. * max(eta2 - vm**2 * eta0, 0.) * CONV_V2 / (eta0 * CONV_V0)
    assert r8 / r1 == pytest.approx(expected, rel=1e-9, abs=1e-15)
    # the old v^2 (instead of (v_T^perp)^2) weight is larger, by a lot near threshold
    assert r8 / r1 < JCHI * (JCHI + 1.) / 3. * eta2 * CONV_V2 / (eta0 * CONV_V0)


@pytest.mark.parametrize("mchi, delta, E", [(1000., 30., 5.), (1000., 30., 15.), (300., -20., 5.), (300., 50., 40.)])
def test_lagrangian_13_inelastic(default_halo, mchi, delta, E):
    '''L13 with delta != 0: c8 O8 + i (m_N delta/q^2) c8 O11 (1409.0536). With only protons,
    rate/rate(O1 with c1 = c8) = (j(j+1)/3) [<v^2 - vmin^2> + delta^2/q^2]
    = (j(j+1)/3) <|v_el^perp|^2>, the elastic-basis result.'''
    m_M = 10.
    vm = vmin_kms(mchi, E, delta)
    eta0, eta2 = _halo(vm, 0), _halo(vm, 2)
    rapidd.reset_coefficients()
    rapidd.set_any_N_EFT_coeff(13, mchi, m_M, nucleon="p", delta=delta)
    c8 = rapidd.Cp_val(8)
    assert rapidd.Cp_val(11, q_power=-2) == pytest.approx(1j * rapidd.core.mproton * delta * 1e-6 * c8)
    r13 = rate(mchi, E, delta=delta)
    set_op(1, c8, 0.)
    r1 = rate(mchi, E, delta=delta)
    q2 = 2. * _mass_nucleus(A_AR, Z_AR) * E * 1e-6
    vperp2 = max(eta2 - vm**2 * eta0, 0.) * CONV_V2 / (eta0 * CONV_V0)
    assert r13 / r1 == pytest.approx(JCHI * (JCHI + 1.) / 3. * (vperp2 + (delta * 1e-6)**2 / q2), rel=1e-9)


@pytest.mark.parametrize("delta", [0., 40.])
def test_rate_depends_on_re_c_cstar(default_halo, delta):
    '''Global phases do not change the rate, coefficients with a relative phase i do not interfere
    (Xe, interfering pairs 1-3 and 4-6)'''
    mchi, E = 500., 10.

    def r(coeffs):
        rapidd.reset_coefficients()
        for op, (cp, cn) in coeffs.items():
            rapidd.set_any_Ncoeff(cp, op, "p")
            rapidd.set_any_Ncoeff(cn, op, "n")
        return rate(mchi, E, target="Xe", basis="iso_GCN5082", delta=delta)

    for i, j in [(1, 3), (4, 6)]:
        a, b = (0.7, -0.2), (1.3, 0.4)
        real = r({i: a, j: b})
        ph = np.exp(0.83j)
        assert r({i: (a[0] * ph, a[1] * ph), j: (b[0] * ph, b[1] * ph)}) == pytest.approx(real, rel=1e-12)
        no_interf = r({i: a}) + r({j: b})
        assert real != pytest.approx(no_interf, rel=1e-6)
        assert r({i: a, j: (1j * b[0], 1j * b[1])}) == pytest.approx(no_interf, rel=1e-12)


def test_complex_coefficients_set_and_reset():
    rapidd.reset_coefficients()
    rapidd.set_any_Ncoeff(1. + 2.j, 7, "p")
    assert rapidd.Cp_val(7) == 1. + 2.j and rapidd.Cn_val(7) == 0.
    rapidd.set_any_Ncoeff(3., 7, "p")              # a real value also clears the imaginary part
    assert rapidd.Cp_val(7) == 3. and isinstance(rapidd.Cp_val(7), float)
    rapidd.set_any_Ncoeff(-1.j, 7, "n", q_power=-2)
    rapidd.set_any_Ncoeff(0.5, 7, "n", q_power=2)
    assert rapidd.Cn_val(7) == 0. and rapidd.Cn_val(7, q_power=-2) == -1.j and rapidd.Cn_val(7, q_power=2) == 0.5
    with pytest.raises(ValueError):
        rapidd.set_any_Ncoeff(1., 7, "p", q_power=1)
    rapidd.reset_coefficients()
    assert all(rapidd.Cp_val(7, k) == 0. and rapidd.Cn_val(7, k) == 0. for k in (-2, 0, 2))


def test_lagrangian_19_inelastic_coefficients():
    '''c4 O4 + c6 O6 of L19 is i (delta/m_N) c14 [S_chi.S_N - (S_chi.q^)(S_N.q^)] (1409.0536)'''
    rapidd.reset_coefficients()
    mchi, m_M, delta = 100., 5., 25.
    rapidd.set_any_N_EFT_coeff(19, mchi, m_M, nucleon="n", delta=delta)
    c14, c4, c6 = rapidd.Cn_val(14), rapidd.Cn_val(4), rapidd.Cn_val(6, q_power=-2)
    mN = rapidd.core.mneutron
    assert c4 == pytest.approx(1j * delta * 1e-6 / mN * c14)
    # O6 = (S_chi.q/m_N)(S_N.q/m_N): c6(q) q^2/m_N^2 = -c4 at every q
    assert c6 / mN**2 == pytest.approx(-c4)
    rapidd.set_any_N_EFT_coeff(19, mchi, m_M, nucleon="n", delta=0.)
    assert rapidd.Cn_val(4) == 0. and rapidd.Cn_val(6, q_power=-2) == 0.
    rapidd.reset_coefficients()


@pytest.mark.parametrize("q_power", [-2, 2])
@pytest.mark.parametrize("E", [2., 20.])
def test_q_power_evaluated_at_isotope_q(default_halo, q_power, E):
    '''On a single isotope, c_k q^k (k = -2, 2) gives the same rate as the q independent
    coefficient c_k q^k evaluated at q^2 = 2 m_T E_R'''
    mchi, a = 100., 0.3 + 0.4j
    q2 = 2. * _mass_nucleus(A_AR, Z_AR) * E * 1e-6
    rapidd.reset_coefficients()
    for nuc in "pn":
        rapidd.set_any_Ncoeff(a, 1, nuc, q_power=q_power)
    r_q = rate(mchi, E)
    set_op(1, a * q2**(q_power // 2))
    assert r_q == pytest.approx(rate(mchi, E), rel=1e-12)


def test_eft_coefficients_independent_of_energy(default_halo):
    '''L6 (q^2 dependent c1): one call of set_any_N_EFT_coeff is valid at all energies, and the
    deprecated nucl and E arguments are ignored'''
    mchi, m_M = 100., 1.
    rapidd.reset_coefficients()
    for nuc in "pn":
        rapidd.set_any_N_EFT_coeff(6, mchi, m_M, nucleon=nuc)
    rates = [rate(mchi, E, target="Xe", basis="iso_GCN5082") for E in (2., 20.)]
    rapidd.reset_coefficients()
    for E, r in zip((2., 20.), rates):
        with pytest.warns(DeprecationWarning):
            for nuc in "pn":
                rapidd.set_any_N_EFT_coeff(6, mchi, m_M, "xe", E, nuc)
        assert rate(mchi, E, target="Xe", basis="iso_GCN5082") == pytest.approx(r, rel=1e-12)
    rapidd.reset_coefficients()
