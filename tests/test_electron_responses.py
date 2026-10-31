'''DM responses R_1 ... R_4 of the C library against an independent transcription of
Eqs. (DM response 1-4) of arXiv:1912.08204 (main.tex lines 449-458).'''
import itertools

import numpy as np
import pytest

import rapidd

from conftest import GEV2_TO_CM2

M_E = rapidd.M_ELECTRON
Q_REF = 7.2973525643e-3 * M_E
OPS = (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
S_CHI_OPS = (4, 5, 6, 8, 9, 11, 12, 13, 14, 15)  # operators with S_chi
INTERFERING = ((4, 6), (12, 15), (13, 14))
# (q [GeV], dE [GeV], mchi [GeV])
KINEMATICS = ((3e-6, 15e-9, 0.01), (2e-5, 80e-9, 0.1), (1e-4, 300e-9, 1.0), (5e-4, 20e-9, 0.05))


def R_paper(c, jchi, q_me, vperp2, qv):
    '''Eqs. (DM response 1-4) of arXiv:1912.08204 for real couplings c[i] (i = 1..15), with
    q_me = |q|/m_e, vperp2 = (v_el^perp)^2 and qv = (q/m_e) . v_el^perp'''
    jj = jchi * (jchi + 1.)
    q2 = q_me**2
    R1 = (c[1]**2 + c[3]**2 / 4 * q2 * vperp2 - c[3]**2 / 4 * qv**2 + c[7]**2 / 4 * vperp2 + c[10]**2 / 4 * q2
          + jj / 12 * (3 * c[4]**2 + c[6]**2 * q2**2 + (4 * c[8]**2 + 2 * c[12]**2) * vperp2
                       + (2 * c[9]**2 + 4 * c[11]**2 + 2 * c[4] * c[6]) * q2
                       + (4 * c[5]**2 + c[13]**2 + c[14]**2 - 2 * c[12] * c[15]) * q2 * vperp2
                       + c[15]**2 * q2**2 * vperp2 - c[15]**2 * q2 * qv**2
                       + (-4 * c[5]**2 + 2 * c[13] * c[14] + 2 * c[12] * c[15]) * qv**2))
    R2 = qv * (-c[7]**2 / 2 / q2 - jj / 6 * ((4 * c[8]**2 + 2 * c[12]**2) / q2 + (c[13] + c[14])**2))
    R3 = (c[3]**2 / 4 * q2 + c[7]**2 / 4
          + jj / 12 * (4 * c[8]**2 + 2 * c[12]**2 + (4 * c[5]**2 + c[13]**2 + c[14]**2 - 2 * c[12] * c[15]) * q2
                       + c[15]**2 * q2**2))
    R4 = -c[3]**2 / 4 + jj / 12 * (-4 * c[5]**2 - c[15]**2 * q2 + 2 * c[12] * c[15] + 2 * c[13] * c[14])
    return np.array(np.broadcast_arrays(R1, R2, R3, R4))


def kinematics(q, dE, mchi, v):
    '''(v_el^perp)^2 and (q/m_e).v_el^perp at k = 0 from App. C of arXiv:1912.08204'''
    mu = mchi * M_E / (mchi + M_E)
    vperp2 = v**2 + q**2 / (4 * mu**2) * (mchi - M_E) / (mchi + M_E) - dE / mu
    qv = (dE - q**2 / (2 * M_E)) / M_E
    return vperp2, qv


def couplings(coeffs, q):
    '''c_i(q) = c_i^s + (q_ref/q)^2 c_i^l from {(op, interaction): value}'''
    c = np.zeros(16)
    for (op, inter), val in coeffs.items():
        c[op] += val * ((Q_REF / q)**2 if inter == "long-range" else 1.)
    return c


def R_rapidd(q, dE, mchi, jchi):
    out = [rapidd.dm_response_electron(mchi, jchi, resp=j, q=q, dE=dE, m_e=M_E, return_v2=True) for j in (1, 2, 3, 4)]
    return np.array([r[0] for r in out]) / GEV2_TO_CM2, out[0][1] / GEV2_TO_CM2


def set_coeffs(coeffs):
    rapidd.reset_coefficients()
    for (op, inter), val in coeffs.items():
        rapidd.set_any_Ecoeff(val, op, interaction=inter)


def check(coeffs, jchi):
    set_coeffs(coeffs)
    for q, dE, mchi in KINEMATICS:
        c = couplings(coeffs, q)
        R0, Rv2 = R_rapidd(q, dE, mchi, jchi)
        v = 7e-4  # any speed: R_1 is linear in v^2
        vp2, qv = kinematics(q, dE, mchi, v)
        ref = R_paper(c, jchi, q / M_E, vp2, qv)
        scale = np.max(np.abs(ref)) + 1e-300
        got = R0.copy()
        got[0] += Rv2 * v**2
        np.testing.assert_allclose(got, ref, rtol=1e-10, atol=1e-12 * scale)
        # v independent part of R_1 (v = 0 in the kinematics)
        vp2_0, _ = kinematics(q, dE, mchi, 0.)
        np.testing.assert_allclose(R0[0], R_paper(c, jchi, q / M_E, vp2_0, qv)[0], rtol=1e-10, atol=1e-12 * scale)
        # the coefficient of (v_el^perp)^2 in R_1 is R_3
        np.testing.assert_allclose(Rv2, R0[2], rtol=1e-12, atol=1e-300)


@pytest.mark.parametrize("jchi", (0.5, 0.0))
@pytest.mark.parametrize("interaction", ("contact", "long-range"))
@pytest.mark.parametrize("op", OPS)
def test_single_operator(op, interaction, jchi):
    check({(op, interaction): 0.7}, jchi)
    rapidd.reset_coefficients()


@pytest.mark.parametrize("jchi", (0.5, 0.0))
@pytest.mark.parametrize("ops", INTERFERING)
def test_interference(ops, jchi):
    for signs in ((1., 1.), (1., -1.)):
        check({(ops[0], "contact"): 0.8 * signs[0], (ops[1], "long-range"): 1.3 * signs[1]}, jchi)
        check({(ops[0], "contact"): 0.8 * signs[0], (ops[1], "contact"): 1.3 * signs[1]}, jchi)
    rapidd.reset_coefficients()


def test_all_operators_random():
    rng = np.random.default_rng(1912)
    for _ in range(5):
        coeffs = {(op, inter): rng.normal() for op in OPS for inter in ("contact", "long-range")}
        check(coeffs, 0.5)
    rapidd.reset_coefficients()


@pytest.mark.parametrize("op", OPS)
def test_long_range_at_q_ref(op):
    '''c^l at q = q_ref is the same as c^s'''
    for q_dE_m in ((Q_REF, 30e-9, 0.1), (Q_REF, 200e-9, 1.0)):
        q, dE, m = q_dE_m
        set_coeffs({(op, "contact"): 0.9})
        Rs = R_rapidd(q, dE, m, 0.5)
        set_coeffs({(op, "long-range"): 0.9})
        Rl = R_rapidd(q, dE, m, 0.5)
        np.testing.assert_allclose(Rl[0], Rs[0], rtol=1e-12)
        np.testing.assert_allclose(Rl[1], Rs[1], rtol=1e-12)
    rapidd.reset_coefficients()


@pytest.mark.parametrize("op", OPS)
def test_scalar_dm(op):
    '''With jchi = 0 the operators with S_chi do not contribute'''
    set_coeffs({(op, "contact"): 1.0, (op, "long-range"): 0.5})
    for q, dE, m in KINEMATICS:
        R0, Rv2 = R_rapidd(q, dE, m, 0.0)
        if op in S_CHI_OPS:
            assert np.all(R0 == 0.) and Rv2 == 0.
        else:
            assert np.any(R0 != 0.)
    rapidd.reset_coefficients()


def test_legacy_si_sd_value():
    '''MadDM: R_1 = (c_1^2 + jchi(jchi+1)/4 c_4^2) * gev_to_pb * pb_to_cm2 at q = 0'''
    set_coeffs({(1, "contact"): 3e-3, (4, "contact"): 7e-3})
    for jchi in (0., 0.5):
        assert rapidd.dm_response_electron(0.1, jchi) == (3e-3**2 + jchi * (jchi + 1.) / 4. * 7e-3**2) * 3.893793721e+8 * 1.e-36
        assert rapidd.dm_response_electron(0.1, jchi, resp=2) == 0.
    rapidd.reset_coefficients()
