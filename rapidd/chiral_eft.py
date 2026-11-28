"""
Spin-independent and spin-dependent DM-nucleus rates with the chiral EFT structure factors of

  SI: M. Hoferichter, P. Klos, J. Menendez, A. Schwenk, PRD 99, 055031 (2019), arXiv:1812.05617 (HKMS),
      with the xenon responses of PRD 94, 063505 (2016), arXiv:1605.08043. Fit coefficients from the
      ChiralEFT4DM notebook v1.0 (https://theorie.ikp.physik.tu-darmstadt.de/strongint/ChiralEFT4DM.html).
  SD: M. Hoferichter, J. Menendez, A. Schwenk, PRD 102, 074018 (2020), arXiv:2007.08529 (HMS), App. E
      (nominal SD structure factors of the LZ limits), and
      B. S. Hu et al., PRL 128, 072502 (2022), arXiv:2109.00193 (Hu2022): ab initio structure factors of three
      chiral interactions, whose spread (with the density of the two-body currents) is the nuclear uncertainty
      band of the LZ SD limits.

Equation numbers refer to the arXiv versions (HKMS v2, HMS v2, Hu2022 v2).

Only the nuclei of the RAPIDD targets are included: 19F, 40Ar, 70-76Ge, 128-136Xe (SI) and the odd isotopes
19F, 73Ge, 129Xe, 131Xe (SD; the other isotopes have J = 0).

Cross section for a spin-1/2 DM particle (HKMS Eq. 1 without the velocity-suppressed terms; HMS Eqs. 80-82, 86):

  dsigma/dq^2 = 1/(4 pi v^2) [ |F_SI(q)|^2 + |F_SD(q)|^2 ],

  F_SI = sum_{I=+,-} (c_I^M - q^2/m_N^2 cdot_I^M) F_I^M + c_pi F_pi + c_b F_b + q^2/(2 m_N^2) sum_I c_I^Phi'' F_I^Phi''
         (plus the L = 2 multipole of 131Xe, added incoherently),

  |F_SD|^2 = 4 pi/(2J+1) j_chi(j_chi+1)/12 S(q),
  S(q) = sum_L [ (a_p + a_n) F_+^{Sigma'_L} + ((a_p - a_n) + delta'(q) (a^A_p - a^A_n)) F_-^{Sigma'_L} ]^2
       + sum_L [ (a_p + a_n) F_+^{Sigma''_L} + ((a_p - a_n) + delta''(q) (a^A_p - a^A_n)) F_-^{Sigma''_L} ]^2,

with F_+- = F_p +- F_n, a_N = c4^A_N + c4^T_N the coefficients of O_4 = S_chi.S_N in the RAPIDD (NREFT)
normalization (sigma_SD^N = 3 mu_N^2 a_N^2 / (16 pi)) and c4^A (axial-vector, AA) and c4^T (tensor, TT) the two
contributions:
  - AA: HMS Eqs. 80-86 (Hu2022 Eqs. 2-3), with the corrections delta'(q) (transverse, Sigma') and delta''(q)
    (longitudinal, Sigma'') of HMS Eq. 86: axial radius, induced pseudoscalar (pion pole) and two-body currents
    delta a, delta a^P (HMS Eqs. 77-78, normal-ordering approximation). They act on the isovector coupling only.
  - TT: the same one-body transverse and longitudinal responses without delta' and delta'' (HMS Eqs. 108-109:
    bar S_ij = S_ij at delta = 0). HMS derive the tensor cross section for neutrinos and state that the
    two-body corrections of the tensor current would take a different form (they are not included); the
    nuclear side is the same for DM, while the coefficient is the WIMP-tensor matching, c4 = 8 sum_q C_q^TT g_T^q
    at leading order (O_4 only, no O_6 and no delta; Bishara et al., arXiv:1707.06998, App. A.7). The Hu2022
    structure factors include the two-body axial currents and apply to AA only: with sd_model = "Hu2022" a
    tensor coupling raises a ValueError.

The two-body SI term c_b F_b is written as in HKMS Eq. 1: F_b is the structure factor of HKMS Eq. 24
(normalized to F_b(0) = -2 E_b/M_pi, E_b the binding energy), c_b = f_pi^theta and c_pi = f_pi + 2 f_pi^theta
(HKMS Eq. 7 for a Dirac particle without the spin-2 couplings), with the pion couplings f_pi, f_pi^theta of HKMS
Eq. 8; the trace-anomaly response is F_pi^theta = 2 F_pi + F_b (HKMS Eq. 25). The fit files of F_b are called
*_wimps_eb_fit.dat in ChiralEFT4DM and B in
lib/nuclear_responses/chiral_eft/HKMS2019_SI.txt.

The DM-nucleus rate per unit target mass is dR/dE_R = rho/(2 pi m_chi) [|F_SI|^2 + |F_SD|^2] eta(v_min) (HKMS
App. D), with the halo table currently loaded in RAPIDD (read_halo). The couplings are in GeV^-2 (c^M = f_N of
RAPIDD's O_1 coefficient, sigma_SI^N = mu_N^2 c^2 / pi).

equivalent_sigma gives the DM-nucleon cross sections equivalent to a model in LZ (same expected number of events
as the reference model of each published LZ limit), to be compared directly with the published curves.

Not included: the velocity-suppressed terms xi_5, xi_8, xi_11 of HKMS Eq. 1 (O_5, O_8, O_11: for the SS, VV,
AA and TT couplings only O_5 from VV appears, suppressed by (q v / m_chi)^2 < 1e-6), the spin-2 couplings,
the coherent q^2/(m_chi m_N) pieces of the tensor coupling (HKMS Eq. 7) and the two-body currents of the tensor
current.
"""

import functools
import os
from dataclasses import dataclass, replace

import ctypes
import numpy as np

from .core import _crapidd, base_dir, nucleus_mass, mproton, mneutron, mnucleon
from .alpha_map import gev_to_pb, pb_to_cm2

########################################################################################
# Constants
########################################################################################
HBARC = 0.1973269804            # GeV fm
GEV2_FM2 = HBARC**2             # 1 fm^2 = 1/GEV2_FM2 GeV^-2
C_KMS = 299792.458              # speed of light [km/s]
# GeV/cm^3 1/GeV^5 (km/s)^-1 -> events/(kg day keV), the same factor of lib/source/difRateGen.c
CONV_RATE = 1.698603e+14
GEV2_TO_CM2 = gev_to_pb * pb_to_cm2   # GeV^-2 -> cm^2

# isospin-averaged masses of HKMS (Sec. IV) and HMS (Sec. III) [GeV]
M_N = 0.93892
M_PI = 0.138
M_P = 0.938272
M_NEUTRON = 0.939565
# nucleon mass of the harmonic-oscillator length of the ChiralEFT4DM fits [GeV]
M_N_OSC = 0.938919

# --- SI hadronic input of HKMS 2019 (Sec. II B and App. A; values of the ChiralEFT4DM notebook) ---
SIGMA_DOT = 0.27        # scalar radius sigma-dot [GeV^-1]
SIGMA_DOT_S = 0.3       # GeV^-1
XI_UD = 0.37            # (m_d - m_u)/(m_d + m_u), HKMS Eq. 9
F_PI_Q = {"u": 0.32, "d": 0.68}                       # f_q^pi = m_q/(m_u + m_d), HKMS Eq. 9 (f_s^pi ~ 0)
F_Q_N = {"p": {"u": 20.8e-3, "d": 41.1e-3, "s": 43e-3},  # scalar couplings f_q^N (ChiralEFT4DM)
         "n": {"u": 18.9e-3, "d": 45.1e-3, "s": 43e-3}}
F_Q_HEAVY = 68e-3                                     # f_Q^N per heavy quark (c, b, t)
KAPPA_P, KAPPA_N, KAPPA_S = 1.792847356, -1.91304272, 0.006     # kappa_s: HKMS App. A
R2_E_P, R2_E_N, R2_E_S = 0.7071 / GEV2_FM2, -0.1161 / GEV2_FM2, 0.0012 / GEV2_FM2   # GeV^-2 (r_E,s: HKMS App. A)

# --- axial-vector and tensor charges (HMS 2020, Table 1) ---
G_A = 1.27641
G_A_Q = {"p": {"u": 0.842, "d": -0.427, "s": -0.085}, "n": {"u": -0.427, "d": 0.842, "s": -0.085}}
G_T_Q = {"p": {"u": 0.784, "d": -0.204, "s": -0.0027}, "n": {"u": -0.204, "d": 0.784, "s": -0.0027}}
R2_A = 0.46 / GEV2_FM2          # axial radius <r_A^2> [GeV^-2] (HMS Sec. III)
F_PION = 0.09228                # pion decay constant [GeV] (HMS Sec. III)
G_PI_NN = np.sqrt(4. * np.pi * 13.7)    # g_piNN^2/(4 pi) = 13.7 (HMS Sec. III)

# --- low-energy constants of the two-body axial currents (HMS 2020, Table 5) ---
LEC_C1, LEC_C3, LEC_C4, LEC_C6 = -1.20, -4.45, 2.96, 5.01   # GeV^-1
LAMBDA_CHI = 0.7                                            # GeV

# (nuclear density [fm^-3], two-body correction delta a(0)) of the HMS band (Sec. III D): rho = 0.09 ... 0.11 fm^-3,
# 20% ... 30% reduction of the one-body isovector axial coupling at q = 0, c_D fixed by this condition. The extreme
# c_D = -6.08 (0.30) of HMS Table 5 is reached for rho = 0.09 fm^-3 and 30% (rho = 0.11 fm^-3 and 20%).
TWO_BODY_SETTINGS = {
    "central": (0.10, -0.25),
    "rho0.09_30": (0.09, -0.30),
    "rho0.09_20": (0.09, -0.20),
    "rho0.11_30": (0.11, -0.30),
    "rho0.11_20": (0.11, -0.20),
}
TWO_BODY_BAND = ("rho0.09_30", "rho0.09_20", "rho0.11_30", "rho0.11_20")

# --- Hu et al. 2022: the three chiral interactions and the densities of the two-body currents ---
SD_MODELS = ("HMS2020", "Hu2022")
BANDS = ("central", "min", "max")
HU2022_INTERACTIONS = ("EM1.8_2.0", "N4LO_LNL2", "DNNLOgo")
# Hu2022 (text after Eq. 3 and caption of Fig. 3): rho = 0.08 ... 0.12 fm^-3 (the notebook default is 0.09 ... 0.12).
# The band is the envelope over the three interactions and this grid of densities; the central value (not given by
# Hu2022) is defined here as the mean of the three interactions at rho = 0.10 fm^-3.
HU2022_RHO_BAND = (0.08, 0.09, 0.10, 0.11, 0.12)
HU2022_RHO_CENTRAL = 0.10
# constants of the Hu2022 notebook (ChEFT_SpinDependent/DM_constants.py of https://github.com/BaishanHu/AbInitio4DarkMatter,
# commit e10cccb), MeV and fm: they differ from those of HMS (e.g. g_A = 1.29, F_pi = 92.4 MeV, <r_A^2> = 0.47 fm^2)
HU2022_CONSTANTS = {"HBARC": 197.3269718, "M_PROTON": 938.2720813, "M_NEUTRON": 939.5654133, "gA": 1.29,
                    "gpinn": 13.1066, "Fpi": 92.4, "mpi": 138.04, "LambdaChi": 700.0, "r2A_fm2": 0.47}
# low-energy constants of the two-body currents of each interaction (cell 6 of DarkMatter_ChEFT_SD.ipynb): c1, c3, c4
# [GeV^-1], cD, and hat c6 (dimensionless, the same for the three interactions; it enters as (1 + hat c6)/m_N)
HU2022_LECS = {"EM1.8_2.0": {"c1": -0.81, "c3": -3.2, "c4": 5.4, "cD": 1.264, "c6": 3.7040},
               "N4LO_LNL2": {"c1": -1.1, "c3": -5.54, "c4": 4.17, "cD": -1.8, "c6": 3.7040},
               "DNNLOgo": {"c1": -0.74, "c3": -3.622246, "c4": 2.446123, "cD": 0.081, "c6": 3.7040}}

data_dir = os.path.join(base_dir, '..', 'lib', 'nuclear_responses', 'chiral_eft')

########################################################################################
# Couplings
########################################################################################


@dataclass(frozen=True)
class ChiralCouplings:
    """Couplings of the chiral EFT rate [GeV^-2], index 0 = isoscalar '+', 1 = isovector '-' for the SI ones
    (c_+- = (c_p +- c_n)/2) and 0 = proton, 1 = neutron for the SD ones.

    cM, cMdot, cPhi, cpi, cb : SI couplings c_+-^M, cdot_+-^M, c_+-^Phi'', c_pi, c_b of HKMS Eq. 1 (definitions in
               HKMS Eq. 7; c_b multiplies the two-body structure factor F_b)
    c4A, c4T : coefficients of O_4 = S_chi.S_N (RAPIDD normalization) from the axial-vector (AA, with the
               corrections delta', delta'') and tensor (TT, without them) quark currents
    """
    cM: tuple = (0., 0.)
    cMdot: tuple = (0., 0.)
    cPhi: tuple = (0., 0.)
    cpi: float = 0.
    cb: float = 0.
    c4A: tuple = (0., 0.)
    c4T: tuple = (0., 0.)

    def __add__(self, other):
        add2 = lambda a, b: (a[0] + b[0], a[1] + b[1])
        return ChiralCouplings(add2(self.cM, other.cM), add2(self.cMdot, other.cMdot), add2(self.cPhi, other.cPhi),
                               self.cpi + other.cpi, self.cb + other.cb, add2(self.c4A, other.c4A),
                               add2(self.c4T, other.c4T))

    def si_only(self):
        return replace(self, c4A=(0., 0.), c4T=(0., 0.))

    def sd_only(self):
        return ChiralCouplings(c4A=self.c4A, c4T=self.c4T)

    def without_two_body(self):
        """The same couplings without the two-body SI couplings c_pi and c_b (one-body SI response only)"""
        return replace(self, cpi=0., cb=0.)

    @property
    def has_si(self):
        return any(x != 0. for x in self.cM + self.cMdot + self.cPhi + (self.cpi, self.cb))

    @property
    def has_sd(self):
        return any(x != 0. for x in self.c4A + self.c4T)

    @property
    def has_tensor(self):
        return any(x != 0. for x in self.c4T)


def iso(cp, cn):
    """(c_p, c_n) -> (c_+, c_-) = ((c_p + c_n)/2, (c_p - c_n)/2)"""
    return (0.5 * (cp + cn), 0.5 * (cp - cn))


def couplings_standard_SI(c_p, c_n=None):
    """Standard SI coupling: only c^M (one-body, no radius, Phi'' or two-body terms), c = O_1 coefficient [GeV^-2]"""
    return ChiralCouplings(cM=iso(c_p, c_p if c_n is None else c_n))


def couplings_scalar(f_p, f_n, fdot_p=0., fdot_n=0., f_pi=0., f_pi_theta=0.):
    """Scalar (SS) couplings of HKMS Eq. 7 (Dirac DM, zeta = 1, no spin-2 couplings) from the nucleon and pion
    matrix elements of HKMS Eqs. 8 and A7: f_N [GeV^-2], fdot_N [GeV^-4], f_pi, f_pi^theta [GeV^-2].
    c_pi = f_pi + 2 f_pi^theta, c_b = f_pi^theta."""
    return ChiralCouplings(cM=iso(f_p, f_n), cMdot=tuple(M_N**2 * x for x in iso(fdot_p, fdot_n)),
                           cpi=f_pi + 2. * f_pi_theta, cb=f_pi_theta)


def couplings_vector(f1_p, f1_n, f2_p, f2_n, f1dot_p, f1dot_n, mchi):
    """Vector (VV) couplings of HKMS Eq. 7 (Dirac DM) from the nucleon matrix elements f_1^V, f_2^V [GeV^-2] and
    fdot_1^V [GeV^-4] of HKMS Eq. A7"""
    mu_N = mchi * M_N / (mchi + M_N)
    f1, f2, f1dot = iso(f1_p, f1_n), iso(f2_p, f2_n), iso(f1dot_p, f1dot_n)
    return ChiralCouplings(cM=f1,
                           cMdot=tuple(M_N**2 * (f1dot[I] + f2[I] / (4. * M_N**2)) for I in (0, 1)),
                           cPhi=tuple(f2[I] + 0.5 * (1. + mu_N / mchi) * f1[I] for I in (0, 1)))


def vector_nucleon_couplings(C, V_p=None, V_n=None):
    """f_1^V, f_2^V, fdot_1^V of the proton and of the neutron (HKMS Eq. A7) from the vector Wilson
    coefficients C = {'u': C_u, 'd': C_d, 's': C_s} [GeV^-2]. V_p, V_n: quark vector charges of the nucleons
    (default {'u': 2, 'd': 1} and {'u': 1, 'd': 2}); the strange vector charge is zero."""
    V_p = V_p or {"u": 2., "d": 1.}
    V_n = V_n or {"u": 1., "d": 2.}
    f1p = sum(C.get(q, 0.) * V_p.get(q, 0.) for q in ("u", "d", "s"))
    f1n = sum(C.get(q, 0.) * V_n.get(q, 0.) for q in ("u", "d", "s"))
    Cs = C.get("u", 0.) + C.get("d", 0.) + C.get("s", 0.)
    # F_2^{q,N} = V^p_q kappa_N + V^n_q kappa_N' + kappa_s: f_2^{V,p} = f1p kappa_p + f1n kappa_n + kappa_s sum_q C_q
    f2p = f1p * KAPPA_P + f1n * KAPPA_N + Cs * KAPPA_S
    f2n = f1n * KAPPA_P + f1p * KAPPA_N + Cs * KAPPA_S
    rp = R2_E_P / 6. - KAPPA_P / (4. * M_P**2)
    rn = R2_E_N / 6. - KAPPA_N / (4. * M_NEUTRON**2)
    rs = R2_E_S / 6. - KAPPA_S / (4. * M_N**2)
    f1dotp = f1p * rp + f1n * rn + Cs * rs
    f1dotn = f1n * rp + f1p * rn + Cs * rs
    return (f1p, f1n), (f2p, f2n), (f1dotp, f1dotn)


def couplings_from_quarks(mchi, C_SS=None, C_VV=None, C_AA=None, C_TT=None):
    """
    Couplings of a Dirac DM particle from the Wilson coefficients [GeV^-2, GeV^-3 for SS] of HKMS Eq. 3
    (Lambda = 1):

      L = sum_q [ C_q^SS m_q chibar chi qbar q + C_q^VV chibar gamma^mu chi qbar gamma_mu q
                + C_q^AA chibar gamma^mu gamma_5 chi qbar gamma_mu gamma_5 q + C_q^TT chibar sigma^munu chi qbar sigma_munu q ]

    Each argument is a dict {'u': ..., 'd': ..., 's': ..., ('c', 'b', 't' for SS)}. Hadronic input of HKMS 2019
    (SI) and HMS 2020 Table 1 (axial-vector and tensor charges). AA: O_4 coefficient c4 = -4 g_A^N with
    g_A^N = sum_q C_q^AA Delta q^N (HKMS Eq. A5); TT: c4 = 8 g_T^N with g_T^N = sum_q C_q^TT g_T^{q,N}.
    The heavy quarks enter through C'_g = -sum_Q C_Q/(12 pi) (HKMS Eq. 4): f_N gets m_N f_Q sum_Q C_Q (HKMS Eq. A7),
    f_pi gets -M_pi 2/27 sum_Q C_Q and f_pi^theta = M_pi 2/27 sum_Q C_Q (HKMS Eq. 8).
    """
    c = ChiralCouplings()
    if C_SS:
        f = {N: M_N * (sum(C_SS.get(q, 0.) * F_Q_N[N][q] for q in ("u", "d", "s"))
                       + F_Q_HEAVY * sum(C_SS.get(Q, 0.) for Q in ("c", "b", "t"))) for N in ("p", "n")}
        fdot = (C_SS.get("u", 0.) * (1. - XI_UD) / 2. * SIGMA_DOT + C_SS.get("d", 0.) * (1. + XI_UD) / 2. * SIGMA_DOT
                + C_SS.get("s", 0.) * SIGMA_DOT_S)
        heavy = sum(C_SS.get(Q, 0.) for Q in ("c", "b", "t"))
        f_pi = M_PI * (sum(C_SS.get(q, 0.) * F_PI_Q[q] for q in ("u", "d")) - 2. / 27. * heavy)
        f_pi_theta = M_PI * 2. / 27. * heavy
        c = c + couplings_scalar(f["p"], f["n"], fdot, fdot, f_pi, f_pi_theta)
    if C_VV:
        f1, f2, f1dot = vector_nucleon_couplings(C_VV)
        c = c + couplings_vector(*f1, *f2, *f1dot, mchi)
    if C_AA:
        gA = [sum(C_AA.get(q, 0.) * G_A_Q[N][q] for q in ("u", "d", "s")) for N in ("p", "n")]
        c = c + ChiralCouplings(c4A=(-4. * gA[0], -4. * gA[1]))
    if C_TT:
        gT = [sum(C_TT.get(q, 0.) * G_T_Q[N][q] for q in ("u", "d", "s")) for N in ("p", "n")]
        c = c + ChiralCouplings(c4T=(8. * gT[0], 8. * gT[1]))
    return c


########################################################################################
# Nuclear responses
########################################################################################

# RAPIDD targets: Z, mass numbers and mole fractions (the same as lib/source/phys_consts.h)
TARGETS = {
    "Xe": (54, (128, 129, 130, 131, 132, 134, 136), (0.01910, 0.26401, 0.04071, 0.21232, 0.26909, 0.10436, 0.08857)),
    "Ge": (32, (70, 72, 73, 74, 76), (0.2052, 0.2745, 0.0776, 0.3652, 0.0775)),
    "Ar": (18, (40,), (1.0,)),
    "F": (9, (19,), (1.0,)),
}
_SI_FILE_TAG = {"Xe": "xe", "Ge": "ge", "Ar": "ar", "F": "f"}
SI_RESPONSES = ("m+", "m-", "phi+", "phi-", "pi", "eb")


def oscillator_length(A):
    """Harmonic-oscillator length b [GeV^-1], hbar omega = (45 A^-1/3 - 25 A^-2/3) MeV (HKMS Eq. C15, HMS App. E)"""
    hw = (45. * A**(-1. / 3.) - 25. * A**(-2. / 3.)) * 1e-3
    return 1. / np.sqrt(M_N_OSC * hw)


def _fit(c, u):
    """Fit function exp(-u/2) sum_i c_i u^i (HKMS Eq. C14, HMS Eq. E1, Hu2022 text before Fig. 3)"""
    return np.exp(-u / 2.) * np.polynomial.polynomial.polyval(u, c)


def _read_table(name):
    """Rows (lists of strings) of the text table lib/nuclear_responses/chiral_eft/<name>, without the comments"""
    with open(os.path.join(data_dir, name)) as f:
        return [line.split() for line in f if line.strip() and not line.lstrip().startswith('#')]


# response names of HKMS2019_SI.txt (O, tau) -> keys of si_responses
_SI_KEYS = {("M", "+"): "m+", ("M", "-"): "m-", ("Ppp", "+"): "phi+", ("Ppp", "-"): "phi-", ("Pi", "0"): "pi",
            ("B", "0"): "eb"}
_SD_KEYS = {"Sp": "Sigma'", "Spp": "Sigma''"}


@functools.lru_cache(maxsize=None)
def _si_fits(target):
    """{response: {column: coefficients}} of the HKMS fits of target (HKMS2019_SI.txt); the columns are the mass
    numbers, and '131L2' is the L = 2 multipole of 131Xe"""
    Z = TARGETS[target][0]
    out = {r: {} for r in SI_RESPONSES}
    for row in _read_table('HKMS2019_SI.txt'):
        A, Z_, O, tau, L = int(row[0]), int(row[1]), row[2], row[3], int(row[4])
        if Z_ == Z:
            out[_SI_KEYS[(O, tau)]]["%d" % A + ("L%d" % L if L else "")] = np.array(row[5:], dtype=float)
    return out


@functools.lru_cache(maxsize=None)
def _sd_data(sd_model="HMS2020"):
    """{isotope: entry} of the SD tables HMS2020_SD.txt, Hu2022_SD.txt and spin_expectation.txt: entry = {'A', 'Z',
    'J', 'b_fm', and 'S_p', 'S_n', "Sigma'", "Sigma''" ({nucleon or isospin: {L: coefficients}}) for HMS2020, or
    'interactions' ({interaction: {'S_p', 'S_n', "Sigma'", "Sigma''"}}) for Hu2022}"""
    if sd_model not in SD_MODELS:
        raise ValueError("sd_model must be one of %s" % (SD_MODELS,))
    names = {"HMS2020": "HMS2020_SD.txt", "Hu2022": "Hu2022_SD.txt"}
    out = {}
    for row in _read_table(names[sd_model]):
        A, Z, O, N, L = int(row[0]), int(row[1]), _SD_KEYS[row[2]], row[3], row[4]
        if sd_model == "Hu2022":
            inte, J, b, coeffs = row[5], float(row[6]), float(row[7]), row[8:]
        else:
            inte, J, b, coeffs = None, float(row[5]), float(row[6]), row[7:]
        name = "%d%s" % (A, {9: "F", 32: "Ge", 54: "Xe"}[Z])
        e = out.setdefault(name, {"A": A, "Z": Z, "J": J, "b_fm": b})
        if inte is not None:
            e = e.setdefault("interactions", {}).setdefault(inte, {})
        e.setdefault(O, {}).setdefault(N, {})[L] = [float(x) for x in coeffs]
    for row in _read_table('spin_expectation.txt'):
        A, Z, model, inte = int(row[0]), int(row[1]), row[2], row[3]
        if model != sd_model:
            continue
        e = next(v for v in out.values() if v["A"] == A and v["Z"] == Z)
        if sd_model == "Hu2022":
            e = e["interactions"][inte]
        e["S_p"], e["S_n"] = float(row[4]), float(row[5])
    return out


def _sd_isotope_data(sd_model="HMS2020"):
    return _sd_data(sd_model)


def sd_isotopes(sd_model="HMS2020"):
    """Isotopes with SD responses: {name: (A, Z, J)}"""
    return {k: (v["A"], v["Z"], v["J"]) for k, v in _sd_isotope_data(sd_model).items()}


def _sd_entry(A, Z, sd_model="HMS2020"):
    for v in _sd_isotope_data(sd_model).values():
        if v["A"] == A and v["Z"] == Z:
            return v
    return None


def si_responses(target, A, q):
    """HKMS responses F_+^M, F_-^M, F_+^Phi'', F_-^Phi'', F_pi, F_b of the isotope A of target at q [GeV]
    (dict of arrays; 'L2' adds the L = 2 multipoles of 131Xe, when available)"""
    fits = _si_fits(target)
    u = 0.5 * (np.asarray(q, dtype=float) * oscillator_length(A))**2
    out = {r: _fit(fits[r][str(A)], u) for r in SI_RESPONSES}
    l2 = "%dL2" % A
    if l2 in fits["m+"]:
        out["L2"] = {r: _fit(fits[r][l2], u) for r in ("m+", "m-", "phi+", "phi-")}
    return out


def si_amplitude(target, A, q, c):
    """F_SI(q) of HKMS Eq. 1 [GeV^-2] and the L = 2 amplitude of 131Xe (0 otherwise)"""
    q = np.asarray(q, dtype=float)
    F = si_responses(target, A, q)
    x = q**2 / M_N**2

    def one_body(R):
        return ((c.cM[0] - x * c.cMdot[0]) * R["m+"] + (c.cM[1] - x * c.cMdot[1]) * R["m-"]
                + 0.5 * x * (c.cPhi[0] * R["phi+"] + c.cPhi[1] * R["phi-"]))

    amp = one_body(F) + c.cpi * F["pi"] + c.cb * F["eb"]
    amp_L2 = one_body(F["L2"]) if "L2" in F else np.zeros_like(q)
    return amp, amp_L2


def structure_factor_SI(target, A, q, c):
    """|F_SI(q)|^2 [GeV^-4], 4 pi v^2 dsigma/dq^2 of the SI part"""
    amp, amp_L2 = si_amplitude(target, A, q, c)
    return amp**2 + amp_L2**2


def sd_responses(A, Z, q):
    """One-body HMS responses {(response, nucleon, L): F(q)} of the isotope (A, Z), or None if J = 0 / not tabulated"""
    entry = _sd_entry(A, Z)
    if entry is None:
        return None
    b = entry["b_fm"] / HBARC
    u = 0.5 * (np.asarray(q, dtype=float) * b)**2
    return {(r, N, L): _fit(np.array(entry[r][N][L]), u)
            for r in ("Sigma'", "Sigma''") for N in ("p", "n") for L in entry[r][N]}


def sd_responses_hu(A, Z, q, interaction):
    """One-body Hu2022 responses {(response, '+' or '-', L): F_+-(q)} (F_+- = F_p +- F_n) of the isotope (A, Z) for
    one of HU2022_INTERACTIONS, or None if not tabulated. b = sqrt(hbar/(m_N omega)) with the hbar omega = 16 MeV
    of the IMSRG basis (column b of Hu2022_SD.txt)."""
    entry = _sd_entry(A, Z, "Hu2022")
    if entry is None:
        return None
    if interaction not in entry["interactions"]:
        raise ValueError("Hu2022 interaction must be one of %s" % (HU2022_INTERACTIONS,))
    data = entry["interactions"][interaction]
    u = 0.5 * (np.asarray(q, dtype=float) * entry["b_fm"] / HBARC)**2
    return {(r, I, L): _fit(np.array(data[r][I][L]), u)
            for r in ("Sigma'", "Sigma''") for I in ("+", "-") for L in data[r][I]}


def spin_expectation(A, Z, sd_model="HMS2020", interaction=None):
    """Shell-model <S_p>, <S_n> of HMS Table 4 (sd_model = "HMS2020") or the ab initio values of the Hu2022 notebook
    for one of HU2022_INTERACTIONS"""
    e = _sd_entry(A, Z, sd_model)
    if not e:
        return (0., 0.)
    if sd_model == "Hu2022":
        e = e["interactions"][interaction]
    return (e["S_p"], e["S_n"])


def nuclear_spin(A, Z):
    e = _sd_entry(A, Z)
    return e["J"] if e else 0.


########################################################################################
# Two-body currents and form-factor corrections of the axial-vector SD response (HMS Eqs. 77-79, 86)
########################################################################################
_GL_X, _GL_W = np.polynomial.legendre.leggauss(64)


def fermi_integrals(rho, Q, m_pi=M_PI, hbarc=HBARC):
    """I_1^sigma, I_2^sigma, I^P, I_c6 (Klos et al., arXiv:1304.7684, App. B; HMS Eq. 79) at density rho [fm^-3] and
    momentum Q [GeV], by Gauss-Legendre quadrature of their integral representations (k_F = (3 pi^2 rho/2)^1/3)"""
    Q = np.atleast_1d(np.asarray(Q, dtype=float))[:, None, None]
    kF = (1.5 * np.pi**2 * rho * hbarc**3)**(1. / 3.)
    p = 0.5 * kF * (_GL_X + 1.)[None, :, None]
    wp = 0.5 * kF * _GL_W[None, :, None]
    c = _GL_X[None, None, :]
    wc = _GL_W[None, None, :]
    den = m_pi**2 + p**2 + Q**2 / 4. - p * Q * c
    w = wp * wc / den
    I1 = 9. / (4. * kF**3) * np.sum(w * p**4 * (1. - c**2), axis=(1, 2))
    I2 = 3. / (2. * kF**3) * np.sum(w * (p**4 + p**2 * Q**2 / 4. - p**3 * Q * c), axis=(1, 2))
    IP = 9. / (8. * kF**3) * np.sum(w * (p**4 * (6. * c**2 - 2.) - 4. * p**3 * Q * c + p**2 * Q**2), axis=(1, 2))
    Ic6 = 9. / (2. * kF**3) * np.sum(w * (p**3 * Q * c - p**2 * Q**2 / 2.), axis=(1, 2))
    return I1, I2, IP, Ic6


def c_D_from_reduction(rho, delta_a0):
    """Contact coupling c_D such that the two-body currents give delta a(q = 0) = delta_a0 at density rho [fm^-3]"""
    I1, I2, _, _ = fermi_integrals(rho, 0.)
    rho_gev = rho * HBARC**3
    A = LEC_C4 / 3. * (3. * I2[0] - I1[0]) - (LEC_C3 - 1. / (4. * M_N)) / 3. * I1[0]
    return (A + delta_a0 * F_PION**2 / rho_gev) * 4. * G_A * LAMBDA_CHI


def delta_a(q, rho, c_D):
    """Two-body corrections delta a(q^2) and delta a^P(q^2) of HMS Eqs. 77-78 (P = 0 normal-ordering approximation)"""
    q = np.atleast_1d(np.asarray(q, dtype=float))
    I1, I2, IP, Ic6 = fermi_integrals(rho, q)
    rho_gev = rho * HBARC**3
    q2, m2 = q**2, M_PI**2
    cD = c_D / (4. * G_A * LAMBDA_CHI)
    da = -rho_gev / F_PION**2 * (LEC_C4 / 3. * (3. * I2 - I1) - (LEC_C3 - 1. / (4. * M_N)) / 3. * I1
                                 - LEC_C6 / 12. * Ic6 - cD)
    daP = rho_gev / F_PION**2 * (-2. * (LEC_C3 - 2. * LEC_C1) * m2 * q2 / (m2 + q2)**2
                                 + (LEC_C3 + LEC_C4 - 1. / (4. * M_N)) / 3. * IP
                                 - (LEC_C6 / 12. - 2. / 3. * LEC_C1 * m2 / (m2 + q2)) * Ic6
                                 - q2 / (m2 + q2) * (LEC_C3 / 3. * (I1 + IP) + LEC_C4 / 3. * (I1 + IP - 3. * I2))
                                 - cD * q2 / (m2 + q2))
    return da, daP


_Q_GRID = np.linspace(0., 0.6, 601)


@functools.lru_cache(maxsize=None)
def _delta_table(rho, delta_a0):
    da, daP = delta_a(_Q_GRID, rho, c_D_from_reduction(rho, delta_a0))
    q2 = _Q_GRID**2
    d1 = -q2 * R2_A / 6. + da
    d2 = -G_PI_NN * F_PION / (G_A * M_N) * q2 / (q2 + M_PI**2) + da + daP
    return d1, d2


def _resolve_two_body(two_body):
    if two_body is None or two_body == "none":
        return None
    if isinstance(two_body, str):
        if two_body not in TWO_BODY_SETTINGS:
            raise ValueError("two_body must be 'none', one of %s or a tuple (rho [fm^-3], delta a(0))"
                             % list(TWO_BODY_SETTINGS))
        return TWO_BODY_SETTINGS[two_body]
    rho, da0 = two_body
    return float(rho), float(da0)


def sd_deltas(q, two_body="central"):
    """delta'(q^2) (transverse) and delta''(q^2) (longitudinal) of HMS Eq. 86: axial radius, induced
    pseudoscalar and two-body currents. two_body: 'none' (0, 0), a key of TWO_BODY_SETTINGS or (rho, delta a(0))"""
    q = np.asarray(q, dtype=float)
    setting = _resolve_two_body(two_body)
    if setting is None:
        return np.zeros_like(q), np.zeros_like(q)
    d1, d2 = _delta_table(*setting)
    if np.any(q > _Q_GRID[-1]):
        da, daP = delta_a(q, setting[0], c_D_from_reduction(*setting))
        q2 = q**2
        return (-q2 * R2_A / 6. + da,
                -G_PI_NN * F_PION / (G_A * M_N) * q2 / (q2 + M_PI**2) + da + daP)
    return np.interp(q, _Q_GRID, d1), np.interp(q, _Q_GRID, d2)


def hu_delta_a(q, interaction, rho):
    """delta a(q^2), delta a^P(q^2) of Hu2022 as implemented in their notebook (DM_func.delta_2BC; the expressions of
    HMS Eqs. 77-78 with c_4 -> c_4 + 1/(4 m_N), c_6 -> (1 + hat c_6)/m_N) with the constants of DM_constants.py and
    the low-energy constants of each interaction (HU2022_CONSTANTS, HU2022_LECS)"""
    K, L = HU2022_CONSTANTS, HU2022_LECS[interaction]
    hbarc, mN = K["HBARC"] * 1e-3, 0.5 * (K["M_PROTON"] + K["M_NEUTRON"]) * 1e-3
    mpi, Fpi, gA, Lchi = K["mpi"] * 1e-3, K["Fpi"] * 1e-3, K["gA"], K["LambdaChi"] * 1e-3
    c1, c3, c4, cD, c6 = L["c1"], L["c3"], L["c4"], L["cD"], L["c6"]
    q = np.atleast_1d(np.asarray(q, dtype=float))
    I1, I2, IP, Ic6 = fermi_integrals(rho, q, m_pi=mpi, hbarc=hbarc)
    rho_gev = rho * hbarc**3
    q2, m2 = q**2, mpi**2
    c4t = c4 + 1. / (4. * mN)
    c6t = (1. + c6) / (12. * mN)
    cDt = cD / (4. * gA * Lchi)
    da = -rho_gev / Fpi**2 * (c4t / 3. * (3. * I2 - I1) + (-c3 + 1. / (4. * mN)) / 3. * I1 - c6t * Ic6 - cDt)
    daP = rho_gev / Fpi**2 * (-2. * (c3 - 2. * c1) * q2 * m2 / (m2 + q2)**2 + (c3 + c4) / 3. * IP
                              - (c6t * Ic6 - 2. / 3. * c1 * m2 / (m2 + q2) * Ic6)
                              - q2 / (m2 + q2) * (c3 / 3. * (I1 + IP) + c4t / 3. * (I1 + IP - 3. * I2))
                              - q2 / (m2 + q2) * cDt)
    return da, daP


def _hu_corrections(q, interaction, rho):
    """delta', delta'' of Hu2022 Eq. 3 (rho = None: 0, 0) with the constants of their notebook"""
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if rho is None:
        return np.zeros_like(q), np.zeros_like(q)
    K = HU2022_CONSTANTS
    mN, mpi = 0.5 * (K["M_PROTON"] + K["M_NEUTRON"]) * 1e-3, K["mpi"] * 1e-3
    r2A = K["r2A_fm2"] / (K["HBARC"] * 1e-3)**2
    da, daP = hu_delta_a(q, interaction, rho)
    q2 = q**2
    return (-q2 * r2A / 6. + da,
            -K["gpinn"] * K["Fpi"] * 1e-3 / (mN * K["gA"]) * q2 / (q2 + mpi**2) + da + daP)


@functools.lru_cache(maxsize=None)
def _hu_delta_table(interaction, rho):
    return _hu_corrections(_Q_GRID, interaction, rho)


def hu_deltas(q, interaction, rho):
    """delta'(q^2), delta''(q^2) of Hu2022 Eq. 3 for one of HU2022_INTERACTIONS at density rho [fm^-3]
    (rho = None: no corrections, the one-body structure factors)"""
    q = np.asarray(q, dtype=float)
    if rho is None:
        return np.zeros_like(q), np.zeros_like(q)
    if np.any(q > _Q_GRID[-1]):
        d1, d2 = _hu_corrections(q.ravel(), interaction, rho)
        return d1.reshape(q.shape), d2.reshape(q.shape)
    d1, d2 = _hu_delta_table(interaction, float(rho))
    return np.interp(q, _Q_GRID, d1), np.interp(q, _Q_GRID, d2)


def _sd_S(Fpm, d1, d2, aA, aT):
    """S(q) of HMS Eqs. 80-82 (Hu2022 Eq. 2) from the isospin responses {(response, '+'/'-', L): F}"""
    a0 = aA[0] + aT[0] + aA[1] + aT[1]
    a1 = aA[0] + aT[0] - aA[1] - aT[1]
    a1A = aA[0] - aA[1]
    S = 0.
    for resp, d in (("Sigma'", d1), ("Sigma''", d2)):
        for L in sorted({k[2] for k in Fpm if k[0] == resp}):
            S = S + (a0 * Fpm[(resp, "+", L)] + (a1 + d * a1A) * Fpm[(resp, "-", L)])**2
    return S


def sd_structure_function(A, Z, q, aA=(0., 0.), aT=(0., 0.), two_body="central", sd_model="HMS2020"):
    """S(q) of HMS Eqs. 80-82 generalized to the couplings a_N = aA_N + aT_N (N = p, n); the corrections
    delta', delta'' multiply the isovector part of the axial-vector couplings aA only.
    With aA = (1, 0) (aA = (0, 1)) it is S_p (S_n) of HMS Eqs. 81-82: S_N(0) = (2J+1)(J+1)/(pi J) <S_N>^2.

    sd_model = "HMS2020": two_body is 'none', a key of TWO_BODY_SETTINGS or (rho, delta a(0)).
    sd_model = "Hu2022": two_body = (interaction, rho [fm^-3] or None for no corrections); tensor couplings aT
    are not allowed (the Hu2022 two-body currents are those of the axial-vector current)."""
    q = np.asarray(q, dtype=float)
    if sd_model == "HMS2020":
        F = sd_responses(A, Z, q)
        if F is None:
            return np.zeros_like(q)
        Fpm = {}
        for (r, N, L) in F:
            Fpm[(r, "+", L)] = F[(r, "p", L)] + F[(r, "n", L)]
            Fpm[(r, "-", L)] = F[(r, "p", L)] - F[(r, "n", L)]
        d1, d2 = sd_deltas(q, two_body)
    elif sd_model == "Hu2022":
        if any(x != 0. for x in aT):
            raise ValueError("The Hu2022 structure factors include the two-body axial-vector currents and apply to "
                             "axial-vector (AA) couplings only: use sd_model='HMS2020' for tensor (TT) couplings")
        interaction, rho = two_body
        Fpm = sd_responses_hu(A, Z, q, interaction)
        if Fpm is None:
            return np.zeros_like(q)
        d1, d2 = hu_deltas(q, interaction, rho)
    else:
        raise ValueError("sd_model must be one of %s" % (SD_MODELS,))
    return np.zeros_like(q) + _sd_S(Fpm, d1, d2, aA, aT)


def structure_factor_SD(A, Z, q, c, two_body="central", jchi=0.5, sd_model="HMS2020"):
    """|F_SD(q)|^2 = 4 pi/(2J+1) j_chi(j_chi+1)/12 S(q) [GeV^-4] (4 pi v^2 dsigma/dq^2 of the SD part);
    two_body as in sd_structure_function"""
    J = nuclear_spin(A, Z)
    q = np.asarray(q, dtype=float)
    if J == 0 or not c.has_sd:
        return np.zeros_like(q)
    S = sd_structure_function(A, Z, q, c.c4A, c.c4T, two_body, sd_model)
    return 4. * np.pi / (2. * J + 1.) * jchi * (jchi + 1.) / 12. * S


########################################################################################
# Rates
########################################################################################
_halo_c = _crapidd.halo
_halo_c.argtypes = [ctypes.c_double, ctypes.c_int]
_halo_c.restype = ctypes.c_double


def halo_eta(vmin, i=0):
    """Halo integral of the table loaded with read_halo, eta_0(vmin) = int f(v)/v d^3v [(km/s)^-1] for i = 0
    (vmin in km/s; 0 outside the table)"""
    vmin = np.asarray(vmin, dtype=float)
    return np.vectorize(lambda v: _halo_c(float(v), int(i)), otypes=[float])(vmin)


def isotope_mass_fractions(target):
    """[(A, Z, mass fraction)] of the target, as in lib/source/phys_consts.c (isotope_mass_fractions)"""
    if target not in TARGETS:
        raise ValueError("chiral EFT target must be one of %s" % list(TARGETS))
    Z, As, f = TARGETS[target]
    m = np.array([nucleus_mass(A, Z) for A in As])
    f = np.array(f)
    w = f * m * f.sum() / (f * m).sum()
    return [(A, Z, wi) for A, wi in zip(As, w)]


def _kinematics(rhoDM, mDM, ER, target):
    """[(A, Z, q [GeV], prefactor rho/(2 pi m_chi) eta w_A CONV_RATE)] of the isotopes of target"""
    out = []
    for A, Z, w in isotope_mass_fractions(target):
        mT = nucleus_mass(A, Z)
        mu = mT * mDM / (mT + mDM)
        q = np.sqrt(2. * mT * ER * 1e-6)
        eta = halo_eta(C_KMS * q / (2. * mu), 0)
        out.append((A, Z, q, w * rhoDM / (2. * np.pi * mDM) * eta * CONV_RATE))
    return out


def _sd_members(sd_model, two_body, band):
    """Settings of the SD structure factors and how to combine their rates ('mean', 'min', 'max')"""
    if sd_model not in SD_MODELS:
        raise ValueError("sd_model must be one of %s" % (SD_MODELS,))
    if band not in BANDS:
        raise ValueError("band must be one of %s" % (BANDS,))
    if sd_model == "HMS2020":
        if band == "central":
            return [two_body], "mean"
        return ["central"] + list(TWO_BODY_BAND), band
    if two_body in (None, "none"):
        # one-body Hu2022: the spread of the three interactions
        return [(i, None) for i in HU2022_INTERACTIONS], ("mean" if band == "central" else band)
    if two_body != "central":
        raise ValueError("with sd_model='Hu2022' two_body must be 'central' (two-body currents) or 'none'")
    if band == "central":
        return [(i, HU2022_RHO_CENTRAL) for i in HU2022_INTERACTIONS], "mean"
    return [(i, r) for i in HU2022_INTERACTIONS for r in HU2022_RHO_BAND], band


def difrate_chiral(rhoDM, mDM, ER, couplings, target="Xe", two_body="central", jchi=0.5, parts=False,
                   sd_model="HMS2020", band="central"):
    """
    dR/dE_R [events/kg/day/keV] of the chiral EFT SI + SD structure factors (elastic scattering).

    Parameters
    ----------
    rhoDM : local DM density [GeV/cm^3]
    mDM : DM mass [GeV]
    ER : recoil energies [keV] (scalar or array)
    couplings : ChiralCouplings
    target : 'Xe', 'Ge', 'Ar' or 'F'
    two_body : two-body currents of the axial-vector SD response. HMS2020: 'central' (default), 'none', a key of
        TWO_BODY_SETTINGS or (rho [fm^-3], delta a(0)), used with band = 'central'. Hu2022: 'central' (with the
        two-body currents) or 'none' (one-body structure factors)
    jchi : DM spin (the SD part scales with j_chi (j_chi+1); 0 for scalar DM)
    parts : if True return {'SI': ..., 'SD': ..., 'total': ...}
    sd_model : SD structure factors, 'HMS2020' (default, nominal LZ) or 'Hu2022' (axial-vector couplings only)
    band : 'central', or the 'min' / 'max' envelope of the SD rate over the nuclear uncertainty:
        HMS2020: the c_D / rho settings TWO_BODY_BAND and 'central'; Hu2022: the three chiral interactions and the
        densities HU2022_RHO_BAND (central: mean of the three interactions at rho = 0.10 fm^-3).
        The SI part has no band.

    The halo integral is the one of the table loaded with read_halo (SHM annual average by default).
    """
    ER = np.asarray(ER, dtype=float)
    scalar = ER.ndim == 0
    ER = np.atleast_1d(ER)
    if sd_model == "Hu2022" and couplings.has_tensor:
        raise ValueError("The Hu2022 SD structure factors apply to axial-vector (AA) couplings only: the model has "
                         "tensor (TT) couplings, use sd_model='HMS2020'")
    members, combine = _sd_members(sd_model, two_body, band)
    kin = _kinematics(rhoDM, mDM, ER, target)
    rate_SI = np.zeros_like(ER)
    if couplings.has_si:
        for A, Z, q, pref in kin:
            rate_SI += pref * structure_factor_SI(target, A, q, couplings)
    rate_SD = np.zeros_like(ER)
    if couplings.has_sd:
        rates = []
        for m in members:
            r = np.zeros_like(ER)
            for A, Z, q, pref in kin:
                r += pref * structure_factor_SD(A, Z, q, couplings, m, jchi, sd_model)
            rates.append(r)
        rates = np.array(rates)
        rate_SD = {"mean": np.mean, "min": np.min, "max": np.max}[combine](rates, axis=0)
    total = rate_SI + rate_SD
    if scalar:
        rate_SI, rate_SD, total = rate_SI[0], rate_SD[0], total[0]
    if parts:
        return {"SI": rate_SI, "SD": rate_SD, "total": total}
    return total


def difrate_chiral_band(rhoDM, mDM, ER, couplings, target="Xe", jchi=0.5, sd_model="HMS2020"):
    """Central rate and the envelope (min, max) of the SD nuclear uncertainty (see difrate_chiral, band):
    HMS2020 c_D / rho band or Hu2022 band (interactions and densities). The SI part has no band."""
    return tuple(difrate_chiral(rhoDM, mDM, ER, couplings, target, "central", jchi, sd_model=sd_model, band=b)
                 for b in BANDS)


########################################################################################
# Equivalent DM-nucleon cross sections for the LZ limits
########################################################################################
LZ24_EFF_PATH = os.path.join(base_dir, '..', 'lib', 'efficiency_tables', 'LZ_NR_2024.dat')
# reference cross section of the equivalent cross sections (sigma_eq does not depend on it)
SIGMA_REF = 1e-40   # cm^2
EQUIVALENT_SIGMA_KEYS = ("SI", "SD_n", "SD_p")


def helm_form_factor_sq(q, A):
    """Helm form factor squared (Lewin & Smith, Astropart. Phys. 6, 87 (1996), Eqs. 4.7 and 4.11: c = 1.23 A^1/3 -
    0.60 fm, a = 0.52 fm, s = 0.9 fm), the SI form factor of the LZ limits (arXiv:2410.17036); q in GeV"""
    q = np.asarray(q, dtype=float)
    c, a, s = 1.23 * A**(1. / 3.) - 0.60, 0.52, 0.9
    rn = np.sqrt(c**2 + 7. / 3. * np.pi**2 * a**2 - 5. * s**2) / HBARC      # GeV^-1
    x = np.maximum(q * rn, 1e-8)
    j1 = np.sin(x) / x**2 - np.cos(x) / x
    return np.where(q * rn > 1e-8, (3. * j1 / x)**2, 1.) * np.exp(-(q * s / HBARC)**2)


def difrate_helm_SI(rhoDM, mDM, ER, c, target="Xe"):
    """dR/dE_R [events/kg/day/keV] of the standard isoscalar SI interaction (O_1 coefficient c [GeV^-2] for protons
    and neutrons) with the Helm form factor: the signal model of the LZ SI limit"""
    ER = np.atleast_1d(np.asarray(ER, dtype=float))
    rate = np.zeros_like(ER)
    for A, Z, q, pref in _kinematics(rhoDM, mDM, ER, target):
        rate += pref * (c * A)**2 * helm_form_factor_sq(q, A)
    return rate


@functools.lru_cache(maxsize=None)
def _lz24_efficiency():
    t = np.loadtxt(LZ24_EFF_PATH)
    E = np.linspace(max(t[0, 0], 0.1), t[-1, 0], 400)
    eff = np.interp(E, t[:, 0], t[:, 1], left=0., right=0.)
    sel = eff > 0
    return E[sel], eff[sel]


def lz_counts(spectrum_fn):
    """Integral of spectrum_fn(E [keV]) times the LZ 2024 NR efficiency (lib/efficiency_tables/LZ_NR_2024.dat):
    proportional to the expected LZ signal (exposure and density cancel in the ratios)"""
    E, eff = _lz24_efficiency()
    return np.trapezoid(spectrum_fn(E) * eff, E)


def _reference_couplings(mDM, sigma_ref=SIGMA_REF):
    """O_1 / O_4 coefficients [GeV^-2] of the reference models at sigma_ref [cm^2]: isoscalar SI (nucleon mass
    (m_p + m_n)/2), neutron-only and proton-only SD (sigma_SI = mu_N^2 c^2/pi, sigma_SD = 3 mu_N^2 c^2/(16 pi))"""
    s = sigma_ref / GEV2_TO_CM2
    mu = lambda m: mDM * m / (mDM + m)
    return (np.sqrt(np.pi * s) / mu(mnucleon), np.sqrt(16. * np.pi * s / 3.) / mu(mneutron),
            np.sqrt(16. * np.pi * s / 3.) / mu(mproton))


def equivalent_sigma(couplings, mDM, jchi=0.5, sigma_ref=SIGMA_REF, halo=None, reference_halo=None):
    """
    DM-nucleon cross sections [cm^2] equivalent to the model in LZ, one per published LZ limit:

      sigma_eq = sigma_ref N_model / N_ref(sigma_ref),

    with N the expected LZ signal (spectrum times the LZ 2024 NR efficiency, lz_counts) at the same DM mass and
    density. N_model uses the halo table halo (the table loaded with read_halo if None), N_ref the halo table
    reference_halo (the same as N_model if None): the published curves assume the standard SHM (core.halo_path),
    so with reference_halo = core.halo_path sigma_eq includes the effect of a different halo of the model.
      N_model : the full model (SI + SD, protons and neutrons, all the xenon isotopes), chiral EFT with HKMS (SI)
                and the HMS central structure factors (SD);
      N_ref   : the reference model of each published curve at the cross section sigma_ref:
                'SI'   isoscalar O_1 with the Helm form factor (sigma = mu_N^2 c^2/pi, mu_N with (m_p + m_n)/2),
                'SD_n' neutron-only axial-vector coupling, HMS central with two-body currents (mu_n),
                'SD_p' proton-only axial-vector coupling, HMS central with two-body currents (mu_p),
                with the same j_chi as the model.
    sigma_eq is compared directly with the published curve: the model is excluded if sigma_eq > sigma_lim.

    sigma_eq does not depend on sigma_ref nor on the density. N_model is not split into SI and SD: LZ sees the
    sum, so a model with SI and SD couplings has an SI sigma_eq that includes the SD signal (and vice versa),
    which is correct and conservative with respect to splitting the signal. The proton-neutron interference
    (e.g. the cancellation of the SI signal at f_n/f_p = -Z/(A-Z)) and the two-body currents are included.
    The comparison is exact when the recoil spectrum of the model has the same shape as the reference one (the LZ
    likelihood then depends only on the normalization): this holds to a good approximation for standard SI and SD
    couplings, while for q- or v-dependent operators and light mediators it is an approximation (still better
    than comparing sigma_p and sigma_n separately).

    Returns {'SI': ..., 'SD_n': ..., 'SD_p': ...}: 0 for a model without signal, nan if the reference model has no
    signal in LZ (DM mass below threshold, or j_chi = 0 for SD).
    """
    from .core import read_halo
    rho = 0.3
    if halo is not None:
        read_halo(halo)
    n_model = lz_counts(lambda E: difrate_chiral(rho, mDM, E, couplings, "Xe", jchi=jchi))
    c_si, c_n, c_p = _reference_couplings(mDM, sigma_ref)
    if reference_halo is not None:
        read_halo(reference_halo)
    try:
        n_ref = {"SI": lz_counts(lambda E: difrate_helm_SI(rho, mDM, E, c_si, "Xe")),
                 "SD_n": lz_counts(lambda E: difrate_chiral(rho, mDM, E, ChiralCouplings(c4A=(0., c_n)), "Xe", jchi=jchi)),
                 "SD_p": lz_counts(lambda E: difrate_chiral(rho, mDM, E, ChiralCouplings(c4A=(c_p, 0.)), "Xe", jchi=jchi))}
    finally:
        # leave loaded the halo of the model
        if reference_halo is not None and halo is not None:
            read_halo(halo)
    out = {}
    for k in EQUIVALENT_SIGMA_KEYS:
        if not n_ref[k] > 0.:
            out[k] = np.nan
        elif not n_model > 0.:
            out[k] = 0.
        else:
            out[k] = sigma_ref * n_model / n_ref[k]
    return out
