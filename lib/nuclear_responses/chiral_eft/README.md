# Chiral EFT structure factors (rapidd.chiral_eft)

Fit coefficients of the chiral EFT responses used by `rapidd/chiral_eft.py` (not read by `read_nuclear_responses`,
which reads only the `*.txt` files of the parent directory). The tables follow the conventions of
`lib/nuclear_responses`: a header describing the source, the fit function and the columns, one section per
nucleus (`# <A><element> (Z = ...)`), aligned columns separated by spaces, comments starting with `#`.
Only the nuclei of the RAPIDD targets (F, Ar, Ge, Xe) are included. Equation numbers refer to the arXiv versions.

| file | content | source | nuclei |
|---|---|---|---|
| `HKMS2019_SI.txt` | SI responses M, Phi'' (isospin +, -), F_pi, F_b | arXiv:1812.05617 (HKMS), ChiralEFT4DM notebook v1.0 | 19F, 40Ar, 70-76Ge, 128-136Xe |
| `HMS2020_SD.txt` | one-body SD responses Sigma', Sigma'' (p, n) | arXiv:2007.08529 (HMS), App. E | 19F, 73Ge, 129Xe, 131Xe |
| `Hu2022_SD.txt` | one-body SD responses Sigma', Sigma'' (+, -) of three chiral interactions | arXiv:2109.00193, AbInitio4DarkMatter notebook | 19F, 73Ge, 129Xe, 131Xe |
| `spin_expectation.txt` | <S_p>, <S_n> of the SD tables | HMS Table 4, Hu2022 notebook | 19F, 73Ge, 129Xe, 131Xe |

All the fits have the form F(u) = exp(-u/2) sum_{k=0}^{5} a_k u^k, u = (q b)^2/2 (note: u, not the
y = (q b/2)^2 of the NREFT tables).

## Formats

    HKMS2019_SI.txt       A  Z  O  tau  L  a0 ... a5
    HMS2020_SD.txt        A  Z  O  N    L  J  b  a0 ... a5
    Hu2022_SD.txt         A  Z  O  tau  L  interaction  J  b  a0 ... a5
    spin_expectation.txt  A  Z  model  interaction  S_p  S_n

- O: `M` = F^M, `Ppp` = F^Phi'', `Pi` = F_pi, `B` = F_b (SI); `Sp` = Sigma', `Spp` = Sigma'' (SD).
- tau: isospin `+`, `-` with F_+- = F_p +- F_n (**not** divided by 2, unlike the iso_* NREFT tables); `0` for the
  two-body SI responses F_pi, F_b. N: `p`, `n`.
- L: multipole. SI: 0 for the responses of HKMS Eq. 1, 2 for the L = 2 multipole of 131Xe (added incoherently).
- J: spin of the nucleus; b: harmonic-oscillator length [fm]. For the SI table b is not tabulated:
  b = 1/sqrt(m_N hbar omega), hbar omega = (45 A^-1/3 - 25 A^-2/3) MeV, m_N = 938.919 MeV (HKMS Eqs. C14-C15).

## SI: `HKMS2019_SI.txt`

M. Hoferichter, P. Klos, J. Menéndez, A. Schwenk, PRD 99, 055031 (2019), arXiv:1812.05617 (HKMS), with the xenon
responses of PRD 94, 063505 (2016), arXiv:1605.08043.

The coefficients are copied with the same digits from the `<el>_wimps_<response>_fit.dat` files of the ChiralEFT4DM
notebook v1.0 (https://theorie.ikp.physik.tu-darmstadt.de/strongint/ChiralEFT4DM.html, `ChiralEFT4DM.tar.gz`,
downloaded on 2026-11-17), for `<el>` = f, ar, ge, xe (the Si files are not used), with the responses
`m+`, `m-` -> `M +`, `M -`; `phi+`, `phi-` -> `Ppp +`, `Ppp -`; `pi` -> `Pi 0`; `eb` -> `B 0` and the column `131L2`
-> `L = 2`. F^M_+(0) = A and F^M_-(0) = Z - N. F_b is the two-body structure factor of HKMS Eq. 24, normalized to
F_b(0) = -2 E_b/M_pi with the binding energy E_b; its coupling is c_b = f_pi^theta (HKMS Eq. 7), and the
trace-anomaly response is F_pi^theta = 2 F_pi + F_b (HKMS Eq. 25).

Cross-check with the tables of arXiv:1605.08043: the xenon M and Phi'' coefficients (Table 4) and F_pi (Table 7)
coincide, except F_pi of 132Xe, whose c_1 ... c_5 differ by 0.3-3% (2.11819 vs 2.11087 for c_1); the notebook values
are used. F_b is the response of HKMS, not the F_pi^theta of Table 7 of arXiv:1605.08043 (which did not include the
contact terms).

## SD, nominal: `HMS2020_SD.txt`

M. Hoferichter, J. Menéndez, A. Schwenk, PRD 102, 074018 (2020), arXiv:2007.08529 (HMS), App. E: one-body
responses F_N^{Sigma'_L}, F_N^{Sigma''_L} (N = p, n; L = 1, 3, ...) of 19F, 129Xe, 131Xe (Table 8) and 73Ge (Table 11);
b from Tables 6-7 (19F, 73Ge) and Table 4 of arXiv:1605.08043 (Xe; all equal to the formula above within 1e-4 fm),
J and <S_p>, <S_n> (in `spin_expectation.txt`) from Table 4 of HMS. The coefficients were extracted from the LaTeX
source of arXiv:2007.08529 (v2). These are the structure factors of the nominal LZ 2024 SD limits.

Normalization (HMS Eq. 69): F_N^{Sigma'_1}(0) = sqrt(2) F_N^{Sigma''_1}(0) = sqrt(2/3) sqrt((2J+1)(J+1)/(4 pi J)) <S_N>,
i.e. S_N(0) = (2J+1)(J+1)/(pi J) <S_N>^2 (checked in tests/test_chiral_eft.py).

The corrections delta'(q^2), delta''(q^2) of the axial-vector response (HMS Eqs. 77-78, 86) are computed in
rapidd/chiral_eft.py: axial radius <r_A^2> = 0.46 fm^2, induced pseudoscalar (g_piNN^2/(4 pi) = 13.7, F_pi = 92.28 MeV,
g_A = 1.27641, M_pi = 138 MeV, m_N = 938.92 MeV) and two-body currents with c_1, c_3, c_4, c_6 of HMS Table 5 and the
Fermi-gas integrals of Klos et al. (arXiv:1304.7684, App. B) evaluated by quadrature. c_D is fixed by the quenching
delta a(0) = -(20 ... 30)% at rho = 0.09 ... 0.11 fm^-3 (central: 25% at 0.10 fm^-3); the extremes reproduce the
range c_D = -6.08 ... 0.30 of HMS Table 5.

## SD, uncertainty band: `Hu2022_SD.txt`

B. S. Hu, J. Padua-Argüelles, S. Leutheusser, T. Miyagi, S. R. Stroberg, J. D. Holt, PRL 128, 072502 (2022),
arXiv:2109.00193 (Hu2022): ab initio (VS-IMSRG) one-body responses F_+-^{Sigma'_L}, F_+-^{Sigma''_L}
(F_+- = F_p +- F_n) of 19F, 73Ge, 129Xe, 131Xe for the chiral interactions 1.8/2.0(EM) (`EM1.8_2.0`), N4LO+3N_lnl
(`N4LO_LNL2`) and Delta NNLO_GO(394) (`DNNLOgo`). The LZ 2024 SD limits use their spread as the nuclear uncertainty.

- Source: the public notebook of the paper, https://github.com/BaishanHu/AbInitio4DarkMatter
  (`ChEFT_SpinDependent`, commit e10cccbbdb4b792a52b1352519ef9454eb1e403d of 2021-10-02, downloaded on 2026-11-24).
  The coefficients are copied with the same digits from the files of `ChEFT_SpinDependent/data_fit` that the
  notebook `DarkMatter_ChEFT_SD.ipynb` reads (its dictionary `Dic_IMSRG`; the file names are in the comment of each
  nucleus); the other nuclei of the notebook (23Na, 27Al, 29Si, 127I) and the convergence studies are not copied.
- b = hbar c/sqrt(m_N hbar omega) with the hbar omega = 16 MeV of the IMSRG basis (as in the notebook; not the
  A-dependent formula of HKMS and HMS) and m_N = (m_p + m_n)/2: b = 1.609950 fm for all the nuclei.
- The constants of `DM_constants.py` (g_A = 1.29, g_piNN = 13.1066, F_pi = 92.4 MeV, M_pi = 138.04 MeV,
  Lambda_chi = 700 MeV, <r_A^2> = 0.47 fm^2) and the low-energy constants of cell 6 of the notebook (c_1, c_3, c_4, c_D
  of each interaction, hat c_6 = 3.704) are `HU2022_CONSTANTS` and `HU2022_LECS` of rapidd/chiral_eft.py.
  `hu_delta_a` implements `DM_func.delta_2BC` (HMS Eqs. 77-78 with c_4 -> c_4 + 1/(4 m_N) and
  c_6 -> (1 + hat c_6)/m_N) and the corrections of Hu2022 Eq. 3; the delta', delta'' of the notebook are reproduced
  to 1e-6 (tests/test_chiral_eft.py).
- <S_p>, <S_n> in `spin_expectation.txt`: the values of each interaction (table of the notebook, Fig. 2 of the
  paper). The one-body S_N(0) of the fits reproduces (2J+1)(J+1)/(pi J) <S_N>^2 within 0.5% (same normalization as
  HMS).
- The structure factors with the two-body axial currents are those of the axial-vector (AA) coupling only. The paper
  and the notebook give no central value and no band files: rapidd.chiral_eft builds the band as the envelope of the
  rates of the three interactions for rho = 0.08, 0.09, 0.10, 0.11, 0.12 fm^-3 (the range quoted in the paper; the
  notebook default is 0.09 ... 0.12 fm^-3), and the central value as the mean of the three interactions at
  rho = 0.10 fm^-3.
