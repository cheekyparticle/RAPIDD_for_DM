# Chiral EFT structure factors (rapidd.chiral_eft)

Fit coefficients of the chiral EFT responses used by `rapidd/chiral_eft.py`, for the nuclei of the RAPIDD targets.

| file | content | source | nuclei |
|---|---|---|---|
| `HKMS2019_SI.txt` | SI responses M, Phi'' (isospin +, -), F_pi, F_b | arXiv:1812.05617, ChiralEFT4DM notebook v1.0 | 19F, 40Ar, 70-76Ge, 128-136Xe |
| `HMS2020_SD.txt` | one-body SD responses Sigma', Sigma'' (p, n) | arXiv:2007.08529, App. E | 19F, 73Ge, 129Xe, 131Xe |
| `Hu2022_SD.txt` | one-body SD responses Sigma', Sigma'' (+, -) of three chiral interactions | arXiv:2109.00193, AbInitio4DarkMatter notebook | 19F, 73Ge, 129Xe, 131Xe |
| `spin_expectation.txt` | <S_p>, <S_n> of the SD tables | arXiv:2007.08529 Table 4, AbInitio4DarkMatter notebook | 19F, 73Ge, 129Xe, 131Xe |

All the fits have the form F(u) = exp(-u/2) sum_{k=0}^{5} a_k u^k with u = (q b)^2/2 (not the y = (q b/2)^2 of the
NREFT tables).

## Formats

    HKMS2019_SI.txt       A  Z  O  tau  L  a0 ... a5
    HMS2020_SD.txt        A  Z  O  N    L  J  b  a0 ... a5
    Hu2022_SD.txt         A  Z  O  tau  L  interaction  J  b  a0 ... a5
    spin_expectation.txt  A  Z  model  interaction  S_p  S_n

- O: `M` = F^M, `Ppp` = F^Phi'', `Pi` = F_pi, `B` = F_b (SI); `Sp` = Sigma', `Spp` = Sigma'' (SD).
- tau: `+`, `-` with F_+- = F_p +- F_n (not divided by 2, unlike the NREFT tables); `0` for F_pi and F_b. N: `p`, `n`.
- L: multipole (the L = 2 term of 131Xe in the SI table is added incoherently).
- J: spin of the nucleus; b: harmonic-oscillator length [fm]. In the SI table b = 1/sqrt(m_N hbar omega) with
  hbar omega = (45 A^-1/3 - 25 A^-2/3) MeV and m_N = 938.919 MeV (arXiv:1812.05617, Eqs. C14-C15).

## Sources

- `HKMS2019_SI.txt`: Hoferichter, Klos, Menéndez, Schwenk, PRD 99, 055031 (2019). Coefficients copied from the
  `<el>_wimps_<response>_fit.dat` files of the ChiralEFT4DM notebook v1.0
  (https://theorie.ikp.physik.tu-darmstadt.de/strongint/ChiralEFT4DM.html), with `m+`, `m-` -> `M +`, `M -`;
  `phi+`, `phi-` -> `Ppp +`, `Ppp -`; `pi` -> `Pi 0`; `eb` -> `B 0`. F_b is normalized to F_b(0) = -2 E_b/M_pi
  (Eq. 24 of the paper).
- `HMS2020_SD.txt`: Hoferichter, Menéndez, Schwenk, PRD 102, 074018 (2020), Tables 8 and 11 (LaTeX source of v2),
  with b from Tables 6-7 and from Table 4 of arXiv:1605.08043 for Xe. Normalization (Eq. 69):
  S_N(0) = (2J+1)(J+1)/(pi J) <S_N>^2. These are the structure factors of the nominal LZ 2024 SD limits; the
  axial-vector corrections delta', delta'' (Eqs. 77-78, 86) are computed in `rapidd/chiral_eft.py`.
- `Hu2022_SD.txt`: Hu et al., PRL 128, 072502 (2022), interactions 1.8/2.0(EM) (`EM1.8_2.0`), N4LO+3N_lnl
  (`N4LO_LNL2`) and Delta NNLO_GO(394) (`DNNLOgo`). Coefficients copied from `ChEFT_SpinDependent/data_fit` of
  https://github.com/BaishanHu/AbInitio4DarkMatter (commit e10cccb), with b = 1.609950 fm (hbar omega = 16 MeV, as
  in the notebook). The LZ 2024 SD limits use the spread of the three interactions as the nuclear uncertainty:
  `rapidd.chiral_eft` takes the envelope of their rates for rho = 0.08 ... 0.12 fm^-3 as the band and their mean at
  0.10 fm^-3 as the central value.
