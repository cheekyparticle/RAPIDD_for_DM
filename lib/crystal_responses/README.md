# Silicon crystal response W_1 and ionization probabilities

Data of the DM-electron rate in silicon (`target="Si"` in `rapidd.electron`, `difrate_electron_crystal` in
`lib/source/difRateGen.c`, tables read by `lib/source/crystal_resp.c`).

## Files

| File | Content |
|---|---|
| `Si_W1.txt` | W_1(q, Delta E), unscreened: Delta E q^2 V_cell Im eps(q, Delta E) / (pi^2 alpha) |
| `Si_W1_screened.txt` | W_1(q, Delta E), RPA screened: Delta E q^2 V_cell Im[-1/eps(q, Delta E)] / (pi^2 alpha) |
| `Si_W1_nolfe.txt` | W_1(q, Delta E), independent-particle response: Delta E q^2 V_cell Im eps_noLFE(q, Delta E) / (pi^2 alpha) |
| `Si_ionization_p100K.dat` | probability p_Q(Delta E) of Q electron-hole pairs at 100 K |
| `generate/qcdark2_to_rapidd.py` | conversion of the QCDark2 dielectric functions into the three W_1 tables |

W_1 is the first crystal response of arXiv:2105.02233 (Eq. R_crystal_2D), dimensionless (natural units), with
|f_crystal|^2 = q^3 W_1/(8 Delta E alpha m_e^2) (arXiv:2105.02233, App. B) and
|f_crystal|^2 = q^5 V_cell Im eps/(8 pi^2 (alpha m_e)^2) (arXiv:2104.12786, Eq. 16). eps is the composite dielectric
function (local field effects, LFE, for q <= 8 alpha m_e), eps_noLFE the one without LFE (diagonal G = G' elements),
whose imaginary part is the independent-particle response. In the RAPIDD rate, following arXiv:2402.06817
(non-spin-polarised material), the screened table is used for the couplings to the electron density (c_1, c_11) and
`Si_W1_nolfe.txt` for the couplings to the electron spin (c_4, c_6, c_9, c_10), which are not screened and whose
response is the single-particle one (arXiv:2402.06817 neglects local-field corrections); with `screening=False` the
density couplings use `Si_W1.txt`.

## Format of the W_1 tables

Header lines start with `//` (source, commit, V_cell, M_cell of QCDark2, grid). Then 501 rows, Delta E = 0.1 i eV
(i = 0 ... 500, 0 ... 50 eV), of 1251 columns, q = 0.01 + 0.02 j alpha m_e (j = 0 ... 1250, up to 25.01 alpha m_e =
93.5 keV), tab separated. Im eps and Im[-1/eps] below 1e-10 times the maximum of their q row (numerical noise below
the gap, |Im eps| <= 1e-16) and negative values are set to zero.

## Sources

- Dielectric functions: QCDark2, Dreyer, Essig, Fernandez-Serra, Hott, Singal, arXiv:2603.12326,
  https://github.com/meganhott/QCDark2 (CC0 1.0), commit 6d22f936bf49a15db73a24d1274dcc63c37780dd,
  `dielectric_functions/nolfe/Si_nolfe.h5` (`Si_cc-pv(t+d)z_pbe_8k_25q_nolfe_qs0.066286`, for `Si_W1_nolfe.txt`) and
  `dielectric_functions/composite/Si_comp.h5` (`Si_cc-pv(t+d)z_pbe_8k_8qlfe_25qnolfe_qs0.066286`): all-electron
  DFT (PBE, cc-pv(t+d)z, scissor-corrected gap), RPA with local field effects for q <= 8 alpha m_e and without them
  up to 25 alpha m_e, transitions up to 50 eV (no Si 2p). V_cell = 270.107 bohr^3. The unit cell mass of the
  RAPIDD rate is 2 x 28.0855 u (QCDark2 uses M_cell = 52.164 GeV, 0.3% lower: the rates per kg differ by that factor).
- Ionization probabilities: Ramanathan and Kurinsky, Phys. Rev. D 102, 063026 (2020), arXiv:2004.10709,
  supplementary file `p100K.dat` of the arXiv source (`Supplementary/p100K.dat`, distributed with the arXiv
  non-exclusive license; unmodified, sha256 8ed7593e63daaf7758a3a6ff9282857585c1fdeb331ee7d56753e7668e6608f1).
  Tab separated: Delta E [eV] from 1.1 to 50 eV in steps of 0.05 eV, then p_Q for Q = 1 ... 20 (rows normalized
  to 1). Used by QCDark2 (`recoil_spectrum`) and by DAMIC-M (arXiv:2503.14617).

## Regenerating the W_1 tables

Requires h5py (the QCDark2 Python code is not used):

    python generate/qcdark2_to_rapidd.py <QCDark2>/dielectric_functions/composite/Si_comp.h5 \
        <QCDark2>/dielectric_functions/nolfe/Si_nolfe.h5 \
        --commit <QCDark2 commit> --out-dir . --reference ../../tests/data/si_qcdark2_reference.txt

`--reference` also writes the dR/dE of the QCDark2 rate formula used by `tests/test_electron_silicon.py`.
