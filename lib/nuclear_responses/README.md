# Nuclear responses of the DM-nucleus scattering

Fit coefficients of the nuclear responses F_{O1,O2}(E_R) used by `FF` (`lib/source/effFormFact.c`), one
file per nuclear model. They are read by `read_nuclear_responses(dir)` (C), which `rapidd` calls when it
is imported (`rapidd.read_nuclear_responses()`, default directory `lib/nuclear_responses`).

| file | `basis` of `rapidd.difrate_dER` | source | isotopes |
|---|---|---|---|
| `pn_Fitz.txt` | `pn_Fitz` | arXiv:1203.3542 (outdated: the fits diverge at E_R > 100 keV) | 16O, 19F, 23Na, 28Si, 40Ar, 40Ca, 70-76Ge, 127I, 128-136Xe |
| `pn_SN100PN.txt` | `pn_SN100PN` | arXiv:2510.25197 (SN100PN interaction) | 128-136Xe |
| `iso_Fitz.txt` | `iso_Fitz` | DMFormFactor v6 of arXiv:1308.6288v1 | 128-136Xe |
| `iso_GCN5082.txt` | `iso_GCN5082` | DMFormFactor v6 with the one-body density matrices of the GCN5082 shell model ([Elastic](https://github.com/Berkeley-Electroweak-Physics/Elastic)) | 128-136Xe |
| `iso_GCN5082_MEN.txt`, `iso_GCN5082_MEN_WimPy.txt` | none (`FF_iso_GCN5082_MEN`, not among the models of `FF`) | arXiv:1605.08043 (M, Phi''), [WimPy](https://github.com/bradkav/WIMpy_NREFT) (spin-dependent responses of 129Xe, 131Xe) | 40Ar, 128-136Xe |

The coefficients are those that were hardcoded in `effFormFact.c` until 2026-11-23, written with the
shortest decimal representation of the same double-precision numbers (the responses are unchanged
bit by bit).

## Format

Lines starting with `#` are comments (each file starts with a header describing it, and each nucleus
with a `# <A><element> (Z = ...)` line); empty lines are skipped. The columns are separated by spaces.

**Squared responses** (`pn_*.txt`, `iso_*.txt` except `iso_GCN5082_MEN.txt`):

    A  Z  O1  O2  N1  N2  c  a0 ... a10

F_{O1,O2}^{N1,N2}(E_R) = norm * exp(-2y) * sum_{k=0}^{10} a_k y^k, with y = (q b / 2)^2, q^2 = 2 m_T E_R,
b^2 the harmonic-oscillator length squared of the nucleus (`B2` in `effFormFact.c`) and m_T the mass of
the nucleus (A, Z). The basis of the model (prefix of the name) sets the nucleon labels and the column c:

- `pn_*` (proton-neutron basis): N1, N2 = `p`, `n`; c is the prefactor, norm = c (mostly 1; the
  prefactors 4 pi x ... of some light nuclei of `pn_Fitz.txt` are spelled out in the comment of the
  nucleus).
- `iso_*` (isospin basis): N1, N2 = `+`, `-` with +- = (p +- n) / 2; c = j_N is the spin of the
  nucleus, norm = 16 pi / (2 j_N + 1).

**One-nucleon form factors** (`iso_GCN5082_MEN.txt`, arXiv:1605.08043, not squared):

    A  Z  O  tau  a0 ... a5

F_O^tau(E_R) = exp(-y/2) * sum_{k=0}^{5} a_k y^k, with y = (q b)^2 / 2 and
b = 1 / sqrt(m_N (45 A^(-1/3) - 25 A^(-2/3)) MeV), m_N = 938.919 MeV. The response is the product
F_{O1,O2}^{tau1,tau2} = F_O1^tau1 F_O2^tau2; the responses without both factors are taken from
`iso_GCN5082_MEN_WimPy.txt` (squared responses, `iso` format).

**Labels of the responses**: M = M, Sp = Sigma', Spp = Sigma'', D = Delta, Ppp = Phi'', Pp = Phi', Pt
(Pp and Pt enter O_12 and O_15 as Pp + Pt, see `FormFact_v0`); O1 != O2 are the interferences
(M-Ppp, Sp-D).

**Lookup**: `FF(model, A, Z, ...)` returns the response of the row with the same A, Z, O1, O2, N1, N2,
0 if there is none (e.g. 40Ar and 40Ca of `pn_Fitz.txt` are distinct nuclei).

`read_nuclear_responses` checks the number of columns, the nucleon labels, 1 <= A <= 199 (the size
of `B2`) and that no row (A, Z, O1, O2, N1, N2) is repeated, and returns 2 for a missing file and 3
for a malformed or repeated line (with its line number on stderr); on an error the previous tables are
kept.
