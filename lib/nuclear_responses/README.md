# Nuclear responses of the DM-nucleus scattering

Coefficients of the nuclear responses F_{O1,O2}(E_R) used by `FF` (`lib/source/effFormFact.c`), one file per
nuclear model, read by `read_nuclear_responses(dir)` when `rapidd` is imported.

| file (`basis` of `rapidd.difrate_dER`) | basis | source | isotopes |
|---|---|---|---|
| `McElvain2022.txt` (default) | +, - | DMFormFactor v6 (arXiv:1308.6288) with the one-body density matrices of the [Elastic](https://github.com/Berkeley-Electroweak-Physics/Elastic) repository (arXiv:2208.07945); interactions in the header | 12C, 13C, 16O, 18O, 19F, 23Na, 27Al, 28-30Si, 40Ar, 40Ca, 70-76Ge, 127I, 133Cs, 128-136Xe |
| `Fitz2012.txt` | p, n | fits of arXiv:1203.3542 (they diverge at E_R > 100 keV) | 16O, 19F, 23Na, 28Si, 40Ar, 40Ca, 70-76Ge, 127I, 128-136Xe |
| `Anand2013.txt` | +, - | DMFormFactor v6 (arXiv:1308.6288) with its default density matrices | 19F, 23Na, 28-30Si, 70-76Ge, 127I, 128-136Xe |
| `Khaleq2023.txt` | p, n | arXiv:2311.15764, App. E.7 (SN100PN interaction) | 128-136Xe |

`Anand2013.txt` and `McElvain2022.txt` are written by the scripts of `generate/`. The chiral EFT structure factors of
`rapidd.chiral_eft` are in `chiral_eft/`.

`rapidd.check_nuclear_model(target, basis)` rejects a model that does not cover all the isotopes of the target
(targets in `nuclear_targets`, `lib/source/phys_consts.c`):

| target | models |
|---|---|
| C, O, Al, Cs | McElvain2022 |
| F, Na, Ge, I | McElvain2022, Fitz2012, Anand2013 |
| Si | McElvain2022, Anand2013 |
| Ar, Ca, CaWO4 | McElvain2022, Fitz2012 |
| Xe | McElvain2022, Fitz2012, Anand2013, Khaleq2023 |

## Format

Lines starting with `#` are comments; each file starts with a header and each nucleus with a
`# <A><element> (Z = ...)` line. Each data line is

    A  Z  O1  O2  N1  N2  c  a0 ... a10

F_{O1,O2}^{N1,N2}(E_R) = norm * exp(-2y) * sum_k a_k y^k, with y = (q b / 2)^2, q^2 = 2 m_T E_R and b the
harmonic-oscillator length of the nucleus (`B2` in `effFormFact.c`). The basis of the model (table `nuclear_tables`
of `effFormFact.c`) sets the nucleon labels and the column c:

- p, n basis (Fitz2012, Khaleq2023): N1, N2 = `p`, `n`; norm = c, the prefactor.
- isospin basis (Anand2013, McElvain2022): N1, N2 = `+`, `-` with +- = (p +- n)/2; c = j_N is the spin of the
  nucleus and norm = 16 pi / (2 j_N + 1).

Responses: M, Sp = Sigma', Spp = Sigma'', D = Delta, Ppp = Phi'', Pp = Phi', Pt = Phi-tilde' (Pp and Pt enter
O_12 and O_13 as Pp + Pt); O1 != O2 are the interferences M-Ppp and Sp-D. A row that is missing gives 0.
