# RAPIDD

RAPIDD is a tool for computing dark matter direct detection rates and exclusion limits. It combines a C library for the numerical computation of differential recoil rates with a Python interface for setting up dark matter models, generating rates, and deriving experimental limits.

The recent RAPIDD version works both in standalone mode and as a plugin of [MadDM](https://github.com/maddmhep/maddm), allowing direct detection rates to be computed directly from a UFO model.

## Authors

### Current maintainers and authors of the code
- A. Cheek — andrew.cheek1992@gmail.com
- Gian Marco Lucchetti -- gianmarco.lucchetti@unibo.it



### Past contributors
- Ellen Sandford E.M.Sandford@liverpool.ac.uk
- D. G. Cerdeño — davidg.cerdeno@gmail.com
- E. Reid
- H. Schulz 
- M. Peiro
- E. Gerstmayr

## Citation

If you use RAPIDD, please cite

- D. G. Cerdeño, A. Cheek, E. Reid, H. Schulz, *Surrogate Models for Direct Dark Matter Detection*, JCAP 08 (2018) 011, [arXiv:1802.03174](https://arxiv.org/abs/1802.03174).

RAPIDD uses formalisms, codes and data of other authors: please cite also the works below for the parts you use (the data files and the code give the details).

**DM-nucleus scattering**
- Non-relativistic EFT and nuclear responses: A. L. Fitzpatrick, W. Haxton, E. Katz, N. Lubbers, Y. Xu, JCAP 02 (2013) 004, [arXiv:1203.3542](https://arxiv.org/abs/1203.3542) (also the `Fitz2012` responses); N. Anand, A. L. Fitzpatrick, W. C. Haxton, Phys. Rev. C 89, 065501 (2014), [arXiv:1308.6288](https://arxiv.org/abs/1308.6288) (DMFormFactor, `Anand2013` and `McElvain2022` responses).
- `McElvain2022` density matrices: W. C. Haxton, E. Rule, K. McElvain, M. J. Ramsey-Musolf, Phys. Rev. C 107, 035504 (2023), [arXiv:2208.07945](https://arxiv.org/abs/2208.07945) ([Elastic](https://github.com/Berkeley-Electroweak-Physics/Elastic) repository).
- `Khaleq2023` responses: R. Abdel Khaleq, G. Busoni, C. Simenel, A. E. Stuchbery, Phys. Rev. D 109, 075036 (2024), [arXiv:2311.15764](https://arxiv.org/abs/2311.15764).
- Inelastic scattering: G. Barello, S. Chang, C. A. Newby, Phys. Rev. D 90, 094027 (2014), [arXiv:1409.0536](https://arxiv.org/abs/1409.0536).
- Chiral EFT: M. Hoferichter, P. Klos, J. Menéndez, A. Schwenk, Phys. Rev. D 99, 055031 (2019), [arXiv:1812.05617](https://arxiv.org/abs/1812.05617) and Phys. Rev. D 94, 063505 (2016), [arXiv:1605.08043](https://arxiv.org/abs/1605.08043); M. Hoferichter, J. Menéndez, A. Schwenk, Phys. Rev. D 102, 074018 (2020), [arXiv:2007.08529](https://arxiv.org/abs/2007.08529); B. S. Hu, J. Padua-Argüelles, S. Leutheusser, T. Miyagi, S. R. Stroberg, J. D. Holt, Phys. Rev. Lett. 128, 072502 (2022), [arXiv:2109.00193](https://arxiv.org/abs/2109.00193); P. Klos, J. Menéndez, D. Gazit, A. Schwenk, Phys. Rev. D 88, 083516 (2013), [arXiv:1304.7684](https://arxiv.org/abs/1304.7684); F. Bishara, J. Brod, B. Grinstein, J. Zupan, JHEP 11 (2017) 059, [arXiv:1707.06998](https://arxiv.org/abs/1707.06998).

**DM-electron scattering**
- Non-relativistic EFT and atomic responses: R. Catena, T. Emken, N. A. Spaldin, W. Tarantino, Phys. Rev. Research 2, 033195 (2020), [arXiv:1912.08204](https://arxiv.org/abs/1912.08204); the atomic response tables are computed with DarkART, T. Emken, [doi:10.5281/zenodo.6046225](https://doi.org/10.5281/zenodo.6046225); sign of W_2: J.-H. Liang, Y. Liao, X.-D. Ma, H.-L. Wang, Phys. Rev. D 110, L091701 (2024), [arXiv:2405.04855](https://arxiv.org/abs/2405.04855).
- Silicon: crystal responses, R. Catena, T. Emken, M. Matas, N. A. Spaldin, E. Urdshals, Phys. Rev. Research 3, 033149 (2021), [arXiv:2105.02233](https://arxiv.org/abs/2105.02233); screening of general interactions, R. Catena, N. A. Spaldin, Phys. Rev. Research 6, 033230 (2024), [arXiv:2402.06817](https://arxiv.org/abs/2402.06817); dielectric functions, QCDark2, C. Dreyer, R. Essig, M. Fernandez-Serra, M. Hott, A. Singal, Phys. Rev. D 114, 035003 (2026), [arXiv:2603.12326](https://arxiv.org/abs/2603.12326), [github.com/meganhott/QCDark2](https://github.com/meganhott/QCDark2); normalization, DarkELF, S. Knapen, J. Kozaczuk, T. Lin, Phys. Rev. D 105, 015014 (2022), [arXiv:2104.12786](https://arxiv.org/abs/2104.12786); ionization yield, K. Ramanathan, N. Kurinsky, Phys. Rev. D 102, 063026 (2020), [arXiv:2004.10709](https://arxiv.org/abs/2004.10709).
- XENON10 analysis and DM-electron rate conventions: R. Essig, T. Volansky, T.-T. Yu, Phys. Rev. D 96, 043017 (2017), [arXiv:1703.00910](https://arxiv.org/abs/1703.00910); R. Essig, M. Fernandez-Serra, J. Mardon, A. Soto, T. Volansky, T.-T. Yu, JHEP 05 (2016) 046, [arXiv:1509.01598](https://arxiv.org/abs/1509.01598).

**Halo**
- D. Baxter et al., Eur. Phys. J. C 81, 907 (2021), [arXiv:2105.00599](https://arxiv.org/abs/2105.00599) (standard halo conventions); C. Savage, K. Freese, P. Gondolo, Phys. Rev. D 74, 043531 (2006), [astro-ph/0607121](https://arxiv.org/abs/astro-ph/0607121); C. McCabe, Phys. Rev. D 82, 023530 (2010), [arXiv:1005.0579](https://arxiv.org/abs/1005.0579); SHM + LMC: N. L. Rodd, B. R. Safdi, T. R. Slatyer, W. L. Xu, [arXiv:2609.04175](https://arxiv.org/abs/2609.04175).

**Experimental results** (`lib/experiments`)
- XENON1T: E. Aprile et al. (XENON), Phys. Rev. Lett. 123, 251801 (2019), [arXiv:1907.11485](https://arxiv.org/abs/1907.11485) ([data release](https://github.com/XENON1T/s2only_data_release)).
- XENONnT: E. Aprile et al. (XENON), Phys. Rev. Lett. 137, 051003 (2026), [arXiv:2601.11296](https://arxiv.org/abs/2601.11296).
- DAMIC-M: K. Aggarwal et al. (DAMIC-M), Phys. Rev. Lett. 135, 071002 (2025), [arXiv:2503.14617](https://arxiv.org/abs/2503.14617) ([limits](https://github.com/DAMIC-M/LBC_2025_HSBenchmark_Limits)); the rescaled limit also uses QCDark, C. E. Dreyer, R. Essig, M. Fernandez-Serra, A. Singal, C. Zhen, Phys. Rev. D 109, 115008 (2024), [arXiv:2306.14944](https://arxiv.org/abs/2306.14944), and the screening model of G. Cappellini, R. Del Sole, L. Reining, F. Bechstedt, Phys. Rev. B 47, 9892 (1993).
- DarkSide-50: P. Agnes et al. (DarkSide-50), Phys. Rev. Lett. 130, 101002 (2023), [arXiv:2207.11968](https://arxiv.org/abs/2207.11968); Phys. Rev. Lett. 121, 111303 (2018), [arXiv:1802.06998](https://arxiv.org/abs/1802.06998).
- LZ: J. Aalbers et al. (LZ), Phys. Rev. Lett. 135, 011802 (2025), [arXiv:2410.17036](https://arxiv.org/abs/2410.17036).
- PICO-60: C. Amole et al. (PICO), Phys. Rev. D 100, 022001 (2019), [arXiv:1902.04031](https://arxiv.org/abs/1902.04031).
- XLZD: J. Aalbers et al. (XLZD), Eur. Phys. J. C 85, 1192 (2025), [arXiv:2410.17137](https://arxiv.org/abs/2410.17137) (curves from [DirectDetectionPlots](https://github.com/cajohare/DirectDetectionPlots), C. A. J. O'Hare); J. Aalbers et al., J. Phys. G 50, 013001 (2023), [arXiv:2203.02309](https://arxiv.org/abs/2203.02309).

## Requirements

Before installing, make sure the following are available on your system:

| Requirement   | Purpose                                   |
|---------------|--------------------------------------------|
| CMake (≥ 3.10)| Builds the RAPIDD C library                |
| GSL           | Numerical routines used by the C library   |
| pkg-config    | Used by CMake to locate GSL                |
| Python (≥ 3.8)| Runs the Python interface                  |
| pip           | Installs the Python package                |

On Debian/Ubuntu these can be installed with:

```bash
sudo apt-get install cmake libgsl-dev pkg-config
```

On macOS, using Homebrew:

```bash
brew install cmake gsl pkg-config
```

The Python dependencies (`numpy`, `scipy`) are installed automatically by `pip` in the step below.

## Installation

Clone the repository:

```bash
git clone https://github.com/cheekyparticle/RAPIDD_for_DM.git
cd RAPIDD_for_DM
```

### 1. Build the C library

From the repository root, run the provided build script:

```bash
sh instructions_cmake.sh
```

This creates `lib/build/` and compiles the shared library (`libRAPIDD.so` on Linux, `libRAPIDD.dylib` on macOS) together with two standalone executables, `rpd` and `rpd3`.

### 2. Install the Python package

From the repository root:

```bash
pip install -e .
```

The `-e` flag installs RAPIDD in editable mode, so changes to the Python source are picked up without reinstalling.

## Verifying the installation

Open a Python interpreter and import the package:

```python
>>> import rapidd
```

If this completes without errors, RAPIDD is installed correctly and the compiled C library was found and loaded.

## Basic usage

The Python interface loads the compiled C library through `ctypes` and exposes it as a set of higher-level modules:

| Module               | Contents                                                        |
|-----------------------|-------------------------------------------------------------------|
| `rapidd.core`         | Low-level bindings to the C library (halo tables, form factors, coefficients) |
| `rapidd.halo`         | Analytic dark matter velocity distributions (e.g. Standard Halo Model) |
| `rapidd.experiments`  | Detector efficiency curves and experiment-specific response functions |
| `rapidd.alpha_map`    | Mapping of effective operator coefficients onto nucleon-level couplings |
| `rapidd.calc_dRdE`    | Computation and export of differential recoil rates              |
| `rapidd.chiral_eft`   | DM-nucleus SI and SD rates with the chiral EFT structure factors (two-body currents) and the equivalent cross sections of the LZ limits |
| `rapidd.electron`     | DM-electron scattering: NREFT ionization rates, S2 signals, models and limits |
| `rapidd.stats`        | Statistical routines for deriving exclusion limits                |

Worked examples are in `examples/evaluate_NR_rates.ipynb` (nuclear recoils) and `examples/evaluate_ER_rates.ipynb` (DM-electron scattering), with the RAPIDD defaults. The validation (comparisons with the literature, with the signal models of the LZ limits, with wimprates and with the XENON1T/XENONnT data releases) is in `examples/checks_NR_rates.ipynb` and `examples/checks_ER_rates.ipynb`.

### DM-nucleus scattering: targets and nuclear models

`rapidd.difrate_dER(rhoDM, mDM, ER, target=..., basis=...)` gives the DM-nucleus rate per kg of target for the NREFT operators:

- targets (`rapidd.NUCLEAR_TARGET_ISOTOPES`): `Xe`, `Ge`, `Ar`, `F`, `C`, `O`, `Na`, `Al`, `Si`, `Ca`, `I`, `Cs` and `CaWO4` (without W);
- nuclear models (`basis`): `McElvain2022` (default), `Fitz2012`, `Anand2013` and `Khaleq2023`. A model that does not cover all the isotopes of the target is rejected with an error. Sources and coverage are in `lib/nuclear_responses/README.md`.

### DM-nucleus scattering: chiral EFT

`rapidd.chiral_eft` computes the SI and SD rates of a spin-1/2 DM particle with the chiral EFT structure factors, which include the two-body currents: SI from [arXiv:1812.05617](https://arxiv.org/abs/1812.05617), SD from [arXiv:2007.08529](https://arxiv.org/abs/2007.08529) (the nominal structure factors of the LZ SD limits) and [arXiv:2109.00193](https://arxiv.org/abs/2109.00193) (the nuclear uncertainty band, `difrate_chiral_band`). The couplings (`ChiralCouplings`) are built from the Wilson coefficients of the DM-quark operators (`couplings_from_quarks`) or from the MadDM alphas (`alpha_map.chiral_couplings_mdm`). `calc_dRdE.DDrate_save` uses these rates by default (`nuclear_model="chiral_eft"`; `"nreft"` for the one-body $\mathcal O_1$, $\mathcal O_4$ rates).

`equivalent_sigma` (`calc_dRdE.equivalent_sigma_mdm` for MadDM) gives, for each LZ 2024 limit, the cross section of the reference model with the same expected LZ signal as the full model, which can be compared directly with the published curve.

### DM-electron scattering

RAPIDD computes the ionization rate of Xenon and Argon (`target="Xe"` or `"Ar"`) and the excitation rate of silicon (`target="Si"`, see below) for the non-relativistic effective theory of DM-electron interactions of [arXiv:1912.08204](https://arxiv.org/abs/1912.08204) (Catena, Emken, Spaldin, Tarantino): operators $\mathcal{O}_1, \mathcal{O}_3, \dots, \mathcal{O}_{15}$ with contact and long-range couplings, $\mathcal{M} = \sum_i (c_i^s + (q_\text{ref}/q)^2 c_i^\ell) \langle\mathcal{O}_i\rangle$ with $q_\text{ref} = \alpha m_e$, and the atomic responses of `lib/atomic_responses` (computed with [DarkART](https://github.com/temken/DarkART)).

- **$q_\text{max} = 1$ MeV**: the $q$ integral stops at the end of the tables, beyond which the non-relativistic treatment is not valid. For $m_\chi > 0.3$ GeV and a contact coefficient of $\mathcal{O}_3$, $\mathcal{O}_5$, $\mathcal{O}_6$, $\mathcal{O}_{13}$, $\mathcal{O}_{14}$ or $\mathcal{O}_{15}$, `difrate_electron` issues a `QMaxTruncationWarning` with the fraction of the rate that is lost. $\mathcal{O}_1$ and $\mathcal{O}_4$, the only operators used by MadDM, are not affected.
- **Sign of $W_2$**: by default RAPIDD uses the sign corrected in [arXiv:2405.04855](https://arxiv.org/abs/2405.04855); `rapidd.set_w2_sign("Catena")` restores the sign of arXiv:1912.08204.
- **Experiments**: XENON10 and XENON1T (S2-only data release). The XENON1T response is evaluated at $E_e$ plus the energy of the de-excitation photon of the shell; `rapidd.set_er_energy_convention("deposited")` uses the total deposited energy $E_e + E_b$ instead.

```python
import rapidd
rapidd.read_halo('lib/halo_table/halo_table_SHM_1912.08204.dat')  # halo of the paper, with i = 2
rapidd.read_atomic_responses()                                    # W_1 ... W_4

rapidd.reset_coefficients()
rapidd.set_any_Ecoeff(1e-3, 8)                           # c_8^s
rapidd.set_any_Ecoeff(1e-3, 15, interaction="long-range")  # c_15^l
rate = rapidd.difrate_electron(0.4, 0.1, jchi=0.5, shells=rapidd.XE_SHELLS_PAPER)  # E [keV], dRdE, dRdlnE, per shell
only_R2 = rapidd.difrate_electron(0.4, 0.1, resp=2)       # single term R_2 W_2 (can be negative)
total = rapidd.total_rate_electron(0.4, 0.1)              # events/kg/day

rapidd.reset_coefficients()
rapidd.set_any_E_EFT_coeff("magnetic_dipole", 0.1, 1e-5)  # also "anapole", "electric_dipole"
signal = rapidd.electron_signal("Xenon10", rhoDM=0.4, mDM=0.1)
limit = rapidd.electron_limit(signal, coupling=1e-5)     # 90% CL on g/Lambda
```

| Function | Purpose |
|----------|---------|
| `set_any_Ecoeff(c, op, interaction)`, `Ce_val`, `Ce_LR_val`, `reset_coefficients` | contact / long-range coefficients |
| `difrate_electron`, `total_rate_electron` | spectra (total, per shell, per term $R_j W_j$) and total rate; `target="Xe"` (default), `"Ar"` or `"Si"` |
| `QMaxTruncationWarning`, `qmax_lost_fraction`, `QMAX_LOST_FRACTION` | truncation of the $q$ integral at 1 MeV for heavy DM |
| `electron_signal`, `difrate_dS2` | XENON10 / XENON1T S2 spectra and signals |
| `set_dark_photon_coeff`, `set_any_E_EFT_coeff`, `electron_model_coeffs` | dark photon (from $\bar\sigma_e$), anapole, magnetic and electric dipole ($g_e = 2$) |
| `sigma_e_from_coeff`, `coeff_from_sigma_e` | reference cross section $\bar\sigma_e = \mu^2 c^2/(16\pi m_\chi^2 m_e^2)$ |
| `electron_limit`, `electron_limits_vs_mass` | 90% CL limits (XENON10: strongest bin; XENON1T: ROI with background) |
| `set_w2_sign`, `w2_sign` | sign convention of $W_2$ ("Liang", default, or "Catena") |
| `set_er_energy_convention`, `er_energy_convention` | energy of the XENON1T response: $E_e$ + de-excitation photon ("electron", default) or $E_e + E_b$ ("deposited") |
| `dm_response_electron` | DM responses $R_j$ (legacy normalization, used by MadDM for `sigma_e`) |

#### Silicon (crystal target)

`target="Si"` gives the DM-electron rate in silicon with the crystal formalism of [arXiv:2105.02233](https://arxiv.org/abs/2105.02233) (the atomic formula with the unit cell as target and the deposited energy $\Delta E$), using the first crystal response $W_1 = \Delta E\, q^2 V_\text{cell} \,\text{Im}\,\epsilon/(\pi^2\alpha)$ obtained from the RPA dielectric function with local field effects of [QCDark2](https://github.com/meganhott/QCDark2) ([arXiv:2603.12326](https://arxiv.org/abs/2603.12326)); data and sources in `lib/crystal_responses/README.md`.

- **Couplings**: only $c_1, c_4, c_6, c_9, c_{10}, c_{11}$ (contact and long-range), whose DM response is $R_1$ alone; the others need the crystal responses $W_2 \dots W_5$, not available, and raise `UnsupportedCrystalCouplingError`.
- **Screening**, following Catena and Spaldin ([arXiv:2402.06817](https://arxiv.org/abs/2402.06817)): the couplings to the electron density ($c_1$, $c_{11}$) are screened in the RPA, $W_1$ from $\text{Im}[-1/\epsilon]$ with local field effects (as in QCDark2); the spin couplings ($c_4$, $c_6$, $c_9$, $c_{10}$) are not screened in a non-spin-polarised material and use the independent-particle response, $W_1$ from $\text{Im}\,\epsilon$ without local field effects (neglected in arXiv:2402.06817). `screening=False` uses $\text{Im}\,\epsilon$ with local field effects for the density couplings.
- **Range**: $\Delta E \le 50$ eV (no Si 2p), $q \le 25\,\alpha m_e$; `CrystalTruncationWarning` when more than 2% of the rate comes from the edge of the tables.
- **Ionization**: `difrate_dQ` gives the rate per number of electron-hole pairs $Q = 1 \dots 10$ with the probabilities of [Ramanathan and Kurinsky](https://arxiv.org/abs/2004.10709) at 100 K.

```python
rapidd.read_halo(rapidd.halo_path)
rapidd.reset_coefficients()
rapidd.set_any_Ecoeff(rapidd.coeff_from_sigma_e(1e-38, 0.01), 1)    # heavy mediator, sigma_e = 1e-38 cm^2
rate = rapidd.difrate_electron(0.3, 0.01, target="Si")              # Delta E [keV], dRdE [events/kg/day/keV]
per_Q = rapidd.difrate_dQ(rate)                                     # Q, dRdQ [events/kg/day]
```

The QCDark2 heavy- and light-mediator rates are reproduced to $10^{-5}$ with the same quadrature, and within 0.4% (2.6% at 1 MeV) for the total rate and the rate per electron-hole pair with the log-log interpolation in $q$ of RAPIDD (`tests/test_electron_silicon.py`). `examples/damicm_dme_limit_rapidd.py` rescales the DAMIC-M limit (QCDark1 with analytic screening) to this response (`lib/experiments/DAMIC-M/DAMIC-M_2025_DMe_heavy_QCDark_RAPIDD.dat`).

MadDM only uses the SI/SD path ($c_1^s$, $c_4^s$ at $q = 0$) through `rapidd.calc_dRdE.DDrate_electron_save` (targets `Xe`, `Ar`, `Si`; for silicon `dRdE_electron_Si.dat` and `dRdQ_electron_Si.dat`). The comparisons with arXiv:1912.08204 and with the XENON1T and XENONnT data releases are in `examples/checks_ER_rates.ipynb` and in the tests.

## Tests

```bash
pip install pytest
python -m pytest tests
```

`tests/data/extract_paper_figures.py` re-extracts the data of the figures of arXiv:1912.08204 from the PDFs of its source (requires PyMuPDF).

## Project structure

```
RAPIDD_for_DM/
├── instructions_cmake.sh   # convenience script to build the C library
├── lib/                    # C source code, CMake build system, data tables
│   ├── source/              # core C routines (halo, form factors, cross sections, ...)
│   ├── nuclear_responses/   # fit coefficients of the nuclear responses of each nuclear model (DM-nucleus scattering)
│   ├── atomic_responses/    # Xenon and Argon atomic responses W_1 ... W_4 (DM-electron scattering)
│   ├── crystal_responses/   # silicon crystal response W_1 (QCDark2) and ionization probabilities
│   ├── efficiency_tables/   # detector efficiency data
│   ├── experiments/         # experiment-specific data
│   ├── halo_table/          # halo integral tables
│   └── build/               # created after building; contains libRAPIDD and executables
├── rapidd/                 # Python package (interface to the C library)
├── examples/
│   ├── evaluate_NR_rates.ipynb   # nuclear recoil rates (tutorial)
│   ├── checks_NR_rates.ipynb     # nuclear models, modulation, chiral EFT and LZ limits, 19F responses
│   ├── evaluate_ER_rates.ipynb   # DM-electron scattering rates, signals and limits (tutorial)
│   ├── checks_ER_rates.ipynb     # arXiv:1912.08204, XENON1T release, wimprates and XENONnT limits
│   ├── xenonnt_dme_limit_rapidd_ff.py  # XENONnT DM-electron limit rescaled to the RAPIDD atomic responses
│   └── damicm_dme_limit_rapidd.py      # DAMIC-M DM-electron limit rescaled to the RAPIDD silicon response
└── tests/                  # pytest suite
    └── data/                # regression references, data extracted from arXiv:1912.08204
```
