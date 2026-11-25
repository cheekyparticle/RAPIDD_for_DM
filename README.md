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

RAPIDD computes the ionization rate of Xenon and Argon (`target="Xe"` or `"Ar"`) for the non-relativistic effective theory of DM-electron interactions of [arXiv:1912.08204](https://arxiv.org/abs/1912.08204) (Catena, Emken, Spaldin, Tarantino): operators $\mathcal{O}_1, \mathcal{O}_3, \dots, \mathcal{O}_{15}$ with contact and long-range couplings, $\mathcal{M} = \sum_i (c_i^s + (q_\text{ref}/q)^2 c_i^\ell) \langle\mathcal{O}_i\rangle$ with $q_\text{ref} = \alpha m_e$, and the atomic responses of `lib/atomic_responses` (computed with [DarkART](https://github.com/temken/DarkART)).

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
| `difrate_electron`, `total_rate_electron` | spectra (total, per shell, per term $R_j W_j$) and total rate; `target="Xe"` (default) or `"Ar"` |
| `QMaxTruncationWarning`, `qmax_lost_fraction`, `QMAX_LOST_FRACTION` | truncation of the $q$ integral at 1 MeV for heavy DM |
| `electron_signal`, `difrate_dS2` | XENON10 / XENON1T S2 spectra and signals |
| `set_dark_photon_coeff`, `set_any_E_EFT_coeff`, `electron_model_coeffs` | dark photon (from $\bar\sigma_e$), anapole, magnetic and electric dipole ($g_e = 2$) |
| `sigma_e_from_coeff`, `coeff_from_sigma_e` | reference cross section $\bar\sigma_e = \mu^2 c^2/(16\pi m_\chi^2 m_e^2)$ |
| `electron_limit`, `electron_limits_vs_mass` | 90% CL limits (XENON10: strongest bin; XENON1T: ROI with background) |
| `set_w2_sign`, `w2_sign` | sign convention of $W_2$ ("Liang", default, or "Catena") |
| `set_er_energy_convention`, `er_energy_convention` | energy of the XENON1T response: $E_e$ + de-excitation photon ("electron", default) or $E_e + E_b$ ("deposited") |
| `dm_response_electron` | DM responses $R_j$ (legacy normalization, used by MadDM for `sigma_e`) |

MadDM only uses the SI/SD path ($c_1^s$, $c_4^s$ at $q = 0$) through `rapidd.calc_dRdE.DDrate_electron_save`. The comparisons with arXiv:1912.08204 and with the XENON1T and XENONnT data releases are in `examples/checks_ER_rates.ipynb` and in the tests.

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
│   └── xenonnt_dme_limit_rapidd_ff.py  # XENONnT DM-electron limit rescaled to the RAPIDD atomic responses
└── tests/                  # pytest suite
    └── data/                # regression references, data extracted from arXiv:1912.08204
```
