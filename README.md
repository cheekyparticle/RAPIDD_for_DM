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
| `rapidd.electron`     | DM-electron scattering: NREFT ionization rates, S2 signals, models and limits |
| `rapidd.stats`        | Statistical routines for deriving exclusion limits                |

Worked examples are in `examples/evaluate_NR_rates.ipynb` (nuclear recoils) and `examples/evaluate_ER_rates.ipynb` (DM-electron scattering).

### DM-electron scattering

RAPIDD computes the ionization rate of Xenon for the non-relativistic effective theory of DM-electron interactions of [arXiv:1912.08204](https://arxiv.org/abs/1912.08204) (Catena, Emken, Spaldin, Tarantino):

- **Operators**: $\mathcal{O}_1, \mathcal{O}_3, \dots, \mathcal{O}_{15}$ of Tab. I of the paper (DM spin 1/2), with their interferences ($\mathcal{O}_4$-$\mathcal{O}_6$, $\mathcal{O}_{12}$-$\mathcal{O}_{15}$, $\mathcal{O}_{13}$-$\mathcal{O}_{14}$). $\mathcal{O}_2$ and $\mathcal{O}_{16}$ are excluded as in the paper, $\mathcal{O}_{17}, \dots, \mathcal{O}_{20}$ (spin-1 DM) do not enter its DM responses.
- **Couplings**: $\mathcal{M} = \sum_i (c_i^s + (q_\text{ref}/q)^2 c_i^\ell) \langle\mathcal{O}_i\rangle$ with dimensionless contact ($c_i^s$) and long-range ($c_i^\ell$) coefficients and $q_\text{ref} = \alpha m_e$; e.g. $c_1^s = 4 m_\chi m_e a_\text{SI}$ for a contact coupling $a_\text{SI}$ [GeV$^{-2}$]. The DM spin `jchi` is 1/2 for fermions and 0 for scalars (MadDM `dm_spin` $= 2 j_\chi + 1$): the operators with $\mathbf{S}_\chi$ vanish for `jchi = 0`.
- **Rate**: $dR/d\ln E_e = \frac{n_\chi}{128\pi m_\chi^2 m_e^2}\int dq\, q \sum_{j=1}^4 \langle R_j\rangle W_j(k',q)$, with the DM responses $R_1, \dots, R_4$ evaluated at $\mathbf{k} = 0$ (the $\mathbf{k}$ dependence is in the atomic responses) and the tabulated atomic responses $W_1, \dots, W_4$ of the 11 Xenon shells (`lib/atomic_responses`). The terms in $(\mathbf{v}_\text{el}^\perp)^2 = v^2 + \frac{q^2}{4\mu^2}\frac{m_\chi - m_e}{m_\chi + m_e} - \frac{\Delta E}{\mu}$ are integrated with $\int d^3v\, f/v$ (constant part) and $\int d^3v\, f\, v$ ($v^2$ part): the halo table must be written with `i=2` for the operators depending on $\mathbf{v}_\text{el}^\perp$ (an error is raised otherwise).
- **Sign of $W_2$**: by default RAPIDD uses the sign of the vector atomic form factor corrected in [arXiv:2405.04855](https://arxiv.org/abs/2405.04855) (Liang et al.; $R_2$ is the same in the two works). `rapidd.set_w2_sign("Catena")` restores the sign of arXiv:1912.08204, e.g. to reproduce its figures. Only $\mathcal{O}_7, \mathcal{O}_8, \mathcal{O}_{12}, \mathcal{O}_{13}, \mathcal{O}_{14}$ are affected (for $\mathcal{O}_7$ at 10 MeV the rate decreases by a factor ~5); the SI/SD path used by MadDM is not.
- **Experiments**: XENON10 (7 S2 bins, 15 kg day) and XENON1T (S2-only data release response and mass dependent ROI, $20\,\text{MeV} \le m_\chi \le 5\,\text{GeV}$).

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
| `difrate_electron`, `total_rate_electron` | spectra (total, per shell, per term $R_j W_j$) and total rate |
| `electron_signal`, `difrate_dS2` | XENON10 / XENON1T S2 spectra and signals |
| `set_dark_photon_coeff`, `set_any_E_EFT_coeff`, `electron_model_coeffs` | dark photon (from $\bar\sigma_e$), anapole, magnetic and electric dipole ($g_e = 2$) |
| `sigma_e_from_coeff`, `coeff_from_sigma_e` | reference cross section $\bar\sigma_e = \mu^2 c^2/(16\pi m_\chi^2 m_e^2)$ |
| `electron_limit`, `electron_limits_vs_mass` | 90% CL limits (XENON10: strongest bin; XENON1T: ROI with background) |
| `set_w2_sign`, `w2_sign` | sign convention of $W_2$ ("Liang", default, or "Catena") |
| `dm_response_electron` | DM responses $R_j$ (legacy normalization, used by MadDM for `sigma_e`) |

MadDM only uses the SI/SD path ($c_1^s$, $c_4^s$ at $q = 0$) through `rapidd.calc_dRdE.DDrate_electron_save`. The comparison with the figures of the paper is in `tests/test_electron_paper.py` (with the sign of $W_2$ of the paper): the XENON10 limits and the model spectra agree within $\sim 25\%$, while the spectra of Figs. 5-7 and the rates of Fig. 6 of the paper are lower by a kinematics-dependent factor 1.1-1.45 (arXiv:2405.04855 also finds Fig. 7 to be ~3/4 of the formula; and its $\mathcal{O}_3$ bars do not follow its own DM responses); these are marked as expected failures.

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
│   ├── atomic_responses/    # Xenon atomic responses W_1 ... W_4 (DM-electron scattering)
│   ├── efficiency_tables/   # detector efficiency data
│   ├── experiments/         # experiment-specific data
│   ├── halo_table/          # halo integral tables
│   └── build/               # created after building; contains libRAPIDD and executables
├── rapidd/                 # Python package (interface to the C library)
├── examples/
│   ├── evaluate_NR_rates.ipynb   # nuclear recoil rates
│   └── evaluate_ER_rates.ipynb   # DM-electron scattering rates, signals and limits
└── tests/                  # pytest suite
    └── data/                # regression references, data extracted from arXiv:1912.08204
```
