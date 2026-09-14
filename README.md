# RAPIDD

RAPIDD is a tool for computing dark matter direct detection rates and exclusion limits. It combines a C library for the numerical computation of differential recoil rates with a Python interface for setting up dark matter models, generating rates, and deriving experimental limits.

This branch (`maddm`) integrates RAPIDD with [MadDM](https://github.com/maddmhep/maddm), allowing direct detection rates to be computed directly from a MadDM model.

## Authors

- D. G. Cerdeño — davidg.cerdeno@gmail.com
- A. Cheek — andrew.cheek@durham.ac.uk
- E. Reid — elliott.m.reid@durham.ac.uk
- H. Schulz — iamholger@googlemail.com

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

Clone the repository and switch to the `maddm` branch:

```bash
git clone https://github.com/cheekyparticle/RAPIDD_for_DM.git
cd RAPIDD_for_DM
git checkout maddm
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
| `rapidd.stats`        | Statistical routines for deriving exclusion limits                |

## Project structure

```
RAPIDD_for_DM/
├── instructions_cmake.sh   # convenience script to build the C library
├── lib/                    # C source code, CMake build system, data tables
│   ├── source/              # core C routines (halo, form factors, cross sections, ...)
│   ├── efficiency_tables/   # detector efficiency data
│   ├── experiments/         # experiment-specific data
│   └── build/               # created after building; contains libRAPIDD and executables
└── rapidd/                 # Python package (interface to the C library)
```
