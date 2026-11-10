# Xenon and Argon atomic responses W_1 ... W_4

`atomic_response_<j>/Xe_<shell>.txt`: atomic response W_j(k', q) of arXiv:1912.08204 for the 11 Xenon
shells, 100 rows of k' (log steps in [0.1, 100] keV) times 200 columns of q (log steps in [1, 1000] keV),
after 3 comment lines (`//`). `atomic_response_<j>/Ar_<shell>.txt`: the same for the 5 Argon shells
(1s, 2s, 2p, 3s, 3p).

DarkART is the C++ successor of DarkARC (Python, v1.0 of 2019, DOI 10.5281/zenodo.3581334, no longer
maintained), the code used for arXiv:1912.08204 (q in [1, 1000] keV on a 100 x 100 grid, l' <= 7 fixed,
adaptive radial integrals with tolerances 1e-3/1e-2). The DarkART 594ffb5 tables are converged in
r_points (5000 vs 10000: identical) and in l' (l' <= 40 fixed: W_j within 0.5%, transverse part
(q/m_e)^2 W_3 - W_4 within 0.25%): where they differ from DarkARC, the DarkART/RAPIDD results are the
reference. The q range ends at q_max = 1000 keV ~ 2 m_e, beyond which the non-relativistic treatment is
not valid (see the QMaxTruncationWarning of rapidd.electron).

The tables are the output (`<shell>_<j>_Table.txt`) of [DarkART](https://github.com/temken/DarkART)
(Roothaan-Hartree-Fock bound states, hydrogenic continuum states with Z_eff from the binding energy,
sum over l' until convergence), run in `Tabulation` mode with

    element = "Xe"; atomic_responses = (1, 2, 3, 4);
    k_min = 0.1; k_max = 100.0; k_points = 100;
    q_min = 1.0; q_max = 1000.0; q_points = 200;
    tabulate_radial_functions = true; r_points = 5000;

- `Xe_4d.txt` (all four responses) was regenerated on 2026-11-09 with DarkART 594ffb5 (master). The
  previous tables had been produced before the DarkART commit d58015b ("Fix binding energy of
  Xe-4d", Dec. 2022): `data/Xe_4d.txt` contained the binding energy of 5p (0.457 au instead of
  2.778 au), which enters the effective charge of the final (continuum) states, Z_eff = 3.83
  instead of 9.43. The old W_1 of 4d was up to ~4 times smaller at q >= 10 keV; DarkART 594ffb5 with
  the old binding energy reproduces the old tables exactly.
- The other shells are unchanged: DarkART 594ffb5 reproduces them exactly (checked on 10 rows of k' for 5p).

- `Ar_<shell>.txt` (all four responses, 1s, 2s, 2p, 3s, 3p) were generated on 2026-11-09 with DarkART
  594ffb5 (master), the settings above with `element = "Ar"; atomic_shells = ("1s", "2s", "2p", "3s",
  "3p");` (`<shell>_<j>_Table.txt` renamed `atomic_response_<j>/Ar_<shell>.txt`). They coincide exactly
  with an independent run on a k' grid of 34 points (every third row) and with a run for
  q in [1000, 1110] keV at q = 1 MeV. The binding energies used by RAPIDD (`rapidd.AR_BINDING_ENERGIES`)
  are those of the DarkART data files `data/Ar_<shell>.txt` (Roothaan-Hartree-Fock, Bunge et al. 1993),
  times the Hartree energy 27.211386245988 eV: 3227.55 (1s), 335.30 (2s), 260.45 (2p), 34.76 (3s),
  16.08 (3p) eV. The Argon rates reproduce the Argon bars of Fig. 6 and the spectra of Figs. 5 and 7
  of arXiv:1912.08204 within 1% at 10 and 100 MeV and for the long-range couplings
  (tests/test_electron_paper.py).

The effective q_max ~ 1.05 MeV of the paper (explanation of the 1 GeV contact residual of O_3, O_6,
O_13, O_14 in tests/test_electron_paper.py, README.md and examples/checks_ER_rates.ipynb) was obtained
with DarkART 594ffb5 tables extended beyond 1 MeV, which are **not shipped** (q > 1 MeV is outside the
validity of the non-relativistic treatment and RAPIDD never uses it). The numbers quoted there (Argon
at 1 GeV with q <= 1.045 MeV: O_6 1.000, O_13/O_14 0.991, O_3+ 0.981 of the paper; Xenon, 8 shells:
within 0.3%) are written as text and are not recomputed by the tests or the notebook. To reproduce
them, use DarkART 594ffb5 in `Tabulation` mode with two local patches of `src/Atomic_Responses.cpp`
and `src/Radial_Integrator.cpp` (environment variables, no change of the default behaviour):
`DARKART_NO_ONE_CORR=1` skips the negative "one correction" of W_1 (the q = 0 overlap term, irrelevant
at q >= 1 MeV: the extended Argon run coincides exactly with the shipped tables at q = 1 MeV), and
`DARKART_RMAX_A0` sets the end of the radial grid in Bohr radii (default 50). The two runs were:

- Xenon, 8 shells (5p, 5s, 4d, 4p, 4s, 3d, 3p, 3s), `DARKART_RMAX_A0=15 DARKART_NO_ONE_CORR=1`, the
  k' grid of the shipped tables split in two runs (`ext_lo`: k_min = 0.1, k_max = k'_49, 50 points;
  `ext_hi`: k_min = k'_50, k_max = 100.0, 50 points; k'_i = 10^(-1 + 3 i / 99) keV) and

      q_min = 1000.0; q_max = 10^(3 + 48*3/199) = 5291.97...; q_points = 49; r_points = 8000;

  i.e. the q step of the shipped tables continued up to ~5.3 MeV;
- Argon, 5 shells, `DARKART_NO_ONE_CORR=1`, k_points = 34 (every third row of the shipped grid) and

      q_min = 1000.0; q_max = 10^(3 + 3*3/199) = 1109.75...; q_points = 4; r_points = 5000;

The rate is then integrated in ln q from 1 keV up to q = 1.045-1.05 MeV.

The 4d W_1 was also checked against an independent evaluation of Eq. (W_1) of arXiv:1912.08204
(median difference 0.15%), and the 5p W_1 against Fig. 2 of the paper (within 1%).
