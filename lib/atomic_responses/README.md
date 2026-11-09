# Xenon atomic responses W_1 ... W_4

`atomic_response_<j>/Xe_<shell>.txt`: atomic response W_j(k', q) of arXiv:1912.08204 for the 11 Xenon
shells, 100 rows of k' (log steps in [0.1, 100] keV) times 200 columns of q (log steps in [1, 1000] keV),
after 3 comment lines (`//`).

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

The 4d W_1 was also checked against an independent evaluation of Eq. (W_1) of arXiv:1912.08204
(median difference 0.15%), and the 5p W_1 against Fig. 2 of the paper (within 1%).
