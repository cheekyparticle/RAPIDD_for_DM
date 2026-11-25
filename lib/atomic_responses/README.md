# Xenon and Argon atomic responses W_1 ... W_4

`atomic_response_<j>/<El>_<shell>.txt`: atomic response W_j(k', q) of arXiv:1912.08204 for the 11 Xenon shells and
the 5 Argon shells (1s, 2s, 2p, 3s, 3p). Each file has 3 comment lines (`//`), then 100 rows of k' (log steps in
[0.1, 100] keV) times 200 columns of q (log steps in [1, 1000] keV).

The tables are the output (`<shell>_<j>_Table.txt`) of [DarkART](https://github.com/temken/DarkART), commit 594ffb5,
the C++ successor of DarkARC used by the paper, run in `Tabulation` mode with

    element = "Xe"; atomic_responses = (1, 2, 3, 4);
    k_min = 0.1; k_max = 100.0; k_points = 100;
    q_min = 1.0; q_max = 1000.0; q_points = 200;
    tabulate_radial_functions = true; r_points = 5000;

and `element = "Ar"; atomic_shells = ("1s", "2s", "2p", "3s", "3p");` for Argon. The Argon binding energies of RAPIDD
(`rapidd.AR_BINDING_ENERGIES`) are those of the DarkART data files `data/Ar_<shell>.txt`.

`Xe_4d.txt` was regenerated with DarkART 594ffb5: the previous table predated the DarkART fix of the 4d binding
energy (commit d58015b), which entered the effective charge of the continuum states. The other Xenon shells are
unchanged.

The tables stop at q = 1 MeV ~ 2 m_e, beyond which the non-relativistic treatment is not valid (see the
`QMaxTruncationWarning` of `rapidd.electron`). The q_max ~ 1.05 MeV of the paper quoted in `README.md` was found with
DarkART tables extended beyond 1 MeV, which are not shipped.
