# XENONnT limits (verbatim copies of the data releases)

Used by `examples/evaluate_ER_rates.ipynb`, `examples/checks_ER_rates.ipynb` and `examples/checks_NR_rates.ipynb`.

| file | search | source |
|:--|:--|:--|
| `XENONnT_SR012_s2_only_wimp_electron_qedark_after_pcl_limit.csv` | S2-only, DM-electron, heavy mediator (QEDark form factors); same as `lib/experiments/XENONnT/XENONnT_DMe_SI.dat` used by MadDM | [s2_only_data_release](https://github.com/XENONnT/s2_only_data_release), arXiv:2601.11296 |
| `XENONnT_SR012_s2_only_wimp_er_qedark_after_pcl_limit.csv` | S2-only, DM-electron, heavy mediator (QEDark form factors), not shown in the paper: probably the signal with the total deposited energy $E_e + E_b$ (the `electron` file counts the primary electron and the de-excitation photons, as XENON1T): 3.9 times stronger than it at 40 MeV, 1.2 times weaker at 10 GeV. Not used | same |
| `XENONnT_SR012_s2_only_wimp_si_after_pcl_limit.csv` | S2-only, SI DM-nucleon | same |
| `XENONnT_SR012_s2_only_wimp_sdn_after_pcl_limit.csv` | S2-only, SD DM-neutron | same |
| `XENONnT_SR012_s2_only_wimp_si_light_mediator_after_pcl_limit.csv` | S2-only, SI DM-nucleon, light mediator | same |
| `xenonnt_2025_si_wimp.csv`, `xenonnt_2025_sd_neutron_wimp.csv`, `xenonnt_2025_sd_proton_wimp.csv` | WIMP search SR0+SR1 | [wimp_data_release](https://github.com/XENONnT/wimp_data_release), arXiv:2502.18005 |

Limits in cm^2, 90% CL, after the power constraint.
