#ifndef BLABLA
#define BLABLA
/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
Physical constants
##########################################*/

#pragma once

const double static mproton = 0.93827208943; // proton mass in GeV (CODATA 2022)
const double static mneutron = 0.93956542194; // neutron mass in GeV (CODATA 2022)
/* Nucleon mass of the NREFT operators (q/m_N) and of the matching of the relativistic Lagrangians
(Table 1 of arXiv:1308.6288): (m_p + m_n)/2. The same m_N is used for protons and neutrons: it is the
scale that makes the operators dimensionless, so the m_N of a coefficient must be the one of its operator. */
const double static mnucleon = 0.938918755685; // GeV
const double static atomic_mass_unit_GeV = 0.93149410372; // u in GeV (CODATA 2022)
const double static melectron = 0.51099895069e-3; // electron mass in GeV (CODATA 2022)
const double static c = 2.99792458e+5; //Speed of light in km/s
const double static G_F = 1.1663788e-5; //Fermi coupling constant in GeV^-2 (PDG 2024)
const double static higgs_vev = 246.2196402393; //GeV, (sqrt(2) G_F)^(-1/2)
const double static e = 0.3; // e charge
const double static s_w = 0.4895916665957459;
const double static g_n = -3.83;
const double static g_p = 5.59;

// DM-electron scattering (exact double-precision values, as in the MadDM Fortran)
const double static c_light_cms = 29979245800.; // speed of light in cm/s
const double static GeV_to_kg = 1.78266192e-27; // GeV/c^2 -> kg
const double static GeV_to_J = 1.602176634e-10; // GeV -> J
const double static atomic_mass_unit_kg = 1.66053906660e-27; // u -> kg
const double static m_Xe_avg_au = 131.293; // Xenon mass, average over all isotopes (u)
const double static m_Ar_avg_au = 39.948; // Argon mass, standard atomic weight (u), DM-electron rate
const double static gev_to_pb = 3.893793721e+8; // (hbar c)^2: GeV^-2 -> pb
const double static pb_to_cm2 = 1.e-36; // pb -> cm^2
const double static m_to_cm = 100.; // m -> cm
const double static day_to_sec = 86400.; // day -> s
const double static alpha_em = 7.2973525643e-3; // fine-structure constant (CODATA 2022), q_ref = alpha_em m_e

/* DM-nucleus targets (name, isotopes, abundances): table nuclear_targets in phys_consts.c */
#define MAX_TARGET_ISOTOPES 10
int n_nuclear_targets(void);
const char * nuclear_target_name(int i);
const char * nuclear_target_file_name(const char *target);
int nuclear_target_isotopes(const char *target, int *A, int *Z, double *weight);
int nuclear_target_isotope(int i, int k, int *A, int *Z, double *abundance);
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double reduced_mass(double mass_1, double mass_2);
double nucleus_mass(int A, int Z);
void isotope_mass_fractions(const int *A, int Z, const double *abundance, int n, double *mass_fraction);

double sigma_to_f_11(double sigma, double mchi, double mnucleon);
double sigma_to_a_44(double sigma, double mchi, double mnucleon);


int countlines(char *filename);

#endif
