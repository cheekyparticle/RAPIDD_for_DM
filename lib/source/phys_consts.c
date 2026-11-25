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
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include <string.h>

#include "phys_consts.h"


/////Frequently used equations///////////

double reduced_mass(double mass_1, double mass_2){

	return  mass_1*mass_2 / (mass_1 + mass_2);

}

/* Atomic masses in u: AME2020 (Wang et al., Chinese Phys. C 45, 030003 (2021)). The older values of 16O, 19F,
40Ar, 40Ca, Ge and Xe differ from AME2020 by up to 1.6e-5 u (1.3e-7 relative) and are kept unchanged. */
static const struct { int A; int Z; double mass_u; } atomic_masses[] = {
	{16, 8, 15.994914630}, {19, 9, 18.998403224}, {40, 18, 39.962383681}, {40, 20, 39.962590588},
	{12, 6, 12.}, {13, 6, 13.003354835}, {18, 8, 17.999159612}, {23, 11, 22.989769282},
	{27, 13, 26.981538408}, {28, 14, 27.976926535}, {29, 14, 28.976494665}, {30, 14, 29.973770136},
	{127, 53, 126.904472592}, {133, 55, 132.905451961},
	{70, 32, 69.924249656}, {72, 32, 71.922078859}, {73, 32, 72.923462608}, {74, 32, 73.921177394},
	{76, 32, 75.921401565},
	{124, 54, 123.905894239}, {126, 54, 125.904280981}, {128, 54, 127.903531222}, {129, 54, 128.904780103},
	{130, 54, 129.903509426}, {131, 54, 130.905072012}, {132, 54, 131.904144154}, {134, 54, 133.905394991},
	{136, 54, 135.907213743},
	{0, 0, 0.}
};

double nucleus_mass(int A, int Z){
	/*Mass of the nucleus (A, Z) in GeV: atomic mass (AME2020) minus the mass of the Z electrons (their
	binding energy, < 1e-6 of the mass, is neglected). It is the m_T of q^2 = 2 m_T E_R, of the DM-nucleus
	reduced mass and of vmin. The sum of the masses of the free nucleons, Z m_p + (A - Z) m_n, is ~0.9%
	larger (nuclear binding energy).
	For an isotope that is not tabulated the atomic mass is approximated by A u (error < 0.1%).
	Inputs: A (mass number), Z (atomic number).*/
	double mass_u = (double)A;
	for (int i = 0; atomic_masses[i].A != 0; i++) {
		if (atomic_masses[i].A == A && atomic_masses[i].Z == Z) {
			mass_u = atomic_masses[i].mass_u;
			break;
		}
	}
	return mass_u * atomic_mass_unit_GeV - Z * melectron;
}

/* DM-nucleus targets of difrate_dER and nuclear_model_supports_target (one line per target).
Abundances: mole (number) fractions of the CIAAW (ciaaw.org): "Isotopic compositions of the elements 2013" (Meija
et al., Pure Appl. Chem. 88, 293 (2016)) for Xe, Ge, Ca and Ar; for C, O and Si, for which it gives only intervals,
the CIAAW representative compositions of 1997 (C, O) and 2001 (Si), which lie inside them.
For an element the rate per kg of target weights each isotope with its mass fraction (isotope_mass_fractions);
the isotopes that are not listed (no nuclear responses) do not contribute and their mass fraction is not given to
the others:
- C: 12C, 13C (all the stable isotopes).
- O: 17O (0.00038) is not included (it is not in the Elastic repository).
- Si: 28Si, 29Si, 30Si (all the stable isotopes).
- Ar: 40Ar only, with abundance 1 (36Ar 0.003336 and 38Ar 0.000629 are neglected and given to 40Ar).
- Ca: 42Ca (0.00647), 43Ca (0.00135), 44Ca (0.02086), 46Ca (0.00004) and 48Ca (0.00187) are not included
  (42Ca and 44Ca are in the Elastic repository only with other pf-shell interactions than 40Ca).
- Ge: 70-76Ge (all the stable isotopes).
- Xe: 124Xe (0.00095) and 126Xe (0.00089) are not included.
- F, Na, Al, I, Cs: one stable isotope.
CaWO4 (compound, element = 0): the weights are used as given (16O and 40Ca; W has no nuclear form factors). */
static const struct nuclear_target {
	const char *name;       /* target of difrate_dER */
	const char *file_name;  /* name of the output and experiment files of WriteData.c */
	int element;            /* 1: weights = mass fractions (isotope_mass_fractions), 0: abundance as given */
	int n;
	int Z[MAX_TARGET_ISOTOPES];
	int A[MAX_TARGET_ISOTOPES];
	double abundance[MAX_TARGET_ISOTOPES];
} nuclear_targets[] = {
	{"Xe", "XENON", 1, 7, {54, 54, 54, 54, 54, 54, 54}, {128, 129, 130, 131, 132, 134, 136},
		{0.01910, 0.26401, 0.04071, 0.21232, 0.26909, 0.10436, 0.08857}},
	{"Ge", "GERMANIUM", 1, 5, {32, 32, 32, 32, 32}, {70, 72, 73, 74, 76}, {0.2052, 0.2745, 0.0776, 0.3652, 0.0775}},
	{"Ar", "ARGON", 1, 1, {18}, {40}, {1.0}},
	{"F", "FLUORINE", 1, 1, {9}, {19}, {1.0}},
	{"C", "CARBON", 1, 2, {6, 6}, {12, 13}, {0.9893, 0.0107}},
	{"O", "OXYGEN", 1, 2, {8, 8}, {16, 18}, {0.99757, 0.00205}},
	{"Na", "SODIUM", 1, 1, {11}, {23}, {1.0}},
	{"Al", "ALUMINIUM", 1, 1, {13}, {27}, {1.0}},
	{"Si", "SILICON", 1, 3, {14, 14, 14}, {28, 29, 30}, {0.92223, 0.04685, 0.03092}},
	{"Ca", "CALCIUM", 1, 1, {20}, {40}, {0.96941}},
	{"I", "IODINE", 1, 1, {53}, {127}, {1.0}},
	{"Cs", "CAESIUM", 1, 1, {55}, {133}, {1.0}},
	{"CaWO4", "CaWO4", 0, 2, {8, 20}, {16, 40}, {0.616, 0.0}},   /* {0.616, 0.384} with 40Ca */
};
#define N_NUCLEAR_TARGETS ((int)(sizeof(nuclear_targets) / sizeof(nuclear_targets[0])))

static const struct nuclear_target * find_nuclear_target(const char *target){
	for (int t = 0; t < N_NUCLEAR_TARGETS; t++) {
		if (strcmp(target, nuclear_targets[t].name) == 0) return &nuclear_targets[t];
	}
	return NULL;
}

int n_nuclear_targets(void){
	/*Number of DM-nucleus targets.*/
	return N_NUCLEAR_TARGETS;
}

const char * nuclear_target_name(int i){
	/*Name of the i-th DM-nucleus target (NULL if i is out of range).*/
	return (i >= 0 && i < N_NUCLEAR_TARGETS) ? nuclear_targets[i].name : NULL;
}

const char * nuclear_target_file_name(const char *target){
	/*Name of the output and experiment files of the target (e.g. "XENON"), NULL for an unknown target.*/
	const struct nuclear_target *t = find_nuclear_target(target);
	return (t != NULL) ? t->file_name : NULL;
}

int nuclear_target_isotopes(const char *target, int *A, int *Z, double *weight){
	/*Isotopes of the target and their weights in the rate per kg of target: the mass fractions for an
	element, the abundances as given for a compound (CaWO4).
	Outputs: A, Z, weight (arrays of at least MAX_TARGET_ISOTOPES elements).
	Returns the number of isotopes, 0 for an unknown target.*/
	const struct nuclear_target *t = find_nuclear_target(target);
	if (t == NULL) return 0;
	for (int l = 0; l < t->n; l++) {
		A[l] = t->A[l];
		Z[l] = t->Z[l];
		weight[l] = t->abundance[l];
	}
	if (t->element) isotope_mass_fractions(t->A, t->Z[0], t->abundance, t->n, weight);
	return t->n;
}

int nuclear_target_isotope(int i, int k, int *A, int *Z, double *abundance){
	/*k-th isotope of the i-th target: A, Z and abundance (mole fraction). Returns 1, 0 if out of range.*/
	if (i < 0 || i >= N_NUCLEAR_TARGETS || k < 0 || k >= nuclear_targets[i].n) return 0;
	*A = nuclear_targets[i].A[k];
	*Z = nuclear_targets[i].Z[k];
	*abundance = nuclear_targets[i].abundance[k];
	return 1;
}

void isotope_mass_fractions(const int *A, int Z, const double *abundance, int n, double *mass_fraction){
	/*Mass fractions of the isotopes of an element from their mole (number) fractions:
	w_i = f_i m_i / m_avg, with m_avg = sum_j f_j m_j / sum_j f_j over the n listed isotopes.
	The differential rate of each isotope is per kg of that isotope, so the rate per kg of the element
	is sum_i w_i R_i. sum_i w_i = sum_i f_i: isotopes that are not listed (no form factors) do not
	contribute to the rate, but their mass fraction is not given to the others.
	Inputs: A[n] (mass numbers), Z, abundance[n] (mole fractions). Output: mass_fraction[n].*/
	double sum_f = 0., sum_fm = 0.;
	for (int l = 0; l < n; l++) {
		sum_f += abundance[l];
		sum_fm += abundance[l] * nucleus_mass(A[l], Z);
	}
	for (int l = 0; l < n; l++) {
		mass_fraction[l] = (sum_fm > 0.) ? abundance[l] * nucleus_mass(A[l], Z) * sum_f / sum_fm : 0.;
	}
}

////////Function to calculate couplings f (f_n, f_p) for SI interaction by using the cross section sigma//////////////////////////
double sigma_to_f_11(double sigma, double mchi, double mnucleon){

	double mun = reduced_mass(mchi, mnucleon); //reduced mass DM particle and nucleon (mproton or mneutron)
	return sqrt(sigma*M_PI) / (2.*mun);

}

////////Function to calculate couplings a (a_n, a_p) for SD interaction by using the cross section sigma//////////////////////////
double sigma_to_a_44(double sigma, double mchi, double mnucleon){

	double mun = reduced_mass(mchi, mnucleon); //reduced mass DM particle and nucleon (mproton or mneutron)
	return sqrt(sigma*M_PI) / (sqrt(24.)*mun*G_F);

}


int countlines(char *filename){
  // count the number of lines in the file called filename
  FILE *fp = fopen(filename,"r");
  int ch=0;
  int lines=0;

  if (fp == NULL){
  return 0;
  }


  while(!feof(fp))
  {
    ch = fgetc(fp);
    if(ch == '\n')
    {
      lines++;
    }
  }

  fclose(fp);
	//printf("attempt %i\n", lines );
  return lines;
}
